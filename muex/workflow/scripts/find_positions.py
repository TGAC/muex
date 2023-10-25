import pysam as ps
import csv

def gen_strands_dict(reads_with_inserts_path):
	strands = {}
	with open(reads_with_inserts_path,"rt") as csv_file:
		f = csv.reader(csv_file,delimiter="\t")
		for row in f:
			if row[0] in strands:
				pass
			else:
				# strand, transcript ID, gene name
				strands[row[0]] = (row[2],row[1],row[-1])
	return strands

def find_genomic_positions(alignment_path, out_bed_path, out_supp_info_path, strands, upstream_flank=2, downstream_flank=2, min_size=6, max_size=21, leeway=0):
	"""
	# TODO docstring
	"""
	print("ingesting and parsing alignments")
	samfile = ps.AlignmentFile(alignment_path)
	
	print("finding genomic position of each insertion of interest")
	microexons = {}
	for read in samfile:
		if not read.is_unmapped:
			readname = read.query_name
			chr = read.reference_name
			cigar = read.cigartuples
			query_pos_counter = read.reference_start + 1
			ref_pos_counter = read.reference_start + 1

			operation_positions = []
			# Walks the CIGAR string, noting start position of each op in ref and query
			for t in cigar:
				operation_positions.append((t,ref_pos_counter,query_pos_counter))
				if t[0] == 0 or t[0] == 7 or t[0] == 8: # match (M), seq match (=), seq mismatch (X)
					ref_pos_counter += t[1]
					query_pos_counter += t[1]
				elif t[0] == 1 or t[0] == 4: # ins (I), softclip (S) 
					query_pos_counter += t[1]
				elif t[0] == 2 or t[0] == 3: # del (D), skip (N)
					ref_pos_counter += t[1]

			# Extracts the start and end positions of sequence of operations that fall between two introns (N) and that include an expected amount of matches ([min_size-leeway..max_size+leeway], where leeway accounts for long-read noise)
			skip_next_n = False # ensures N operations treated as pairs, avoids idx out of range error
			sequences = []
			for idx,op in enumerate(operation_positions):
				if skip_next_n:
					skip_next_n = False
					continue
				if op[0][0] == 3: # TODO space for sanity check: which operations either side of N?
					next_op_idx = idx+1
					next_op = operation_positions[next_op_idx]
					seq_start = next_op[1]
					inner_matches = 0
					# Scans ahead until it encounters another N, or end of list
					while next_op[0][0] != 3 and inner_matches <= 23:
						if next_op[0][0] == 0:
							# Keeps track of how many Match ops it has encountered, to determine whether to keep the sequence
							inner_matches += next_op[0][1]
						next_op_idx += 1
						if next_op_idx >= len(operation_positions):
							break
						next_op = operation_positions[next_op_idx]
					if inner_matches >= min_size-leeway and inner_matches <= max_size+leeway and next_op[0][0] == 3:
						seq_end = next_op[1]-1
						sequences.append((seq_start,seq_end))
						skip_next_n = True

			for s in sequences:
				# 1-based inclusive genomic position of microexon
				microexon_start = s[0] 
				microexon_end = s[1]

				strand, transc_ID, gene_name = strands[readname]

				muex_id = "_".join(str(x) for x in (chr,strand,microexon_start,microexon_end))
				
				if muex_id in microexons:
					microexons[muex_id][-1].append(readname)
					microexons[muex_id][-2] += 1
				else:
					microexons[muex_id] = [chr,strand,microexon_start,microexon_end,microexon_end-microexon_start+1,gene_name,transc_ID,1,[readname]]
		
	print("generating output: positions .bed file")
	bedfile = []
	for muex_id,info in microexons.items():
		bedfile.append(f"{info[0]}\t{int(info[2])-int(upstream_flank)-1}\t{int(info[3])+int(downstream_flank)}\t{muex_id}\t0\t{info[1]}\n")
		# start pos altered bc BED pos is half-open 0-indexed eg [0,100) = the first 100 bases = chr1:1-100 = bases 0 to 99
	print("writing output")
	with open(out_bed_path,"wt") as out_fh:
		out_fh.writelines(bedfile)

	print("generating output: microexon info")
	info_file = []
	for muex_id,info in microexons.items():
		info_string = "\t".join([str(x) for x in info[:-1]]) # TODO can amend to output comma-sep list of readnames, if useful
		full_string = f"{muex_id}\t{info_string}\n"
		info_file.append(full_string)
	with open(out_supp_info_path,"wt") as out_fh:
		out_fh.writelines(info_file)

if __name__ == "__main__":
	print("generating readname:strand dictionary")
	strands = gen_strands_dict(snakemake.input["reads_with_inserts"])
	print("finding genomic positions")
	find_genomic_positions(snakemake.input["alignment"],snakemake.output["positions"],snakemake.output["supporting_info"],strands,upstream_flank=snakemake.params["upflank"],downstream_flank=snakemake.params["downflank"],leeway=snakemake.params["leeway"])