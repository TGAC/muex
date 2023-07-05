import pysam as ps
import pickle

def extract_cigar(infile,min=6,max=21):
	"""
	Identifies reads with an insertion in the range [min,max] in their CIGAR strings, and extracts the QNAME, RNAME, CIGAR and POS of these reads
	Input:
	infile (string): filepath to alignment file
	min (int): minimum insert size
	max (int): maximum insert size
	Output: 
	reads (list): list of tuples with relevant read info
	"""
	samfile = ps.AlignmentFile(infile)
	reads = []
	for read in samfile:
		if not read.is_unmapped:
			cigar = read.cigartuples
			for t in cigar:
				# Insertions in the range [min,max]
				if t[0] == 1 and t[1] >= min and t[1] <= max:
					reads.append((read.query_name, read.reference_name, read.cigartuples, read.reference_start))
					break
	return reads

def find_insertions(alignments_path: str, transcript_dict_path: str, out_reads_to_extract_path: str, out_reads_with_inserts_path: str, min_size:int=6, max_size:int=21):
	"""
	Finds insertions in alignment cigar strings that exit between two annotated exons and are between min_size and max_size nt long.
	:param alignments_path: Filepath to alignments
	:param transcript_dict_path: Filepath to transcript dict
	:param out_reads_to_extract_path: Filepath to output reads list
	:param out_reads_with_inserts_path: Filepath to output summary file
	:param min_size: min insertion size (nt)
	:param max_size: max insertion size (nt)
	"""
	print("ingesting and parsing alignments")
	reads = extract_cigar(alignments_path,min=min_size,max=max_size)
	print("finding transcriptomic position of each insertion of interest")
	reads_of_interest = []
	unique_reads = set()
	with open(transcript_dict_path, 'rb') as f:
			transcript_dict = pickle.load(f)
	for read in reads:
		transcript = read[1]
		# This is a hacky workaround. depends on .fa file's sequence headers. could tell users to format their files with just ENST... (or else warn them that this program might be sensitive to this)
		if not transcript.startswith("ENST") or "|" in transcript:
			rname_sections = transcript.split("|")
			for s in rname_sections:
				if s.startswith("ENST") and "." in s: # NB if not "." in s, it picks the last one it sees, which doesn't contain the '.x' suffix. yeah, I did say this was hacky didn't I?
					transcript = s
		# end hacky workaround
		transcript_dict_entry = transcript_dict[transcript]
		transcript_positions = transcript_dict_entry[2]
		strand = transcript_dict_entry[1]
		gene_id, gene_name = transcript_dict_entry[3], transcript_dict_entry[4]
		
		query_pos_counter = read[3] + 1 # +1 to change from 0-ind (pySAM) to 1-ind (SAM) 
		ref_pos_counter = read[3] + 1

		operation_positions = []
		# Assumption: first operation is a Match
		# Walks the CIGAR string, noting start position of each op in ref and query
		for t in read[2]:
			operation_positions.append((t,ref_pos_counter,query_pos_counter))
			if t[0] == 0: # match
				ref_pos_counter += t[1]
				query_pos_counter += t[1]
			elif t[0] == 1: # ins
				query_pos_counter += t[1]
			elif t[0] == 2: # del
				ref_pos_counter += t[1]

		# Keeps only inserts that occur on the boundary between two exons in the transcript
		# Consequence: microexons/insertions at the start or end of a transcript will not be reported
		for o in operation_positions:
			if o[0][0] == 1 and o[0][1] >= min_size and o[0][1] <= max_size:
				insert_pos = o[1] - 1 # -1 to get last base of preceding op
				# Assumption: list of exon positions from dictionary is ordered and consecutive
				for idx,exon in enumerate(transcript_positions):
					if insert_pos > exon[3]:
						# If insert_pos is larger than the transcript end of an exon, keep going
						continue
					elif insert_pos == exon[3]:
						# If insert_pos matches the transcript end of an exon, it's at an exon junction
						# Read Name, Transcript ID, Strand, Insert Length, Pos of Insert in Reference, Pos of Insert in Read, Genomic End of Exon before Insert, Genomic Start of Exon After Insert, Gene ID, Gene Name
						if strand == "+":
							reads_of_interest.append((read[0],transcript,strand,o[0][1],o[1],o[2],exon[1],transcript_positions[idx+1][0],gene_id,gene_name))
							unique_reads.add(read[0])
						elif strand == "-":
							reads_of_interest.append((read[0],transcript,strand,o[0][1],o[1],o[2],exon[0],transcript_positions[idx-1][1],gene_id,gene_name))
							unique_reads.add(read[0])
						else:
							pass
					else: 
						# Therefore if insert_pos gets this far, it can't be between any other exons
						break
	
	print("generating output: list of reads to extract")
	reads_file = []
	for entry in list(unique_reads):
		reads_file.append(f"{entry}\n")
	print("writing output: list of reads to extract")
	with open(out_reads_to_extract_path,"wt") as out_fh:
		out_fh.writelines(reads_file)

	print("generating output: reads with inserts")
	tsv_file = []
	for entry in reads_of_interest:
		tsv_file.append(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\t{entry[4]}\t{entry[5]}\t{entry[6]}\t{entry[7]}\t{entry[8]}\t{entry[9]}\n")
	print("writing output: reads with inserts")
	with open(out_reads_with_inserts_path,"wt") as out_fh:
		out_fh.writelines(tsv_file)

if __name__ == "__main__":
	find_insertions(snakemake.input["alignment"],snakemake.input["transcript_dict"],snakemake.output["readnames"],snakemake.output["reads_with_inserts"],min_size=snakemake.params["min_size"],max_size=snakemake.params["max_size"])
