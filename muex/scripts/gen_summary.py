from class_muex import Muex
import csv

def inspect_splice_sites(muex_fasta_path,upstream_flank=2,downstream_flank=2):
	"""
	Input:
	muex_fasta (filepath): path to fasta file containing microexon sequences
	upflank (int): how many nt upstream of the exon to inspect
	downflank (int): how many nt downstream of the exon to inspect
	"""
	# This function currently expects only 2nt either side of the sequence, eg the first and last two nucleotides of any given sequence will be considered splice sites
	# TODO could check if supplied sequence is >=5nt (2+1+2), for cases where users might want shorter muexs than our defaults

	with open(muex_fasta_path, "r") as fa:
		muex_seqs = {}
		lines = fa.readlines()
		# Assumption: All lines (incl. final) end in \n (final depends on behaviour of `bedtools getfasta`). TODO can be made more robust by stripping whitespace?
		for i, line in enumerate(lines):
			if line[0] == ">":
				# extracts info from seqname header line in .fa
				muex_id = line[1:-4]
				strand = muex_id.split("_")[-3]
				if strand == "-":
					# Swaps to be consistent with behaviour of bedtools getfasta
					upstream_flank,downstream_flank = downstream_flank,upstream_flank
				
				# extracts sequences
				seq_flanked = lines[i+1][:-1] # [:-1] strips \n
				upflank_seq = seq_flanked[:upstream_flank] # gets upstream flank
				downflank_seq = seq_flanked[-downstream_flank:] # gets downstream flank
				seq = seq_flanked[upstream_flank:-downstream_flank] # gets just the sequence  

				if upflank_seq[-2:] == "AG" and downflank_seq[:2] == "GT":
					muex_seqs[muex_id] = (seq,True)
				else:
					muex_seqs[muex_id] = (seq,False)

		# return canon, noncanon
		return muex_seqs

if __name__ == "__main__":
	# id, chr, strand, start pos, end pos, length, supporting reads, gene name, transc ID from microexon_supporting_reads.tsv
	# seq, canon_ss from running splice site analysis function on sequences.fa
	all_muexs = []
	muex_seqs = inspect_splice_sites(snakemake.input["seqs"],upstream_flank=snakemake.params["upflank"],downstream_flank=snakemake.params["downflank"])
	with open(snakemake.input["supporting_info"],"rt") as csv_file:
		f = csv.reader(csv_file,delimiter="\t")
		for row in f:
			muex_id = row[0]
			muex_to_add = Muex(id=muex_id, chr=row[1], strand=row[2], start_pos=int(row[3]), end_pos=int(row[4]), gene=row[6], transcript_id=row[7], seq=muex_seqs[muex_id][0], canon_ss=muex_seqs[muex_id][1], supp_reads=row[8])
			all_muexs.append(muex_to_add)
	
	with open(snakemake.output["summary"],"wt") as f:
		file_lines = []
		for muex in all_muexs:
			file_lines.append(f"{muex.id}\t{muex.chr}\t{muex.strand}\t{muex.start_pos}\t{muex.end_pos}\t{muex.length}\t{muex.gene}\t{muex.transcript_id}\t{muex.seq}\t{muex.canon_ss}\t{muex.supp_reads}\n") # TODO make this cleaner by moving the tabulated print format to __str__/__repr__ attribute of Muex class, or something similar
		f.writelines(file_lines)