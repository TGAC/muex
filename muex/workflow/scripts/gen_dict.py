from gtfparse import read_gtf
import pickle
# import defopt

def gen_dict(infile):
	"""
	Input: 
	infile (string): filepath to annotation file (.gtf)
	Output:
	transcript_exons (dict): dict of relevant info {"Transcript_ID":(chr,strand,[(g_start,g_end,t_start,t_end),...],gene_id,gene_name)...}
	""" 
	# TODO consider adding length of each exon to this dict
	# Ingests GTF input as Pandas df
	df = read_gtf(infile,result_type="pandas") # edit for gtfparse version 2.0.1
	df = df[df["feature"] == "exon"] # keeps only exons

	transcript_exons = {}
	# Generates dictionary of relevant values
	# Assumption: first exon encountered in the annotation file is the first exon in the transcript (regardless of its position in the genome), and the rest are in order.
	for _,row in df.iterrows():
		if row["transcript_id"] in transcript_exons:
			# In transcriptome coordinates, current exon starts one base after previous exon ends
			t_start = transcript_exons[row["transcript_id"]][2][-1][3] + 1
			exon_length = row["end"] - row["start"] + 1
			t_end = t_start + exon_length - 1
			transcript_exons[row["transcript_id"]][2].append((row["start"],row["end"],t_start,t_end))
		else:
			# Assumption: strand is same for all exons in transcript.
			t_start = 1
			t_end = row["end"] - (row["start"]-1)
			transcript_exons[row["transcript_id"]] = (row["seqname"],row["strand"],[(row["start"],row["end"],t_start,t_end)],row["gene_id"],row["gene_name"])
	return transcript_exons

def main(in_gtf_path: str, out_filepath: str):
	"""
	Generates and saves a pickled transcript dictionary from an annotation file.
	:param in_gtf_path: Filepath to input annotation .gtf file
	:param out_filepath: Output filepath for pickled dictionary
	"""
	print("Generating transcript:position dictionary")
	transcript_dict = gen_dict(in_gtf_path)
	print("Pickling and saving dictionary")
	with open(out_filepath,"wb") as out_fh:
		pickle.dump(transcript_dict,out_fh,protocol=pickle.HIGHEST_PROTOCOL)
	
if __name__ == "__main__":
	main(snakemake.input["annotation"],snakemake.output["transcript_dict"])
	# defopt.run(main)