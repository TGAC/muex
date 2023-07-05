class Muex:
		def __init__(self,id,chr,strand,start_pos,end_pos,gene=None,transcript_id=None,seq=None,canon_ss=None,frame=None,supp_reads=None):
			self.id = id
			self.chr = chr
			self.strand = strand
			self.start_pos = start_pos
			self.end_pos = end_pos
			self.length = self.end_pos-self.start_pos+1
			self.gene = gene
			self.transcript_id = transcript_id
			self.seq = seq
			self.canon_ss = canon_ss
			self.frame = frame
			self.supp_reads = supp_reads

		def __str__(self):
			return f"{self.chr} {self.strand} {self.start_pos} {self.end_pos}"