import numpy as np
import os
'''
	Requires:
		1) K562_EMG.bed (or some bed file of bidirs)
	Outputs
		1) STable_to_1B_eRNA_ref_marks_overlap.csv

'''
def load_bed_file(FILE):
	G,FH 	= {},open(FILE,'r')
	for line in FH:
		if line[0]!="#":
			chrom,start,stop 	= line.split("\t")[:3]
			if chrom not in G:
				G[chrom] 		= list()
			G[chrom].append((int(start), int(stop)))
	for chrom in G:
		G[chrom].sort()
	return G
def load_bed_file(FILE,w=None):
	G,FH 	= {},open(FILE,'r')
	for line in FH:
		if line[0]!="#":
			chrom,start,stop 	= line.split("\t")[:3]
			if chrom not in G:
				G[chrom] 		= list()
			if w is None:
				G[chrom].append((int(start), int(stop)))
			else:
				x 		= int((int(stop) + int(start))/2.)
				G[chrom].append((x-w, x+w))

	for chrom in G:
		G[chrom].sort()
	return G
def load_meta(FILE):
	FH 		= open(FILE)
	header 	= True
	M 			= {}
	for line in FH:
		if not header:
			line_array 			= line.split("\t")
			if "eGFP" not in line_array[15]:
				M[line_array[0]] 	= line_array[15].split("-")[0]
		else:
			header 		= False
	return M
def get_marks():
	DIR 		= "../files/marks/K562/"
	FILES 	= [x for x in os.listdir(DIR) if "meta" not in x]
	META 		= load_meta(DIR+"metadata.tsv")

	IDS2 		= [META[f.split(".")[0]] for f in FILES ]
	IDS 		= IDS2
	MARKS 	= [ load_bed_file(DIR+f) for f in FILES]
	return MARKS, IDS
def overlap2(A,B,OVER={}):
	for chrom in A:
		if chrom not in OVER:
			OVER[chrom] 	= [[x[0], x[1], list()] for x in A[chrom]]
		for i in range(len(OVER[chrom])):
			OVER[chrom][i][2].append(0)
		if chrom in B:
			a,b 	 	= A[chrom], B[chrom]

			j,N 		= 0, len(b)
			for i,x  in enumerate(a):
				while j < N and b[j][1] < x[0]:
					j+=1
				if j < N and b[j][0] < x[1]:
					OVER[chrom][i][2][-1]= 1
	return OVER
def get_counts(BIDR, MARKS, IDS,OUT):
	FHW 	= open(OUT, "w")
	FHW.write("chrom,start,stop,RefSeq(TSS)," + ",".join(IDS) + "\n")
	R 		= load_bed_file("../files/TSS.bed",w=100)
	B 		= overlap2(BIDR, R)

	ALL 	= float(sum([len(B[chrom]) for chrom in B]))
	for i,M in enumerate(MARKS):
		B 	= overlap2(B, M,OVER=B)
	for chrom in B:
		for start, stop, D in B[chrom]: 
			FHW.write( chrom + "," + str(start) + "," + str(stop) + ","  + ",".join(map(str, D)).lstrip(",") + "\n"  )

def main():
	BIDIR 		= "../files/K562_EMG.bed"
	OUT 			= "../Tables/STable_to_1B_eRNA_ref_marks_overlap.csv"

	BIDR 			= load_bed_file(BIDIR)
	MARKS,IDS  	= get_marks()

	get_counts(BIDR, MARKS,IDS,OUT)
if __name__ == "__main__":
	main()

