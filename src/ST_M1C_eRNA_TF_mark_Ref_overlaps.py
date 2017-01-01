'''
	This requires:
		1) ../files/ChIP_files/K562/ 
		2) ../files/marks/K562/Histone/ 
		3) "../files/marks/K562/DHS1/ 
		3) ../files/K562_EMG.bed
		4) ../files/TSS.bed
	Outputs:
		1) ../Tables/STable_to_1C_eRNA_TF_mark_ref.csv
	
		for each TF the counts association with a mark
		conditioned on eRNA presence

		2) ../Tables/STable_to_1C_Used_ENCODE_TFs.csv
		3) ../Tables/STable_to_1C_Used_ENCODE_Marks.csv
'''

import os
import matplotlib.pyplot as plt
import numpy as np
from pyliftover import LiftOver
import seaborn as sns,time

def load_bed_file(FILE,AA=None,PAD=750,LO=False):
	IO 		= 	None
	if LO:
		print "...lifting over..."
		IO 	= LiftOver('hg38', 'hg19')

	G,FH 	= {},open(FILE,'r')
	for line in FH:
		if line[0]!="#":
			chrom,start,stop 	= line.split("\t")[:3]
			x 						= 0.5*(float(start)+float(stop))
			if IO is not None:
				info 				= IO.convert_coordinate(chrom, int(x), '-')
				if info is not None and len(info):
					chrom,x 		= info[0][:2]
			if chrom not in G:
				G[chrom] 		= list()
			if AA is None:
				G[chrom].append((int(x-PAD), int(x+PAD)))
			else:
				G[chrom].append([ x-PAD, x+PAD, [0 for x in range(AA)] ])
	for chrom in G:
		G[chrom].sort()
	return G
def overlap(A,B,u,OVER=False, LACK=False):
	for chrom in A:
		if chrom in B:
			a,b 	 	= A[chrom], B[chrom]
			j,N 		= 0, len(b)
			for i,(start, stop,v) in enumerate(a):
				while j < N and b[j][1] < start:
					j+=1
				if j < N and b[j][0] < stop:
					A[chrom][i][2][u]=1
	return A
def load_meta(FILE,I=15):
	FH 		= open(FILE)
	header 	= True
	M 			= {}
	V 			= {}
	AS 		= 0
	for line in FH:
		if not header:
			line_array 			= line.split("\t")
			if "eGFP" not in line_array[I]:
				M[line_array[0]] 	= line_array[I].split("-")[0]
				V[line_array[0]] 	= line_array[AS].split("-")[0]
		else:
			AS 	= [i for i,x in enumerate(line.strip("\n").split("\t")) if x=="Assembly"][0]
			print AS
			header= False
	return M,V
def get_bed_files(DIR,OUT=None, I=15):
	FILES 	= [x for x in os.listdir(DIR)  if "meta" not in x]
	META,V 	= load_meta(DIR+"metadata.tsv",I=I)
	FILES 	= [f for f in FILES if f.split(".")[0] in META and "eGFP" not in META[f.split(".")[0]]]
	IDS 		= [META[f.split(".")[0]] for f in FILES ]
	MARKS 	= [load_bed_file(DIR+f,LO=(V[f.split(".")[0]]!="hg19") ) for f in FILES  ]
	if OUT is not None:
		FHW 		= open(OUT, "w")
		FHW.write("ENCODE_ID,name,peak_total\n")
		for i,f in enumerate(FILES):
			N 		= sum([len(MARKS[i][chrom]) for chrom in MARKS[i] ])
			FHW.write(f.split(".")[0] + "," + IDS[i] + "," + str(N) + "\n" )
		FHW.close()
	return MARKS, IDS

def get_overlap_vectors(MARKS, IDS,B, OUT, SUM=False):
	FHW 		= open(OUT, "w")
	FHW.write('chrom,start,stop,' + ",".join(IDS) + "\n" )
	FHW.flush()

	for i,M in enumerate(MARKS):
		overlap(B,M,u=i)
	for chrom in B:
		for start, stop, v in B[chrom]:
			FHW.write(chrom+"," + str(start) + "," + str(stop) + ",")
			FHW.write(",".join(map(str, v) ) + "\n")
def main():
	TF_ChIP_DIR 		= "../files/ChIP_files/K562/" 		#directory of bed files
	Mark_ChIP_DIR 		= "../files/marks/K562/Histone/" 	#directory of bed files
	DHS1_DIR 			= "../files/marks/K562/DHS1/" 		#directory of bed files
	BIDIR 				= "../files/K562_EMG.bed"
	Ref 					= "../files/TSS.bed"


	OUT1 					= "../Tables/STable_to_1C_Used_ENCODE_DHS1.csv"
	OUT2 					= "../Tables/STable_to_1C_Used_ENCODE_Marks.csv"
	OUT3 					= "../Tables/STable_to_1C_Used_ENCODE_TFs.csv"
	OUT4 					= "../Tables/STable_to_1C_eRNA_TF_mark_ref.csv"
	
	R 						= load_bed_file(Ref)


	MARKS, IDS2 		= get_bed_files( Mark_ChIP_DIR, OUT=OUT2)
	print "loaded marks"

	DHS1, IDS1 			= get_bed_files( DHS1_DIR, OUT=OUT1,I=4)
	print "loaded DHS1"

	TFS,   IDS3 		= get_bed_files( TF_ChIP_DIR, OUT=OUT3)
	print "loaded TFS"
	B 						= load_bed_file(BIDIR, AA=len(MARKS+TFS+DHS1)+1)
	
	get_overlap_vectors([R] + DHS1+ MARKS+TFS, ["TSS"] + IDS1+ IDS2+IDS3,B, OUT4)
if __name__ == "__main__":
	main()

