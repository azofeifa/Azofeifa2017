'''
	Requires:
		1) ../files/HCT116_EMG.bed
		2) ../files/TSS.bed
	Outputs:
		1) STable_to_1_promoter_distances.csv
		Table Fields:
			chromosome of bidir. txn.,origin of bidir. txn. (hg19),nearest promoter,distance
'''
import scipy.spatial as ss
import numpy as np,math
def load_bed(FILE,KD=False):
	FH,G 	=open(FILE), {}
	IDS={}
	for line in FH:
		if line[0]!="#":
			line_array 				= line.strip("\n").split("\t")
			chrom,start, stop 	= line_array[:3]
			if chrom not in G:
				G[chrom] 			= list()
				IDS[chrom] 			= list()
			IDS[chrom].append(line_array[-1])
			G[chrom].append((float(start) + float(stop))/2.0)
	for chrom in G:
		x 				= np.zeros((len(G[chrom]),2))
		x[:,0] 		 	= G[chrom]
		if KD:
			G[chrom] 	= ss.KDTree(x)
		else:
			G[chrom] 	= x
	return G, IDS
def get_nearest(A,B,IDS,OUT=""):
	D 		= list()
	FHW  	= open(OUT, "w")
	FHW.write("chromosome of bidir. txn.,origin of bidir. txn. (hg19),nearest promoter,distance\n")
	for chrom in A:
		if chrom in B:
			a,b 	= A[chrom], B[chrom]
			for i in range(a.shape[0]):
				d,idx 	=  b.query(a[i,:])
				if d > 0:
					D.append(math.log(d,10))
				FHW.write(chrom+ "," + str(int(a[i,0]))  + "," + str(IDS[chrom][idx]) + "," + str(int(d))  + "\n" )

def main():
	FILE1	= "../files/HCT116_EMG.bed"
	FILE2 = "../files/TSS.bed"
	OUT 	= "../Tables/STable_to_1_promoter_distances.csv"
	COMP = False
	(A,IDS),(B,IDS) 	= load_bed(FILE1), load_bed(FILE2,KD=True)
	get_nearest(A,B,IDS,OUT=OUT)


if __name__ == "__main__":
	main()