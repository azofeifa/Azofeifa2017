'''
	Requires:
		1) ../files/K562_EMG.bed ; bed file of bidirectional centers in K562
		2) ../files/ChIP_files/K562/; directory of bed files
		3) ../files/STable_to_1C_eRNA_TF_mark_ref.csv
	Outputs:
		1) STable_to_1D_ChIP_BIC_displacements.tsv
'''

import os
import matplotlib.pyplot as plt
import numpy as np
import NormalLaplace as NL
import math,pandas as pd
def load_bed(FILE, W=0):
	G,FH 	= {},open(FILE, "r")
	for line in FH:
		if line[0]!="#":
			chrom,start, stop 	= line.split("\t")[:3]
			if chrom not in G:
				G[chrom] 			= list()
			x 	= 	0.5*(int(start) + int(stop))
			if not W:
				G[chrom].append(x)
			else:
				G[chrom].append([x-W, x+W, list()])
	for chrom in G:
		G[chrom].sort()
	return G
def get_eRNAs(FILE, W=1500):
	df 	= pd.read_csv(FILE)
	df 	= df[df["TSS"]==0]
	G 		= {}
	for j,i in df.iterrows():
		chrom,start, stop = i.chrom,i.start,i.stop
		x 						= 0.5*(start+stop) 
		if chrom not in G:
			G[chrom] 		= list()
		G[chrom].append((int(x)-W, int(x) + W))
	for chrom in G:
		G[chrom].sort()
	return G

def load_meta(FILE):
	M,FH 	= {},open(FILE)
	for line in FH:
		line_array 	= line.split("\t")
		M[line_array[0]] 	= line_array[15].strip("-human")
	return M
def load_peaks(DIR):
	M,T 			= load_meta(DIR+"metadata.tsv"),{}
	for i,FILE in enumerate(os.listdir(DIR)):
		n 	= FILE.split(".")[0]
		if "bed" in FILE and n in M and "eGFP" not in M[n]  :
			G 			= load_bed(DIR+FILE)
			T[FILE] 	= G
	return T,M
def compute_displacements(B,T,M,OUT="",BINS=100):
	FHW 	= open(OUT, "w")
	D 		= {}
	FHW.write("ENCODE_ID,name,pk_displacement,dBIC,"  + ",".join([str(int(i)) for i in np.linspace(-1500,1500,BINS)]) + "\n"  )

	for TF in T:
		D[TF] 	= list()
		name 				= M[TF.split(".")[0]]
		for chrom in T[TF]:
			if chrom in B:
				a,b 	=  B[chrom],T[TF][chrom]
				j,N 	= 0,len(b)
				for i in range(len(a)):
					while j < N and b[j] < a[i][0]:
						j+=1
					k 	= j
					while k < N and b[k] < a[i][1]:
						d 		= b[k] - (0.5*(a[i][1] + a[i][0]))
						D[TF].append(d)
						k+=1
		counts,edges 	= np.histogram(D[TF],bins=100,range=(-1500,1500))
		XX 				= np.zeros((100,2))
		XX[:,0] 			= (edges[1:] + edges[:-1])/2.
		XX[:,1] 			= counts

		EM 				= NL.EM(k=2,T=200, ct=pow(10,-6))
		EM2 				= NL.EM(k=1,T=200, ct=pow(10,-6))
		EM3 				= NL.EM(k=0,T=200, ct=pow(10,-6))


		N 					= np.sum(XX[:,1])
		rvs,ll   		= EM.fit(XX,2)
		rvs2,ll2   		= EM2.fit(XX,1)
		rvs3,ll3   		= EM3.fit(XX,0)
		B1,B2,B3 		= -2*ll + 8*math.log(N),-2*ll2 + 4*math.log(N),-2*ll3 + 1*math.log(N)
		


		MU 				= sum([abs(rv.mu) for rv in rvs if rv.type=="laplace" ])*0.5
		W 					= sum([abs(rv.w)  for rv in rvs if rv.type=="laplace" ])*0.5
		FHW.write(TF.strip(".bed") + "," + name + "," + str(MU) + "," + str(B2-B1) + "," )
		FHW.write(",".join(map(str, counts)) + "\n") 



if __name__ == "__main__":
	RERUN 	= True
	OUT 		= "../Tables/STable_to_1D_ChIP_BIC_displacements.csv"
	K562 		= "../Tables/STable_to_1C_eRNA_TF_mark_ref.csv"
	DIR 		= "../files/ChIP_files/K562/"
	B 			= get_eRNAs(K562)
	T,M 		= load_peaks(DIR)

	compute_displacements(B,T,M,OUT=OUT)
