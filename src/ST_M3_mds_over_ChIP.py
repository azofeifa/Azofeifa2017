'''
	Requires:
		1) A directory of motif displacement files from MDS program;
			each MDS file relates to the Core 2014 K562 GRO-cap dataset. 
			Moreso, each MDS file was ran with -TSS option for a different
			set of TF ChIP-seq peaks
	Outputs:
		A Table: STable_to_3_mds_over_ChIP.csv
		Fields : "ENCODE(ID)","name","MotifModel" ,"lack_hist", "over_hist"
'''

import os,numpy as np,pandas as pd
def extract_from_file(FILE,B=150,ID=""):
	FH 			= open(FILE,"r")
	G,collect 	= {},False
	for line in FH:
		if "Empiracle" in line:
			break
		elif "Binned" in line:
			collect 		= True
		elif collect:
			motif,x,y,z 	= line.strip("\n").split("\t")
			if ID in motif:
				x,y 				= map(float, x.split(",")),map(float, y.split(","))
				xs 				= range(len(x))
				G[motif] 	= [np.histogram(xs,bins=150,normed=1,weights=x)[0],np.histogram(xs,bins=150,normed=1,weights=y)[0]]
	return G


def load_meta_file(FILE):
	FH 		= open(FILE)
	header 	= True
	M 			= {}
	for line in FH:
		if not header:
			line_array 			= line.split("\t")
			M[line_array[0]] 	= line_array[15].split("-")[0]
		else:
			header= False
	return M
def get_max_difference(G):
	MAX 	= -np.inf
	ARGMAX=None
	for motif in G:
		diff 	= np.sum(np.abs(G[motif][0]-G[motif][1]))
		if diff > MAX:
			MAX=diff;ARGMAX=motif
	if ARGMAX:
		return G[ARGMAX][0],G[ARGMAX][1],ARGMAX
	return None,None,None


def make_table(DIR,OUT):
	df 			= pd.DataFrame(columns=("ENCODE(ID)","name","MotifModel" ,"lack_hist", "over_hist") )
	M 				= load_meta_file(DIR+"metadata.tsv")
	i 				= 0
	for f in os.listdir(DIR):
		if f.split("_")[0] in M:
			if "eGFP" not in M[f.split("_")[0]]:
				G 			 	= extract_from_file(DIR+f, ID=M[f.split("_")[0]])
				if G:
					non,over,motif 	= get_max_difference(G)
					if non is not None and over is not None:
						non,over 	= map(str, non), map(str, over)
						df.loc[i] 	= f.split("_")[0], M[f.split("_")[0]],motif, "|".join(non), "|".join(over)
						i+=1
	df.to_csv(OUT,index=False)
def main():
	DIR 	= "/Users/joazofeifa/Lab/new_motif_distances/over_ChIP_non_promoter/"
	OUT 	= "../Tables/STable_to_3_mds_over_ChIP.csv"
	make_table(DIR, OUT)

if __name__ == "__main__":
	main()








