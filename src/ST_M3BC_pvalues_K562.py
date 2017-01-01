'''
	Requires
		1) An MD score file for the Core 2014 K562 GRO-cap experiment
			(../files/SRR1552480_MDS.tsv)
	Outputs:
		1) A table which lists for each motif the p-value of enrichment
			relative to either stationary or non-stationary null model

			../Tables/STable_to_3BC_pvalues_MDS_K562.csv

'''
import math,numpy as np
import scipy.stats as ss
import pandas as pd
def load_MDS(FILE):
	FH 				= open(FILE, "r")
	NS,D,collect 	= {},{},False
	for line in FH:
		if "Empiracle" in line:
			break
		elif "Binned" in line:
			collect=True
		elif collect:
			line_array 			= line.strip("\n").split("\t")
			D[line_array[0]] 	= map(float, line_array[1].split(","))
		elif line[0]!="#":
			line_array 	= line.strip("\n").split("\t")
			motif 		= line_array[0] 
			EX 			= float(line_array[-2].split(",")[0])
			NS[motif] 	= math.exp(float(line_array[4].split(",")[0])),EX
	return NS, D
def get_pvalue(obs,mu, std):
	NN 	= ss.norm(mu, std)
	return NN.cdf(obs)
def get_binomial_pvalue(d, percent=0.1):
	N 	= sum(d)
	c 	= len(d)/2
	w 	= int(len(d)*(percent/2.))
	S 	= sum(d[c-w:c+w])
	EX = N*percent
	STD= math.sqrt(EX*(1-percent))
	return S/N,N,get_pvalue(S, EX,STD)


def main():
	FILE 					= "../files/SRR1552480_MDS.tsv"
	NS_pvals,D 			= load_MDS(FILE)
	bin_params 			= dict([(motif, get_binomial_pvalue(D[motif])) for motif in D])
	df 					= pd.DataFrame(columns=("MotifModel", "MDS" ,"N" , "Stationary p-value", 
															"Non Stationary p-value", "Null Model non-stationary EX"))
	for i,motif in enumerate(bin_params.keys()):
		MDS, N, stat 	= bin_params[motif]
		stat_non, EX 	= NS_pvals[motif]
		df.loc[i] 		= motif, MDS, N, stat,stat_non, EX
	df.to_csv("../Tables/STable_to_3BC_pvalues_MDS_K562.csv")


if __name__ == "__main__":
	main()

