'''
	Requires:
		1) A directory of md scores labeled by the SRRnumber
			output from the MD score program
	Outputs:
		1) A Table where each row is a motif and model;
			each column is a SRR datasets and each cell is the
			md score 
'''
import os,pandas as pd
def load_MDS(FILE,M):
	FH 	= open(FILE, "r")
	for line in FH:
		if "Binned" in line:
			break
		elif line[0]!="#":
			line_array 	= line.strip("\n").split("\t")
			motif,N,MDS = line_array[0],float(line_array[1].split(",")[0]),float(line_array[2].split(",")[0])		
			if motif not in M:
				M[motif] = list()
			M[motif].append((float(N), float(MDS)))
def load_directory(DIR):
	SRRS,M	= list(),dict()
	for f in os.listdir(DIR):
		load_MDS(DIR+f, M)
		SRRS.append(f.split("_")[0])
	return M,SRRS
def main():
	M,SRRS= load_directory("/Users/joazofeifa/Lab/new_motif_distances/human/300/")
	df 	= pd.DataFrame(columns=["MotifModel"]+ [i + "(MDS)" for i in SRRS]+ [i + "(N)" for i in SRRS] )
	for i,motif in enumerate(M.keys()):
		df.loc[i] 	= [motif] + [mds for n,mds in M[motif]] + [n for n,mds in M[motif]]
	df.to_csv("../Tables/STable_to_3E_MDS_N_across_all.csv",index=False)

if __name__ == "__main__":
	main()


