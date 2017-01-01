import os,numpy as np
'''
	Requires:
		1) Directory of MD score files from MDS program
	Outputs:
		1) STable_to_3D_displacements_motif_by_SRAnumber.csv

			-> each row is a motif / SRR combination MD score, 
				N and displacement histogram (~150 bins) precomputed
'''
def get_MDS(counts,percent=0.1):
	Length, Center 	= len(counts),len(counts)/2
	N 						= sum(counts)+1
	d 						= Length*(percent/2)
	SN 					= sum(counts[Center-d:Center+d])
	return N, float(SN)/float(N)
def get_counts(FILE,M,B=150):
	collect 	= False
	FH 			= open(FILE, "r")
	for line in FH:
		if  "Emp" in line and collect:
			break
		elif "Binned" in line:
			collect 	= True
		elif collect:
			line_array 	= line.strip("\n").split("\t")
			x 				= map(float,line_array[1].split(","))
			counts,x1 	= np.histogram(np.linspace(-1500,1500,len(x)), weights=x,bins=B)
			motif 		= line_array[0]
			if motif not in M:
				M[motif] = list()
			M[motif].append(counts)
	return M
def load_dir_for_counts(DIR):
	bins 	= 150
	labels= [str(int(i)) for i in np.linspace(-1500,1500,bins) ]
	S,M 	= list(),dict()
	for i,FILE in enumerate(os.listdir(DIR)):
		M 	= get_counts(DIR+FILE,M,B=bins)
		S.append(FILE.split("_")[0])
	return S, M,labels
def write_out(S,M,labels, OUT):
	f 				= lambda x: str(int(x))
	percent 		= 0.1
	FHW 			= open(OUT, "w")
	FHW.write("MotifModel,SRAnumber,N,MDS(" + str(percent)+ ")," + ",".join(labels) + "\n")
	for m in M:
		for i in range(len(M[m])):
			N,MDS 	= get_MDS(M[m][i])
			FHW.write(m+"," + S[i] + "," + str(N) + "," + str(MDS) + ",".join(map(f, M[m][i])) + "\n")
	FHW.close()
	pass
def main():
	DIR 			= "/Users/joazofeifa/Lab/new_motif_distances/human/300/"
	S,M,labels	= load_dir_for_counts(DIR)
	write_out(S,M,labels, "../Tables/STable_to_3D_displacements_motif_by_SRAnumber.csv")

if __name__ == "__main__":
	main()



