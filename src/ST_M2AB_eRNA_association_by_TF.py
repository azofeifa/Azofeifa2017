'''
	Requires
		1) A directory of MD score files ran from the MDS
			program where each file is named from the
			ENCODE identifier and -TSS option was the ENCODE
			peak file
	Outputs
		a table that for each ENCODE peak file, the percent
		of those peaks overlapped by an eRNA
'''
import pandas as pd,os
def extract_from_file(FILE):
	FH 			= open(FILE,"r")
	info 			= FH.readlines()[16].split("\t")[-1]
	return map(float, info.split(","))
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
def make_table(DIR,OUT):
	df 			= pd.DataFrame(columns=("ENCODE(ID)","name", "total_peaks", "percent_eRNA") )
	M 				= load_meta_file(DIR+"metadata.tsv")
	i 				= 0
	for f in os.listdir(DIR):
		if f.split("_")[0] in M:
			if "eGFP" not in M[f.split("_")[0]]:
				percent, N 	= extract_from_file(DIR+f)
				df.loc[i] 	= f.split("_")[0], M[f.split("_")[0]], N, percent
				i+=1
	df.to_csv(OUT,index=False)
def main():
	DIR 	= "/Users/joazofeifa/Lab/new_motif_distances/over_ChIP/"
	OUT 	= "../Tables/STable_to_2A_eRNA_association_by_TF.csv"
	make_table(DIR, OUT)

if __name__ == "__main__":
	main()





