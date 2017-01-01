import sys
sys.path.append("/Users/joazofeifa/Lab/motif_displacement_analysis/MDAP/")
import differential_MD as dMD
'''
	Requires:
		1) A Directory of MD score computations to be labeled by SRR 
			with the _MDS suffix
	Outputs
		1) an svg figure showing the "MA Plot" of the four examples

'''

def main():
	DIR     = "/Users/joazofeifa/Lab/new_motif_distances/human/300/"
	def add(b,a=DIR , c="_MDS.tsv"):
		return a+b+c

	MDS1    	= dMD.mds_frame()
	MDS1.load_MD_score_file(map(add,["SRR1105736","SRR1105737"]), "DMSO")
	MDS1.load_MD_score_file(map(add,["SRR1105738","SRR1105739"]), "Nutlin")
	df      	= MDS1.differential_single("DMSO", "Nutlin")
	out 		= "../Tables/STable_to_4A_Allen2014_Differential.csv"
	df.to_csv(out,index=False)


	MDS2    	= dMD.mds_frame()
	MDS2.load_MD_score_file(map(add,["SRR1015583","SRR1015584"]), "DMSO")
	MDS2.load_MD_score_file(map(add,["SRR1015589", "SRR1015590"]), "TNFalpha")
	df      	= MDS2.differential_single("DMSO", "TNFalpha")
	out 		= "../Tables/STable_to_4B_Luo2014_Differential.csv"
	df.to_csv(out,index=False)

	
	MDS3    	= dMD.mds_frame()
	MDS3.load_MD_score_file(map(add,["SRR653421","SRR653422"]), "DMSO")
	MDS3.load_MD_score_file(map(add,["SRR653425","SRR653426" ]), "Estradiol")
	df      	= MDS3.differential_single("DMSO", "Estradiol")
	out 		= "../Tables/STable_to_4C_Hah2013_Differential.csv"
	df.to_csv(out,index=False)

	MDS4    	= dMD.mds_frame()
	MDS4.load_MD_score_file(map(add,["SRR1648886","SRR1648890", "SRR1648896", "SRR1648897", "SRR1648909"]), "DMSO")
	MDS4.load_MD_score_file(map(add,["SRR1648887","SRR1648891", "SRR1648898", "SRR1648899", "SRR1648910" ]), "DHT")
	df      	= MDS4.differential_single("DMSO", "DHT")
	out 		= "../Tables/STable_to_4D_Puc2015_Differential.csv"
	df.to_csv(out,index=False)
main()
