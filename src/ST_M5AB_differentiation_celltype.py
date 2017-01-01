'''
	Requires:
		1) /Users/joazofeifa/Lab/new_motif_distances/mouse/150/
		2) motif_displacement_analysis package (MDAP)
	Outputs:
		1) STable_to_4F_KLA_Treatment.csv
		2) STable_to_4E_Flavopiridol_Treatment.csv
'''


import sys
sys.path.append("/Users/joazofeifa/Lab/motif_displacement_analysis/MDAP/")
import differential_MD as dMD

def main():
	AP    = "_MDS.tsv"
	MD_DIR="/Users/joazofeifa/Lab/new_motif_distances/human/300/"
	f     = lambda x: MD_DIR+x+AP 


	MDS2    = dMD.mds_frame()

	MDS2.load_MD_score_file( [ f("SRR1552482"),f("SRR1552483"), f("SRR1552485") ] , "GM12878"  )
	MDS2.load_MD_score_file( [ f("SRR1552480"),f("SRR1552481"),f("SRR1554311"),f("SRR1554312")    ] , "K562"  )

	df      	= MDS2.differential_single("GM12878", "K562")
	out 		= "../Tables/STable_to_5A_GM12878_vs_K562.csv"
	df.to_csv(out,index=False)




	MDS    = dMD.mds_frame()
	MDS.load_MD_score_file( [ f("SRR1145801")  ] , "hESC"  )
	MDS.load_MD_score_file( [ f("SRR1145808")   ] , "endoderm"  )
	MDS.load_MD_score_file( [ f("SRR1145815")  ] , "primitive gut tube" )
	MDS.load_MD_score_file( [ f("SRR1145822") ] , "posterior foregut"  )
	MDS.load_MD_score_file( [ f("SRR1145829") ] , "pancreatic endoderm" )

	
	df 		= MDS.differential_multiple([	"hESC", "endoderm", "primitive gut tube", 
														"posterior foregut", "pancreatic endoderm"])
	out 		= "../Tables/STable_to_5B_hESC_to_ENDODERM.csv"
	df.to_csv(out,index=False)
if __name__ == "__main__":
	main()

