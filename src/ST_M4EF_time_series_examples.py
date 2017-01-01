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
	AP    	= "_MDS.tsv"
	MD_DIR 	= "/Users/joazofeifa/Lab/new_motif_distances/mouse/150/"
	f     	= lambda x: MD_DIR+x+AP 


	DMSO 	= [["SRR930649", "SRR930650", "SRR930651", "SRR930652" ],
				["SRR930653", "SRR930654"],
				["SRR930655", "SRR930656"],
				["SRR930657", "SRR930658"]]
	TREAT = [["SRR930659", "SRR930660", "SRR930661"],
				["SRR930663", "SRR930664"],
				["SRR930665", "SRR930666"],
				["SRR930667", "SRR930668" ]]
	IDS 	= "1hr", "6hr", "12hr", "24hr"


	MDS    = dMD.mds_frame()
	for i,SRRS in enumerate(DMSO):
		MDS.load_MD_score_file(map(f, SRRS), IDS[i] + "(DMSO)", DISPS=False  )

	for i,SRRS in enumerate(TREAT):
		MDS.load_MD_score_file(map(f, SRRS), IDS[i] + "(TREAT)", DISPS=False  )
	
	df 		= MDS.differential_multiple([ID + "(TREAT)" for ID in IDS],
	                        base=[ID + "(DMSO)" for ID in IDS])
	out 		= "../Tables/STable_to_4F_KLA_Treatment.csv"
	df.to_csv(out,index=False)


	MDS2    = dMD.mds_frame()

	MDS2.load_MD_score_file( [f("SRR" + str(i) ) for i in range(935093, 935097)] , "none" )

	MDS2.load_MD_score_file( [f("SRR" + str(i) )for i in range(935097, 935101)] , "2(min)" )
	MDS2.load_MD_score_file( [f("SRR" + str(i) )for i in range(935101, 935105)] , "5(min)" )
	MDS2.load_MD_score_file( [f("SRR" + str(i) )for i in range(935105, 935109)] , "12.5(min)" )
	MDS2.load_MD_score_file( [f("SRR" + str(i) )for i in range(935109, 935113)] , "25(min)" )
	MDS2.load_MD_score_file( [f("SRR" + str(i) )for i in range(935113, 935117)] , "50(min)" )

	df 		= MDS2.differential_multiple(["none","2(min)", "5(min)", "12.5(min)", "25(min)", "50(min)" ])
	out 		= "../Tables/STable_to_4E_Flavopiridol_Treatment.csv"
	df.to_csv(out,index=False)
main()

