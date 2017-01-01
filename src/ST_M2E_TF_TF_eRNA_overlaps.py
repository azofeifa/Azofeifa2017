import S_Table_1_eRNA_TF_mark_Ref_overlaps as STP
import os,numpy as np
'''
	This requires:
		1) ../files/ChIP_files/K562/ 
		2) ../files/marks/K562/
		3) ../files/K562_EMG.bed
		4) ../files/TSS.bed
	Outputs:
		1) ../Tables/STable_to_2C_TF_mark_overlaps_eRNA_counts.csv
	
		for each TF the counts association with a mark
		conditioned on eRNA presence
	Note:
		!!This will take a long time!!
'''

def compute(DIR,MARKS,IDS,OUT):
	FILES 	= [f for f in os.listdir(DIR)]
	META 		= STP.load_meta(DIR+"metadata.tsv")
	FILES 	= [f for f in FILES if f.split(".")[0] in META and "eGFP" not in META[f.split(".")[0]]]
	FHW 		= open(OUT, "w")		
	f1,f2 	= lambda x: x+"_eRNA",lambda x: x+"_no_eRNA"	
	FHW.write("name,total(non)," + ",".join(map(f1,IDS[2:])) + "," + ",".join(map(f2,IDS[2:])) + "\n")
	for f in FILES:
		G 		= STP.load_bed_file(DIR+f, AA=len(MARKS))
		for i,m in enumerate(MARKS):
			G 	= STP.overlap(G,m,i)
		A 		= np.array([x[2] for chrom in G for x in G[chrom] ])
		A 		= A[A[:,0]==0,:] #filter out the promoters !!! 
		line 	= str(META[f.split('.')[0]]) + "," + str(A.shape[0]) + ","
		row1 	= np.sum(A[A[:,1]==1,:],axis=0) #condition on eRNA presence
		row2 	= np.sum(A[A[:,1]==0,:],axis=0) #the eRNA isn't there!
		line	+= ",".join(map(str, row1[2:])) + ","
		line	+= ",".join(map(str, row2[2:])) + "\n"
		FHW.write(line)
		FHW.flush()

def main():
	TF_ChIP_DIR 		= "../files/ChIP_files/K562/" #directory of bed files
	BIDIR 				= "../files/K562_EMG.bed"
	Ref 					= "../files/TSS.bed"


	OUT 					= "../Tables/STable_to_2E_TF_TF_overlaps_eRNA_counts.csv"
	
	R 						= STP.load_bed_file(Ref)
	B 						= STP.load_bed_file(BIDIR)
	MARKS, IDS 			= STP.get_bed_files(TF_ChIP_DIR, OUT=None)
	IDS 					= ["eRNA", "Ref"] + IDS

	compute(TF_ChIP_DIR, MARKS, IDS,OUT)
if __name__=="__main__":
	main()
