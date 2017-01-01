'''
	Requires:
		1) ../files/new_Allen2014_generated_DB_2_MM.txt !!!! This is an old version!! !!!!
	Outputs:
		1) 
'''
import numpy as np
def load_mds_db(FILE):
	X,FH,collect 	= list(), open(FILE, "r"),False
	for line in FH:
		if "Markov" in line:
			break
		elif "#Estimated" in line and "(~TSS)" in line:
			collect=True
		elif collect:
			X.append(str(len(X)+1)+","+ line.strip('\n').split("\t")[1] )
	return X
def main():
	FILE 	= "../files/new_Allen2014_generated_DB_2_MM.txt"
	TABLE = "../Tables/STable_to_3B_GC_bias.csv"
	X 		= load_mds_db(FILE)
	FHW 	= open(TABLE, "w") 	
	FHW.write("position_from_eRNA_origin,A,C,G,T\n")
	FHW.write("\n".join(X))

if __name__ == "__main__":
	main()
