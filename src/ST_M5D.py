import pandas as pd,numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss
import math,time 
def compute(df):
	'''
		change to either
		1) tissue
		2) general_celltype
		3) specific_celltype
	'''
	ct_type 	= "tissue"
	t_code  	= "treatment_code_2"

	tfs 		= dict([(x,sum(getattr(df, x))) for x in list(df.columns) if "HO"==x.split("_")[0]])

	'''
		0) just general counts
	'''

	cts 		= dict([(ct,float(len(df[getattr(df, ct_type)==ct ]))) for ct in list(set(getattr(df,ct_type)))]) #get unique set
	exps 		= dict([(e,float(len(df[getattr(df, t_code)==e ]))) for e in list(set(getattr(df,t_code)))]) #get unique set
	print ".",
	'''
		1) probability of a cell type give an experiment; precompute this
	'''
	
	ct_e 		= dict([(ct,dict([(e,len(df[ (getattr(df, ct_type)==ct) & (getattr(df, t_code)==e)  ])/exps[e]  ) for e in exps]))  for ct in cts])

	print ".",

	'''
		2) probability of the experiment given the TF; precompute this
	'''

	e_tf 		= dict([(e,dict([(tf, sum(getattr(df[getattr(df, t_code)==e], tf))/tfs[tf]  ) for tf in tfs]) )for e in exps])
	print ".",


	m,n 		= len(cts),len(tfs)
	A 			= np.zeros((m,n)) #going to test the association between motif_i and ct_j

	'''
		according to the model, we need the probability of experiment given a TF
		and the probility of that experiment given the cell type;and then we are going to sum
	'''
	PARAMS 			= list()
	for tf in tfs:
		for ct in cts:
			p 			= 0.0
			for e in exps:
				'''
					1) probability of a cell type give an experiment
				'''
				f1 	= ct_e[ct][e]
				'''
					2) probability of the experiment given the TF
				'''
				f2 	= e_tf[e][tf]

				p+=f1*f2


			n 		= tfs[tf]
			EX 	= p*n
			B 		= ss.norm(EX, EX*(1.0-p) )
			obs  	= sum(getattr( df[getattr(df, ct_type)==ct],tf))
			pv 	= B.cdf(obs)
			ENR 	= obs/float(cts[ct])
			if math.isnan(pv):
				pv = 0.5
			PARAMS.append((tf, ct, p,n,obs,pv,ENR))

	nDF 	= pd.DataFrame(PARAMS,columns=("motif", "ct" , "bin_p","bin_n", "obs", "pvalue", "ENR"))
	nDF.to_csv("../Tables/STable_to_5D_cell_tf_associations.csv",index=False)








def main():
	'''
		READ: This file contains for each SRR experiment and 
		for each TF motif model, the number of times that 
		motfs MD score was significantly increased in favor that SRR
		experiment (this number maxes out at the total number of 
		possible comparisons,i.e. number of SRRs)
	
		Currently this file comes from M_Figure_5C.make_count_table() but should change
	'''
	FILE 			= "/Users/joazofeifa/Lab/Article_drafts/EMG_paper/Tables/STable_to_5_SRR_tf_counts.csv"

	counts 		= pd.read_csv(FILE)
	counts 		= counts.set_index(counts.SRAnumber)
	df 			= pd.read_csv("../files/conditions_table_joey.csv")
	df 			= df.set_index(df.SRAnumber)
	df 			= df[df.organism=="human"]

	df 			= df.merge(counts)
	df 			= df.set_index(df.SRAnumber)
	compute(df)


main()

