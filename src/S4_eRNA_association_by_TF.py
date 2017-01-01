'''
	Requires:
		1) ../Tables/STable_to_2A_eRNA_association_by_TF.csv (ST_M2AB_eRNA_association_by_TF.py)
	Outputs:
		A bar chart of this table
'''
import matplotlib.pyplot as plt, seaborn as sns
import numpy as np,pandas as pd


def main(REP=True):
	TABLE 	= "../Tables/STable_to_2A_eRNA_association_by_TF.csv"
	df 		= pd.read_csv(TABLE)
	if REP:
		names = list(set(df.name))
		y 		= [np.mean(df[df.name.isin([n]) ].percent_eRNA) for n in names]
		labels=names
	else:
		y 			= [min(x,1) for x in list(df.percent_eRNA)]
		labels 	= list(df.name)
	y 				= [(y*100,i)for i,y in enumerate(y)]
	y.sort(reverse=True)
	idx 			= [j for i,j in y]
	labels 		= [labels[i] for i in idx]
	y 				= [i for i,j in y]
	x 			= range(0,len(y))
	sns.set_style("white")
	F 			= plt.figure(figsize=(15,7),tight_layout=True)
	ax 		= plt.gca()
	ax.bar(x,y,edgecolor="white",width=0.5)
	sns.despine(ax=ax)
	ax.set_ylabel("Percent of Binding Events associated with an eRNA",fontsize=15)
	ax.set_xlabel("Regulator (ChIP-seq, ENCODE)",fontsize=15)
	ax.set_xlim(0,130)
	plt.savefig("../svg_final/S1_PercentChIPWitheRNA.svg")
	plt.show()

if __name__ == "__main__":
	main()