import pandas as pd,matplotlib.pyplot as plt, seaborn as sns
import numpy as np
'''
	Requires:
		1) ../Tables/STable_to_1D_ChIP_BIC_displacements.csv
	Outputs:
		1) figure showing stats of this table
'''
def draw(df, F_OUT,REP=True,THRESHOLD=100):
	sns.set_style("white")
	F 			= plt.figure(facecolor="white",figsize=(17,6),tight_layout=True)
	ax 		= plt.gca()
	names 	= list(df.name)
	N 			= float(len(df))
	dBIC 		= zip(df.dBIC, range(len(df.index)))
	disps 	= zip(df.pk_displacement,range(len(df.index)))
	if REP:
		names 	= list(set(df.name))
		dBIC 		= [ (np.mean(df[df.name==n].dBIC),i)  for i,n in enumerate(names) ]
		disps 	= [ (np.mean(df[df.name==n].pk_displacement),i)  for i,n in enumerate(names) ]
	dBIC.sort(reverse=True)
	percent 		= float(len([j for j,(x,i) in enumerate(dBIC) if x > THRESHOLD])) / len(dBIC)
	ax.bar( [j for j,(x,i) in enumerate(dBIC) if x > THRESHOLD], [x for x,i in dBIC if x > THRESHOLD] ,
		width=0.5, color="red",label="Considered Bimodial Binding\n " + str(percent*100)[:5] + "%" )
	ax.bar( [j for j,(x,i) in enumerate(dBIC) if x <= THRESHOLD], [x for x,i in dBIC if x <= THRESHOLD] ,width=0.5)
	ax.legend(loc="best",fontsize=20)
	ax.set_ylabel(r"$\Delta BIC$",fontsize=25)
	if not REP:
		ax.set_xlabel("ChIP Regulator (TF-binding)\nIncluding Replicates",fontsize=15)
	else:
		ax.set_xlabel("ChIP Regulator (TF-binding)",fontsize=15)
	sns.despine(ax=ax)
	plt.savefig(F_OUT)
	plt.show()


def main():
	FILE 	= "../Tables/STable_to_1D_ChIP_BIC_displacements.csv"
	F_OUT = "../svg_final/S1_BIC_DISP.svg"
	df 	= pd.read_csv(FILE)

	draw(df,F_OUT)
if __name__ == "__main__":
	main()
