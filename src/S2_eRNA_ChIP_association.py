'''
	Requires:
		1) ../Tables/STable_to_1C_eRNA_TF_mark_ref.csv (ST_M1B_eRNA_mark_overlap.py)
	Outputs:
		2) S1_eRNA_mark_association.svg
			bar graph showing percent of eRNAs associated with a mark
'''

import pandas as pd
import seaborn as sns,numpy as np
import matplotlib.pyplot as plt
def draw(df,st=4,sp=22,
				F_OUT="../svg_final/S1_eRNA_mark_association.svg",
				fontsize=15,xlabel="",REP=True):

	df 		= df[df["TSS"]==0] #filter out promoter association
	sns.set_style("white")
	F 			= plt.figure(facecolor="white",figsize=(15,6),tight_layout=True)
	ax 		= plt.gca()
	N 			= float(len(df))
	marks  	= df.columns[st:sp]
	vals 		= [(100*sum(getattr(df, col))/N,col)  for col in marks]
	if REP:
		marks2= set([m.split(".")[0] for m in marks])
		vals 	= [(np.max([ vals[i][0] for i,col in enumerate(marks) if col.split(".")[0]==m]),m) for m in marks2 ]
	vals.sort(reverse=True)
	marks 	= [y for x,y in vals]
	ax.bar( range(len(marks)), [x for x,y in vals] ,width=0.5,  )
	if fontsize > 10:
		ax.set_xticks(np.arange(0.25, len(marks) ,1))
		ax.set_xticklabels(marks,rotation=90,fontsize=fontsize)
	ax.set_ylabel("Percent of eRNAs associated with Mark",fontsize=15)
	ax.set_xlim(-0.5,len(marks))
	ax.set_xlabel(xlabel,fontsize=20)
	sns.despine(ax=ax)
	plt.savefig(F_OUT)
	plt.show()



def main():
	df 	= pd.read_csv("../Tables/STable_to_1C_eRNA_TF_mark_ref.csv")
	
	draw(df,st=4,sp=25, 
		F_OUT="../svg_final/S1_eRNA_mark_association.svg",
		xlabel="Chromatin Mark" )
	draw(df,st=22,sp=len(df.columns), 
		F_OUT="../svg_final/S1_eRNA_TF_association.svg" ,fontsize=3,
		xlabel="ChIP Regulator (TF-binding)" )



if __name__ == "__main__":
	main()