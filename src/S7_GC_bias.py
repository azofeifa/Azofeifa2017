'''
	Requires:
		1) ../Tables/STable_to_3B_GC_bias.csv (ST_M3B_GC_bias.py)
	Outputs:
		1) Figure showing this table; nucleotide bias across eRNAs
'''
import pandas as pd,numpy as np
import seaborn as sns,matplotlib.pyplot as plt
def main():
	TABLE 	= "../Tables/STable_to_3B_GC_bias.csv"
	df 		= pd.read_csv(TABLE).as_matrix()[:,1:]
	labels 	= "A", "C", "G","T"
	colors 	= "green", "blue", "purple", "red" 
	sns.set(font_scale=2)
	sns.set_style("white")
	K 			= 5
	positions= np.arange(0,df.shape[0],K)
	F 			= plt.figure(tight_layout=True)
	ax 		= plt.gca()
	bg 		= np.mean(df[:100,:],axis=0)

	for i in range(4):
		ax.plot(positions,[np.mean(df[j:j+K,i])-bg[i] for j in positions], 
					label=labels[i],color=colors[i])
	ax.set_xticks([])
	ax.set_xlabel("eRNA origin\n(3KB Window)" )
	ax.set_ylabel("Difference from\nBackground Expectation" )
	sns.despine(ax=ax)
	ax.legend(loc="best",fontsize=15)
	plt.savefig("../svg_final/S1_GC_bias.svg")
	plt.show()
if __name__ == "__main__":
	main()
