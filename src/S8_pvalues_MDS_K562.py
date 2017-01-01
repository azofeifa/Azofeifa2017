'''
	Requires:
		1) ../Tables/STable_to_3BC_pvalues_MDS_K562.csv (ST_M3BC_pvalues_K562.py)
	Outputs:
		svg figure showing pvalue stats of this table
'''
import pandas as pd,matplotlib.pyplot as plt
import seaborn as sns


def main():
	df 	= pd.read_csv("../Tables/STable_to_3BC_pvalues_MDS_K562.csv")
	sns.set(font_scale=2)
	sns.set_style("white")
	F 			= plt.figure(tight_layout=True,figsize=(15,6))
	ax1,ax2 	= F.add_subplot(1,2,1),F.add_subplot(1,2,2)
	B 			= 50
	ax1.set_title("Stationary Background Model")
	ax1.hist(df["Stationary p-value"],bins=B,edgecolor="white")
	ax2.hist(df["Non Stationary p-value"],bins=B,edgecolor="white")
	ax2.set_title("Non Stationary Background Model")
	ax2.set_ylim(ax1.get_ylim())
	for ax in (ax1,ax2):
		ax.set_xlabel('p-value (one-tailed)' )
		ax.set_ylabel('frequency' )
		sns.despine(ax=ax)
		ax.yaxis.grid()
	plt.savefig("../svg_final/S1_pvalue_histograms_MDS_K562")
	plt.show()
if __name__ == "__main__":
	main()