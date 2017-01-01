import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
'''
Required Files

1) ROC_TF_EMG.tsv
'''


def remove_duplicates(FP,TP):
	nFP, nTP 	= list(),list()
	i=0
	while i < len(FP):
		nFP.append(FP[i])
		nTP.append(TP[i])
		while i < len(FP) and TP[i]==nTP[-1]:
			i+=1

	return nFP, nTP

def main():
	FILE 	= "../files/ROC_TF_EMG.tsv"
	T 		= {}
	f 		= lambda x: float(x)
	H3_FP, H3_TP 		= list(),list()
	H3_DHS1, H3_DHS1 	= list(),list()

	with open(FILE) as FH:
		for line in FH:
			if len(line.strip("\n").split("\t"))==6:
				E, CT, TF, AUC, FP, TP 	= line.strip("\n").split("\t")
				if "OL" not in TF:
					if TF not in T:
						T[TF] 				= list()
					nFP, nTP 	= remove_duplicates(map(f, FP.split(",")) , map(f, TP.split(",")))
					T[TF].append((float(AUC), nFP, nTP ,CT  ))
	for TF in T:
		T[TF] 	= max(T[TF])
	#============================================================
	fontsize 	= 30
	BINS 		= 75
	vs 				=   sns.color_palette("Paired", n_colors=12, desat=.5)
	color 		= vs[9]
	sns.set_style("white")
	F 			= plt.figure(figsize=(16,8))
	ax 		= F.add_axes([0.1,0.1,0.4,0.8])
	ax2 		= F.add_axes([0.52,0.1,0.1,0.8])
	ax3 		= F.add_axes([0.61,0.1,0.01,0.8])

	sns.despine(left=False,bottom=False)
	sns.despine(left=True,bottom=False,ax=ax2)
	sns.despine(left=True,bottom=True,ax=ax3)
	ax2.set_yticks([])
	ax3.set_xticks([])
	ax2.set_ylabel("Area under the Curve\n(AUC)\n", fontsize=fontsize)
	auc_values 	= list()
	for tf in T:
		ax.plot(T[tf][1], T[tf][2], color=color,alpha=0.35,lw=2)
		auc_values.append(T[tf][0])
	
	counts,edges 	= np.histogram(auc_values, bins=BINS,range=(0,1))
	#labels 			= bin_labels(T, edges)
	ax3.yaxis.tick_right()		
	ax3.set_yticks((edges[1:] + edges[:-1])/2.0)
	#ax3.set_yticklabels(labels)
	ax3.set_yticklabels([])
	ax2.barh(edges[:-1], counts, edgecolor='white', color=color, align='center', height=(edges[-1]-edges[1])/BINS,lw=1.)
	ax2.set_ylim(0.0,1.0)

	ax2.set_xticks([0,max(counts)/2, max(counts)])
	ax2.set_xticklabels(["0",str(max(counts)/2 ), str(int(max(counts)))], fontsize=fontsize)
	ax2.set_xlabel("Frequency(N=" + str(len(T)) + ")", fontsize=fontsize)	
	ax.plot([0,1], [0,1], ls="--", lw=2.0, color=vs[3])
	ax.set_xticks([0,0.25,0.5,0.75, 1.0])
	ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=fontsize)
	
	ax.set_yticks([0,0.25,0.5,0.75, 1.0])
	ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=fontsize)
	ax.set_xlabel("False Positive Rate",fontsize=fontsize)
	ax.set_ylabel("True Positive Rate",fontsize=fontsize)
	ax.set_xlim(-0.05,1.0)
	plt.savefig('../svg_final/M2B_ROC.svg')
	plt.show()

main()