'''
	Requires:
		../Tables/STable_to_3_mds_over_ChIP.csv (ST_M3_mds_over_ChIP.py)
	Outputs:
		A figure that displays the motif displacement distribution over eRNAs
		associated with a TF binding site and eRNAs lacking TF binding site;
		boxplot and scatter of the MD scores
'''
import pandas as pd
import seaborn as sns, matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
def compute_mds(d, percent=0.1):
	n 	= len(d)
	c 	= n/2
	w 	= int(n*(percent/2))
	N 	= float(sum(d))
	S 	= sum(d[c-w:c+w])
	return S/N

def main():
	table 	= "../Tables/STable_to_3_mds_over_ChIP.csv"
	df 		= pd.read_csv(table)
	non,over = df.lack_hist,df.over_hist
	over,non = [np.array(map(float, o.split("|"))) for o in over],[np.array(map(float, o.split("|"))) for o in non]
	mdsover 	= np.array(map(compute_mds, over))
	mdsnon 	= np.array(map(compute_mds, non))
	over,non = [o/max(max(o),max(non[i])) for i,o in enumerate(over) ],[o/max(max(o),max(over[i])) for i,o in enumerate(non) ]
	sns.set_style("white")
	
	F 			= plt.figure(figsize=(10,10))
	ax1,ax2 	= F.add_axes([0.1,0.1,0.14,0.8]),F.add_axes([0.25,0.1,0.14,0.8])
	axcbar 	= F.add_axes([0.4,0.1,0.015, 0.8])
	axbox 	= F.add_axes([0.55,0.55,0.39, 0.39])
	axscatter= F.add_axes([0.55,0.1,0.39, 0.39])
	'''
		boxplots
	'''
	df_mds 		= pd.DataFrame(columns=("over", "non"))
	for i in range(len(mdsover)):
		df_mds.loc[i] 	= mdsover[i],mdsnon[i]
	sns.boxplot(data=df_mds,ax=axbox,width=0.2)
	axbox.yaxis.grid()
	axbox.set_xticklabels(["TF-bound eRNAs\n(ChIP-peak)","TF-unbound eRNAs\n(no ChIP-peak)"],fontsize=15)
	sns.despine(ax=axbox)
	axbox.set_ylabel("Motif Displacement Score (MDS)",fontsize=15)
	
	'''
		scatter of MDS (over) vs MDS (non)
	'''
	axscatter.scatter(mdsnon, mdsover)
	axscatter.plot([0,1],[0,1],lw=2,color="g", ls="--",label="y=x")
	axscatter.set_xlim(0,1)
	axscatter.set_ylim(0,1)
	axscatter.set_xlabel('MDS\n(TF-unbound eRNAs)',fontsize=15)
	axscatter.set_ylabel('MDS\n(TF-bound eRNAs)',fontsize=15)
	sns.despine(ax=axscatter)


	'''
		heatmaps of motif displacement distributions
	'''

	im 		= ax1.imshow(over,cmap=cm.GnBu,interpolation="nearest", aspect="auto")
	im 		= ax2.imshow(non,cmap=cm.GnBu,interpolation="nearest", aspect="auto")
	cbar 		= plt.colorbar(im, cax=axcbar)
	ax1.set_title("TF-bound eRNAs\n(ChIP-peak)")
	ax2.set_title("TF-unbound eRNAs\n(no ChIP-peak)")
	ax1.set_ylabel('TF (ChIP-seq, ENCODE) + matched motif Model',fontsize=15)
	for ax in (ax1, ax2):
		sns.despine(ax=ax, left=True,bottom=True)
		ax.set_yticks([])
		ax.set_xticks([])
		ax.set_xlabel("eRNA origin\n(3KB Window)")
	cbar.set_ticks([])
	cbar.set_label("Frequency of Motif near eRNA origins",fontsize=15)
	cbar.outline.set_visible(False)
	plt.savefig("../svg_final/S6_md_over_ChIP.svg")
	plt.show()

if __name__=="__main__":
	main()