'''
	Requires:
		1) ../Tables/STable_to_3E_MDS_N_across_all.csv (ST_M3E_MD_scores_across_all)
	Outputs:
		Tornado plot (histogram from MD score across 499)
		And clustergram
'''
import pandas as pd,numpy as np
import matplotlib.cm as cm
import scipy.cluster.hierarchy as sch

import seaborn as sns,matplotlib.pyplot as plt
def main():
	df 		= pd.read_csv("../Tables/STable_to_3E_MDS_N_across_all.csv")
	sns.set_style("white")
	F 				= plt.figure(figsize=(10,10))
	ax1,ax2 		= F.add_axes([0.1,0.25,0.3,0.6]),F.add_axes([0.6,0.25,0.3,0.6])
	ax1.set_ylabel("Motif Model(n=641)",fontsize=20)
	ax2.set_xlabel("Nascent Transcript Dataset (n=494)",fontsize=20)
	axd1 			= F.add_axes([0.54,0.25,0.05,0.6])

	X 			= df.as_matrix()
	X 			= np.array(X[:,1:494],dtype=float)
	center 	= np.median(X,axis=1)
	scale 	= np.std(X,axis=1)
	X 			= np.array([(X[i,:] - center[i])/scale[i] for i in range(X.shape[0])])

	Y1 		= sch.linkage(X, method='ward')
	with plt.rc_context({'lines.linewidth': 0.5}):
		Z1 		= sch.dendrogram(Y1,  orientation="right",
											color_threshold=0.0)
	idx1 		= Z1["leaves"]

	Y2 		= sch.linkage(X.T, method='ward')
	axd2 		= F.add_axes([0.6,0.86,0.3,0.05])
	with plt.rc_context({'lines.linewidth': 0.5}):
		Z2 		= sch.dendrogram(Y2,  orientation="top",
											color_threshold=0.0)
	idx2 		= Z2["leaves"]

	X 			= X[idx1,:]
	X 			= X[:,idx2]

	B 			= 50
	Y 			= np.array([np.array(np.histogram(map(float, X[i,:]),bins=B)[0],dtype=float) for i in range(X.shape[0])])
	Y 			= [y/float(max(y) )for y in Y]
	cb1 		= ax1.imshow(Y,aspect="auto",interpolation="nearest",cmap=cm.GnBu)
	cb2 		= ax2.imshow(X,aspect="auto",interpolation="nearest",
						cmap=cm.seismic,vmin=-3.5,vmax=3.5)
	axcb1,axcb2 = F.add_axes([0.1,0.15,0.3,0.015]),F.add_axes([0.6,0.15,0.3,0.015])
	cb1 			= plt.colorbar(cb1,orientation="horizontal",cax=axcb1)
	cb1.set_label('Frequency across 494 Datasets\n(Normalized to Max)',fontsize=15)
	ax1.set_xlabel("MD Score\n(centered and scaled)")
	cb2 			= 	plt.colorbar(cb2,orientation="horizontal",cax=axcb2)
	cb2.set_label('MD Score\n(Centered by Mean\nScaled by Standard Deviation)',fontsize=15)
	cb1.outline.set_visible(False);cb2.outline.set_visible(False)
	for i,ax in enumerate((ax1,ax2,axd1,axd2)):
		ax.set_yticks([])
		if i > 0:
			ax.set_xticks([])
		else:
			vals 	= [-3.2,-2.4,-1.6,-0.8,0,0.8,1.6,2.4, 3.2]
			ax.set_xticks(np.linspace(0,B,len(vals) ))
			ax.set_xticklabels(vals)
		sns.despine(ax=ax,left=True,bottom=True)
	plt.savefig("../svg_final/S1_MDS_across_clustergram.svg")
	plt.show()


if __name__ =="__main__":
	main()

