'''
	Required: 
		1) ../Tables/STable_to_1D_ChIP_BIC_displacements.csv 	(ST_M1D_BIC_ChIP_binding.py)
		2) ../Tables/STable_to_1C_eRNA_TF_mark_ref.csv 			(ST_M1C_eRNA_TF_mark_Ref_overlaps.py)
	Outputs:
		1) Figure that shows Presence/Abscence matrix
			and TF peak Displacement from Bidir epicenter
'''
import matplotlib.pyplot as plt
import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.cm as cm
import scipy.cluster.hierarchy as sch
from sklearn.cluster import KMeans
from scipy.spatial.distance import pdist,squareform
import pandas as pd
import time
import seaborn as sns
def Kmeans(D2,K=10):
	print "fitting",
	clf 		= KMeans(init='k-means++', n_clusters=K)
	clf.fit(D2)
	print "done"

	labels 	= clf.labels_
	U 			= {}
	T 			={}
	for i,l in enumerate(labels):
		if l not in T:
			T[l] 	= list()
		T[l].append(i)
	S 		= [  i for l in T for i in T[l]]
	return S



def main():
	sns.set_style("white")
	F 		= plt.figure(figsize=(15,8), facecolor="white")
	ax1 	= F.add_axes([0.1,0.1,0.69,0.8])
	ax2 	= F.add_axes([0.83,0.1,0.1,0.8])
	df 	= pd.read_csv("../Tables/STable_to_1C_eRNA_TF_mark_ref.csv")
	df2 	= pd.read_csv("../Tables/STable_to_1D_ChIP_BIC_displacements.csv")
	df 	= df[(df.TSS == 0) & (df.chrom != "chrY") ]
	nA 	= list()
	for c in df.columns[22:]:
		nDF 	= np.array(df2[df2.name==c.split(".")[0]].as_matrix()[:,4:],dtype=float)
		nA.append(np.mean(nDF,axis=0))
	A 		= np.matrix(df.as_matrix()[:,22:],dtype=float).T
	idx 	= Kmeans(A, K=20)
	nA 	= np.array([nA[i]/np.max(nA[i]) for i in idx])

	idx2 	= Kmeans(A.T, K=20)
	A 		= A[idx,:]
	A 		= A[:,idx2]

	ax1.imshow(A,vmin=0,vmax=1,aspect="auto", cmap=cm.Blues, interpolation="nearest")
	ax2.imshow(nA,vmin=0,vmax=1,aspect="auto", cmap=cm.GnBu, interpolation="nearest")

	sns.despine(ax=ax1,left=True,bottom=True)
	sns.despine(ax=ax2,left=True,bottom=True)
	for ax in (ax1,ax2):
		ax.set_xticks([])
		ax.set_yticks([])
	ax1.set_ylabel("Regulator (ENCODE)",fontsize=25)
	ax1.set_xlabel("eRNA (K562)",fontsize=25)
	ax2.set_xlabel("Displacement of Regulator\nfrom eRNA origin",fontsize=25)

	plt.savefig("../svg_final/M1C_Combinatorics_TF_K562.svg")

	plt.show()
main()
