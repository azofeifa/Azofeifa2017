import M_Figure_5B_triangles as mf5b
import pandas as pd
import matplotlib.pyplot as plt,numpy as np
from matplotlib import cm as cm
import scipy.cluster.hierarchy as sch
import matplotlib as mpl
from sklearn.cluster import KMeans

class color:
	def __init__(self,D,cmap=cm.Set1):
		self.data 	= D
		self.K 		= dict([(d,i) for i,d in enumerate(set(self.data))])
		self.norm 	= mpl.colors.Normalize(vmin=0, vmax=len(self.K))
		self.cmap 	= cmap
		self.m 		= mpl.cm.ScalarMappable(norm=self.norm, cmap=self.cmap)
	def get_color(self, k):
		return self.m.to_rgba(self.K[k])
def Kmeans(D2,K=10):
	clf 		= KMeans(init='k-means++', n_clusters=10, n_init=10)
	clf.fit(D2)

	labels 	= clf.labels_
	U 			= {}
	T 			={}
	for i,l in enumerate(labels):
		if l not in U:
			U[l] 	= 0
			T[l] 	= list()
		U[l]+=float(np.max(D2[i,:]))
		T[l].append(i)
	S 		= [ ( U[l]/float(len(T[l])) , l ) for l in U]
	S.sort(reverse=True)
	U 		= list()
	for x,l in S:
		U+=T[l]
	return U



def make_count_table(OUT,df):
	'''
		change to either
		1) tissue
		2) general_celltype
		3) specific_celltype
	'''
	ct_type 	= "tissue"


	FH 			= open(OUT, "r")
	lines 		= FH.readlines()
	SRRS 			= lines[0].strip("\n").split(',')
	MOTIFS 		= dict([(str(i),m) for i,m in enumerate(lines[1].strip("\n").split(','))])
	cts 			= set(getattr(df, ct_type))
	exps 			= set(getattr(df, "treatment_code"))
	counts_ct 	= dict([ (m, dict([ (ct, 0.0) for ct in cts ])) for m in MOTIFS.values()  ])
	counts_exp 	= dict([ (m, dict([ (exp, 0.0) for exp in exps ])) for m in MOTIFS.values()  ])
	counts_SRR 	= dict([ (m, dict([ (SRR, 0.0) for SRR in SRRS ])) for m in MOTIFS.values()  ])
	for i,l in enumerate(lines[2:]):
		SRRi 			= SRRS[i]
		cti 			= getattr(df.loc[SRRi], ct_type)
		expi 			= getattr(df.loc[SRRi], "treatment_code")
		for j,motifs in enumerate(l.strip("\n").split("\t")):
			'''
				down means differentiall active in favor SRRi
				up means differentiall active in favor SRRj
			'''
			downs,ups= motifs.split("|")
			for u in downs.split(","):
				if u:	
					counts_ct[MOTIFS[u]][cti]+=1
					counts_exp[MOTIFS[u]][expi]+=1
					counts_SRR[MOTIFS[u]][SRRi]+=1
	df1 			= pd.DataFrame(counts_ct)
	df1.index.name= 	"ct"
	df1.to_csv("../Tables/STable_to_5_cell_tf_counts.csv")

	df1 			= pd.DataFrame(counts_exp)
	df1.index.name= 	"exp"
	df1.to_csv("../Tables/STable_to_5_exp_tf_counts.csv")

	df1 			= pd.DataFrame(counts_SRR)
	df1.index.name= 	"SRAnumber"
	df1.to_csv("../Tables/STable_to_5_SRR_tf_counts.csv")
	return True
def visualize_raw_table(df):
	ct 		= "tissue"
	T 			= pd.read_csv("../Tables/STable_to_5_SRR_tf_counts.csv")
	A 			= np.array(T.as_matrix()[:,1:],dtype=float)
	cts 		= [getattr(df.loc[s],ct) for s in T.SRAnumber]
	C 			= color(cts)	

	F 			= plt.figure(figsize=(15,7), facecolor="white")
	ax1 		= F.add_axes([0.18,0.1,0.015,0.8])


	ax2 		= F.add_axes([0.2,0.1,0.4,0.8])
	Y1 		= sch.linkage(A, method='ward')
	Z1 		= sch.dendrogram(Y1,  orientation="right",
										color_threshold=0.0,no_plot=True)
	idx1 		= Z1["leaves"]
	Y2 		= sch.linkage(A.T, method='ward')
	Z2 		= sch.dendrogram(Y2,  orientation="right",
										color_threshold=0.0,no_plot=True)
	idx2 		= Z2["leaves"]
	idx2 		= Kmeans(A.T,K=80)
	cts 		= [(ct,i)for i,ct in enumerate(cts)]
	cts.sort()
	idx1 		= [y for x,y in cts]
	cts 		= [x for x,y in cts]
	A 			= A[idx1,:]
	A 			= A[:,idx2]

	for j,i in enumerate(idx1):
		ax1.barh([A.shape[0]-j],[1],color=C.get_color(cts[j]),edgecolor=C.get_color(cts[j]))
	ax1.set_xlim(0,1)
	ax1.set_ylim(0,A.shape[0])

	ax2.imshow(A,vmin=A.min(),vmax=200, 
						cmap=cm.Blues, aspect="auto", interpolation="nearest")
	mf5b.despine(ax2,left=True,bottom=True)
	mf5b.despine(ax1,left=True,bottom=True)
	ax1.set_xticks([])
	ax2.set_yticks([])
	plt.show()






def main():
	regenerate 	= False #generate pairwise distance matrix; long time
	count_table = False #generate pairwise distance matrix; long time
	OUT 			= "../Tables/STable_to_5_pairwise_differences_all.csv"
	df 			= pd.read_csv("../files/conditions_table_joey.csv")
	df 			= df.set_index(df.SRAnumber)
	df 			= df[df.organism=="human"]
	if regenerate:
		MDS_DIR 							= "/Users/joazofeifa/Lab/new_motif_distances/human/300/"
		MDS,CTS, SRRS, MOTIFS,DF 	= mf5b.load_MDS_files(MDS_DIR,df)
		mf5b.make_distance_matrix(MDS, CTS, SRRS, MOTIFS, df,OUT)
	if count_table:
		make_count_table(OUT,df)
	visualize_raw_table(df)



main()