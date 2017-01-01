import pandas as pd,os,numpy as np,math
from scipy.special import erf
import matplotlib.pyplot as plt,matplotlib as mpl
import time
def despine(ax,left=False, bottom=False):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	if left:
		ax.spines['bottom'].set_visible(False)
	if bottom:
		ax.spines['left'].set_visible(False)
def compute_pvalue(P1,P2,N1,N2, mean=0,threshold=pow(10,-7), diff_threshold=0.1):
	if not N1 or not N2:
		return None
	p       = (P1*N1 + P2*N2)  / (N1+N2)
	SE      = math.sqrt(p*(1-p)*( (1.0/N1) + (1.0/N2)  ))
	if not SE :
		return None
	top 		= (np.mean(P2)-np.mean(P1)-mean)
	Z       	= top / SE
	pval    	= 0.5*(1.0 + erf(Z / math.sqrt(2.0) ))
	sign 		= None
	if pval < threshold and abs(top) > diff_threshold:
		sign 	= 0
	elif pval > 1-threshold and abs(top) > diff_threshold:
		sign 	= 1
	return sign
def load_MDS_files(DIR,df, TYPE=0, MOTIFS=None):
	G 	=   {}
	for i,f in enumerate(os.listdir(DIR)):
		if f.split("_MDS")[0] in df.index:
			M 	= {}
			N 	= {}
			for line in open(DIR+f):
				if "Binned" in line:
					break
				elif line[0]!="#":
					line_array 	= line.split("\t")
					m,n,x 			= line_array[0],line_array[1],line_array[2]
					M[m] 				= map(float, x.split(","))[TYPE]
					N[m] 				= map(float, n.split(","))[TYPE]
			G[f.split("_MDS")[0]]=M,N
		else:
			print f.split("_MDS")[0]
	SRRS 	= G.keys()
	if MOTIFS is None:
		MOTIFS= G[SRRS[0]][0].keys()

	MDS 		= np.zeros((len(SRRS), len(MOTIFS)))
	CTS 		= np.zeros((len(SRRS), len(MOTIFS)))
	for i in range(MDS.shape[0]):
		MDS[i,:] 	= [ G[SRRS[i]][0][m] for m in MOTIFS ]
		CTS[i,:] 	= [ G[SRRS[i]][1][m] for m in MOTIFS ]
	return MDS,CTS, SRRS, MOTIFS,df
def make_distance_matrix(MDS, CTS, SRRS, MOTIFS, df,OUT):
	ST 						= dict([(y,x) for x,y in enumerate(SRRS) ])
	FHW 						= open(OUT, "w")
	FHW.write(",".join(SRRS) + "\n")
	FHW.write(",".join(MOTIFS) + "\n")
	n,G 						= MDS.shape[0], dict([ (SRRi, dict([(SRRj, [list(),list()] ) for SRRj in SRRS]) ) for SRRi in SRRS])
	for i in range(n):
		print i
		for j in range(i+1,n):
			diffs 					= MDS[j,:]-MDS[i,:]
			m 							= np.mean(diffs)
			u,v 						= SRRS[i], SRRS[j]
			for k in range(MDS.shape[1]):
				if abs(diffs[k]) > 0.1:
					d 	= compute_pvalue(MDS[i,k],MDS[j,k],
												CTS[i,k],CTS[j,k], mean=m,
												threshold=pow(10,-6), diff_threshold=0.1)
					'''
						if d is zero than its ia significant shift in favor of i
							0 -> i
						if d is one than its ia significant shift in favor of j
							1 -> j
					'''
					if d is not None:	
						G[u][v][d].append(k)
						G[v][u][1-d].append(k)
	for u in SRRS:
		FHW.write("\t".join([",".join(map(str, G[u][v][0])) + "|" + ",".join(map(str, G[u][v][1]))  for v in  SRRS]) + "\n")
	FHW.close()
def load_dist_table(FILE):
	FH 		= open(FILE, "r")
	lines 	= FH.readlines()
	SRRS	 	= lines[0].strip("\n").split(",")
	MOTIFS 	= dict([(u,v) for u,v in enumerate(lines[1].strip("\n").split(","))])
	G 			= list()
	for l in lines[2:]:

		line_array 	= [ dict([ (MOTIFS[int(u)],int(sign))  for sign,y in enumerate(x.split('|')) for u in y.split(",") if u ]) for x in l.strip('\n').split("\t")]
		G.append(line_array)
	return G, SRRS, MOTIFS

def draw(OUT, df, ax=None):


	motif 			= "HO_NANOG_HUMAN.H10MO.A"
	motif 			= "HO_IRF2_HUMAN.H10MO.C"
	motif 			= "HO_GATA1_HUMAN.H10MO.S"
	reset 			= False

	F 					= plt.figure(facecolor="white",figsize=(15,10))
	ax_matrix 		= F.add_axes([0.1,0.2,0.45,0.7])
	ax_enrichment 	= F.add_axes([0.6,0.2,0.3,0.7])


	DIR 				= "../files/spec_motifs_across_all_profile/"
	if motif +"_matrix" + ".csv"  not in os.listdir(DIR) or reset:
		print "Generating", motif
		G,SRRS, MOTIFS 					= load_dist_table(OUT)
		FHW 				= open(DIR+motif + "_matrix" + ".csv", "w")
		FHW2 				= open(DIR+motif + "_cti" + ".csv", "w")
		FHW3 				= open(DIR+motif + "_expi" + ".csv", "w")
		cts 				= dict([(x,0) for x in df.tissue])
		experiments 	= dict([(x,0) for x in df.treatment_code])
		A 					= np.zeros((len(G), len(G)))
		for i in range(A.shape[0]):
			cti,expi,keyword 	= df.loc[SRRS[i]].tissue,df.loc[SRRS[i]].treatment_code,df.loc[SRRS[i]].keyword
			for j in range(A.shape[0]):
				A[i,j] 			= int(motif in G[i][j]   )
			cts[cti] 				+=np.sum(A[i,:])
			experiments[expi] 	+=np.sum(A[i,:])
		FHW.write(",".join(SRRS) + "\n")
		FHW.write("\n".join([",".join( map(str, A[i,:]))  for i in range(A.shape[0])]))
		FHW2.write("cell_type,diffs\n")
		FHW2.write('\n'.join([x+ "," + str(y) for x,y in zip(cts.keys(), cts.values())])+"\n")
		FHW3.write("experiment_type,diffs\n")
		FHW3.write('\n'.join([x+ "," + str(y) for x,y in zip(experiments.keys(), experiments.values())]))
		FHW.close();FHW2.close();FHW3.close()


	df 				= df[df.organism=="human"]
	ct_counts 		= df.tissue.value_counts()
	exp_counts 		= df.treatment_code.value_counts()

	df1 				= pd.read_csv(DIR+motif+ "_matrix" + ".csv")
	df2 				= pd.read_csv(DIR+motif+ "_cti" + ".csv" )
	df3 				= pd.read_csv(DIR+motif+ "_expi" + ".csv" )
	df2 				= df2.set_index(df2.cell_type)

	A 					= df1.as_matrix()
	cts 				= list(ct_counts.index)
	cts.sort()
	xy 				= zip([df2.loc[ct].diffs/ct_counts[ct]   for ct in cts] ,range(len(ct_counts)) )
	#xy 				= zip([df2.loc[ct].diffs    for ct in cts] ,range(len(ct_counts)) )
	for i,(dff,ct) in enumerate(xy):
		ax_enrichment.plot([i,i], [0,dff],color="blue",lw=2,alpha=0.5)
		ax_enrichment.scatter([i], [dff],color="green",s=40)
	ax_enrichment.set_xlim(-1,len(ct_counts))

	ax_enrichment.set_xticks(range(len(ct_counts)))
	ax_enrichment.set_xticklabels([ cts[ct].replace("_", " ") for dff, ct in xy],rotation=90,fontsize=20)

	
	mask 		= np.tri(A.shape[0], k=0)
	A 			= np.ma.array(A, mask=mask.T)
	A 			=  A[::-1,:]
	ax_matrix.imshow(A,vmin=0, vmax=1,cmap=mpl.cm.Blues,interpolation="nearest", aspect="auto")
	despine(ax_matrix,left=True,bottom=True)
	despine(ax_enrichment)
	plt.savefig("../svg_final/M5B_triangle_" + motif + ".svg")
	plt.show()


def main():
	regenerate 	= False
	OUT 			= "../Tables/STable_to_5_pairwise_differences_all.csv"
	df 			= pd.read_csv("../files/conditions_table_joey.csv")
	df 			= df.set_index(df.SRAnumber)

	if regenerate:
		MDS_DIR 							= "/Users/joazofeifa/Lab/new_motif_distances/human/300/"
		MDS,CTS, SRRS, MOTIFS,DF 	= load_MDS_files(MDS_DIR,df)
		make_distance_matrix(MDS, CTS, SRRS, MOTIFS, df,OUT)
	draw(OUT, df, ax=None)
if __name__ == "__main__":
	main()
























