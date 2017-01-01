import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn3, venn3_circles,venn2

'''
Required:

1) overlap_stats.bed
'''
def main():
	FILE 	= "../files/overlap_stats.bed"
	FH 	= open(FILE)
	D 	= list()
	ALL = {}
	G 	= {}
	P 	= list()
	T 	= None
	S 	= 3.2 * pow(10,9)
	lines 	= [line for line in FH]
	size 	= False
	
	for line in lines:
		lbls, sizes,data 	= line.strip("\n").split("\t")
		lbls 		= lbls.split(",")
		data 		= [int(d) for d in data.split(",")]
		sizes 		= [float(d) for d in sizes.split(",")]
		
		if len(lbls)==2:
			if size:
				G[lbls[0]] 		= sizes[0] /S 

				G[lbls[1]] 		= sizes[1] /S
			else:
				G[lbls[0]] 		= data[0] 
				G[lbls[1]] 		= data[1] 				
			if size:
				P.append([set(lbls), sizes[-1]/S ])
			else:
				P.append([set(lbls), data[-1] ])
				
		else:
			if size:
				T 			= sizes[-1]/S
			else:
				T 			= data[-1]

	P 		= [(x,y-T) for x,y in P]
	for lbl in G:
		for ol,per in P:
			if lbl in ol:
				G[lbl]-=per
		G[lbl]-=T


	plt.figure(figsize=(15,10))
	sns.set_style("white")
	if size:
		R 			= 4

		subsets 	= (round(float(G["Bidir"]),R), round(float(G["TF"] ),R), 
							round(float(P[0][1] ),R), round(float(G["DNase"] ),R), 
							round(float(P[1][1] ),R), round(float(P[2][1] ),R) , round(float(T ),R))
	else:
		subsets 	= ( G["Bidir"]-1000 ,  G["TF"] , 
							 P[0][1]  ,  G["DNase"] -5000  , 
							 P[1][1]+1000 ,  P[2][1]   ,  T+1000 )
		
	v = venn3(subsets=(subsets), set_labels = ("BTE\n(GRO-seq)", "TF Binding\n(ChIP)", "Open Chromatin\n(DNase)"))

	for i in range(3):	 v.set_labels[i].set_fontsize(30)
	for i in range(len(v.subset_labels)):	v.subset_labels[i].set_fontsize(25)
	for i in range(3):	v.set_labels[i].set_ha("right")
	v.set_labels[1].set_va("bottom")

	for i in range(len(v.patches)):	v.patches[i].set_alpha(0.7)

	vs 			=   sns.color_palette("colorblind", n_colors=7, desat=.25)
	j 				= 0

	v.get_patch_by_id('100').set_color("blue")
	v.get_patch_by_id('010').set_color("red")
	v.get_patch_by_id('001').set_color("green")


		
	for i in range(len(v.subset_labels)):	v.subset_labels[i].set_color("black")
	for i in range(len(v.subset_labels)):	v.subset_labels[i].set_ha("center")
	for i in range(len(v.subset_labels)):	v.subset_labels[i].set_va("center")


	c = venn3_circles(subsets=subsets,color="white")
	plt.savefig("../svg_final/M2A_Venn.svg")
	plt.show()
main()
