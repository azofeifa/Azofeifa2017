'''
	Required:
		1) all_processed.tsv
		2) ENCFF000OYM.bigWig
		3) DHS1.bw
		4) H3K4me1.bw
	Outputs: (2 FIGURES)
		1) A figure showing an example super enhacner with GRO-seq, DHS and H3K27ac etc.
		2) meta profile of mark signal across the bidirectional
'''
import pyBigWig,numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from scipy.stats.kde import gaussian_kde

def load_file_meta_mark(FILE):
	G 		= {}
	FH 	= open(FILE, 'r')
	for line in FH:
		if line[0]!="#":
			line_array 	= line.strip("\n").split("\t")
			MARK 			= line_array[1]
			TSS 			= map(float, line_array[2].split(","))
			NON 			= map(float, line_array[3].split(","))
			if not MARK:
				MARK 		= "DHS1"
			if MARK not in G:
				G[MARK] 	= [list(),list()]
			G[MARK][0].append(TSS)
			G[MARK][1].append(NON)
	FH.close()
	return G
def main():
	FILE 		= "../files/all_processed.tsv"
	G 			= load_file_meta_mark(FILE)
	BINS 		= 500
	start 	= 8221176
	stop 		= 8260207 
	GENES 	= list()
	origins 	= [ 0.1*BINS , 0.16*BINS , 0.24*BINS,0.47*BINS,0.69*BINS,0.78*BINS,0.84*BINS  ]
	origins 	= [o+30 for o in origins]

	N 			= 4

	DIR 		= "../files/"
	H3K27ac 	= pyBigWig.open(DIR+"ENCFF000OYM.bigWig")
	DHS1   	= pyBigWig.open(DIR+"DHS1.bw")
	H3K4me1 	= pyBigWig.open(DIR+"H3K4me1.bw")
	GROf 		= pyBigWig.open("../files/DMSO2_3.pos.bw")
	GROr 		= pyBigWig.open("../files/DMSO2_3.neg.bw")


	DHS1_counts   		= DHS1.stats("chr1", int(start), int(stop),nBins=BINS)
	H3K27ac_counts   	= H3K27ac.stats("chr1", int(start), int(stop),nBins=BINS)
	H3K4me1_counts   	= H3K4me1.stats("chr1", int(start), int(stop),nBins=BINS)
	GROf_counts 		= GROf.stats("chr1", int(start), int(stop),nBins=BINS)
	GROr_counts 		= GROr.stats("chr1", int(start), int(stop),nBins=BINS)

	F 						= plt.figure(facecolor="white",figsize=(16,7))
	IDS 					= "DHS1", "H3K27ac\n(Chromatin Mod.)", "H3K4me1\n(Chromatin Mod.)", "GRO-seq\n(Transcription)","GRO-seq\n(Transcription)"
	counts 				= DHS1_counts,H3K27ac_counts,H3K4me1_counts,GROf_counts,GROr_counts

	for i,ID in enumerate(IDS):
		if i < 4:
			ax 				= F.add_subplot(N,1,i+1)

		ax.set_ylabel(IDS[i],rotation=0,fontsize=20)
		edges 			= range(BINS)
		UNROLL 			= list()
		for k,x in enumerate(counts[i]):
			if x:
				for j in range(int(x)):
					UNROLL.append(edges[k])
		my_pdf = gaussian_kde(UNROLL,bw_method=0.0047)
		x = np.linspace(edges[0],edges[-1],1000)
		if i == 4:
			ax.fill_between(x,-my_pdf(x), np.zeros((len(x), )),color="red",alpha=0.4)
		else:
			ax.fill_between(x, np.zeros((len(x), )), my_pdf(x),color="blue",alpha=0.4)
		ax.set_xlim(-30,BINS)	

		ax.set_xticks([])
		ax.set_yticks([])
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.spines['left'].set_visible(False)
	if N > 6:
		ax 				= F.add_subplot(N,1,N)
		ax.set_ylabel("RNAP Origin\n(Tfit)", rotation=0,fontsize=20)
		ax.scatter(origins,np.zeros((len(origins))),s=50 )
		ax.set_ylim(-1,1)
		ax.set_xlim(-30,BINS)
		ax.set_xticks([])
		ax.set_yticks([])
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.spines['left'].set_visible(False)


	# plt.show()
	plt.savefig("../svg_final/M1A_IGV_Chromatin_Mark_Example.svg")

	F 		= plt.figure(figsize=(16,10),facecolor="white")
	step 	= 0.9/len(G)
	MM 	= 0
	for mark in G.keys():
		counts 	= np.array(G[mark][1])[0,:]
		NORM 		= np.sum(counts)
		counts 	/=NORM
		MM 		= max(MM, np.max(counts))
		

	MARKS 		= ("DHS1", "H2AFZ", "H3K4me2", "H3K9ac", "H3K27ac", "H3K4me1", "H3K4me3", "H3K9me1", "H3K36me3", "H3K79me2", "H4K20me1", "H3K9me3", "H3K27me3" )[::-1]
	for i,mark in enumerate(MARKS):
		ax 		= F.add_axes([0.1,0.05+step*i, 0.7,step/2.])
		counts 	= np.array(G[mark][1])[0,:]
		NORM 		= np.sum(counts)
		if mark!="DHS1":
			counts 	/=NORM
			norm 		= mpl.colors.Normalize(vmin=0, vmax=MM)
		else:
			norm 		= mpl.colors.Normalize(vmin=min(counts), vmax=max(counts))
		cmap 		= cm.GnBu
		M			= cm.ScalarMappable(norm=norm, cmap=cmap)

		colors1 	= map(M.to_rgba, counts)
		ax.bar(np.linspace(-1500,1500,len(counts)),np.ones((len(counts), )), color=colors1, edgecolor=colors1,width=3000.0/len(counts) )
		ax.set_yticks([])
		ax.set_ylabel(mark, rotation=0,fontsize=25)

		ax.tick_params(labelsize=20)
		if i > 0:
			ax.set_xticks([])
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.spines['left'].set_visible(False)

		ax.set_xlim(-1500,1500)
	ax 		= F.add_axes([0.9,0.05 , 0.02,0.9])
	norm 		= mpl.colors.Normalize(vmin=0, vmax=1)
	cmap 		= cm.GnBu
	M			= cm.ScalarMappable(norm=norm, cmap=cmap)
	cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
	
	ax.tick_params(labelsize=20)
	plt.savefig("../svg_final/M1B_HeatMaps_Marks.svg")

	plt.show()
if __name__ =="__main__":
	main()