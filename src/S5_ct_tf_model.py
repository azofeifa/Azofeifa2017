import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns,math,numpy as np
def display(df):

	dfF 	= df[(df.bin_n*df.bin_p > 0) & (df.obs > 0) ]
	sns.set_style("white")
	F 	= plt.figure(figsize=(15,6),tight_layout=True)
	ax1=F.add_subplot(1,2,1)
	ax1.set_xscale("log")
	ax1.set_yscale("log",basey=2)
	df1,df2 		= dfF[df.pvalue > 0.999],dfF[df.pvalue < 0.999]
	x1,y1 		= df1.bin_n*df1.bin_p,df1.obs/(df1.bin_n*df1.bin_p)
	x2,y2 		= df2.bin_n*df2.bin_p,df2.obs/(df2.bin_n*df2.bin_p)

	ax1.scatter(x1,y1,edgecolor="",color="red",label="significant TF-ct associations")
	ax1.scatter(x2,y2,edgecolor="",color="blue",alpha=0.2)

	ax1.plot([min(x2), max(x2)],[1,1],lw=2.0,color="green",label="No Change")
	ax1.set_xlabel("\n\n" + r"$N(tf)\cdot\sum_{exp\in EXP}p(ct|exp)p(exp|TF) $",fontsize=20)
	ax1.set_ylabel(r"$N(TF=tf ,CT=ct )$",fontsize=20)
	ax1.legend(loc="best")

	ax2=F.add_subplot(1,2,2)

	tfs 		= set(df.motif)
	ENTROPY 	= [ (df[df.motif==tf].ENR+1) / (sum(df[df.motif==tf].ENR+1))  for tf in tfs]

	NN 		= len(ENTROPY[7].index)
	ENTROPY 	= [ -sum([p*math.log(p,2) for p in e ]) for e in ENTROPY]
	counts,edges 	= np.histogram(ENTROPY,bins=50)
	edges 			= (edges[1:] + edges[:-1])/2.
	width 			= (edges[-1] - edges[0])/50.0
	ax2.bar(edges,counts, width=width,alpha=0.7 )


	ax2.plot([-math.log(1.0/NN, 2),-math.log(1.0/NN, 2) ],[0,max(counts)],lw=3,ls="--",color="green",label="Max Entropy\nNo Information" )
	ax2.legend(loc="best")
	ax2.set_xlabel("\n" +r'$H(CT|TF=tf)$',fontsize=20.0)
	ax2.set_ylabel("Frequency(n=641)",fontsize=20.0)
	sns.despine(ax=ax1)
	sns.despine(ax=ax2)

	plt.savefig("../svg_final/S5_null_ct.svg")
	plt.show()



df 	= pd.read_csv("../Tables/STable_to_5D_cell_tf_associations.csv")
display(df)