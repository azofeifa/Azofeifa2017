import matplotlib.pyplot as plt
import NormalLaplace as NL
import numpy as np, math
import scipy.stats as ss
from scipy.spatial.distance import pdist,squareform
import scipy.cluster.hierarchy as sch
import matplotlib as mpl
from matplotlib import cm

'''
Required Files

1) SRR1552480_MDS.tsv
'''
def pval(sum1,sum2,p=0.1):
	mu 	= p*sum2
	sigma = mu*(1.0-p)
	Z 		= (sum1 - mu)/math.sqrt(sigma)
	N 		= ss.norm(0,1)
	return N.cdf(Z)

def main(WRITE=False):
	f 						= "../files/SRR1552480_MDS.tsv"
	FH 					= open(f,'r')
	X,collect 			= list(), False
	M 						= list()
	SUM 					= list()
	PVS 					= list()
	for line in FH:
		if  "Emp" in line and collect:
			break
		elif "Binned" in line:
			collect 	= True
		elif collect:
			line_array 	= line.strip("\n").split("\t")
			x 				= map(float,line_array[1].split(","))
			counts,x1 		= np.histogram(np.linspace(-1500,1500,len(x)), weights=x,bins=200)
			sum1,sum2 		= sum([counts[j] for j,x in enumerate(x1[1:]) if abs(x) < 150  ]), sum(counts)
			p 				= sum1/sum2
			XX 					= np.zeros((200,2))
			XX[:,0] 			= (x1[1:] + x1[:-1])/2.
			XX[:,1] 			= counts

			PVS.append(pval(sum1,sum2, p=0.1) )
			X.append(counts)
			SUM.append(p)
	X 		= np.array(X)
	for i in range(X.shape[0]):
		X[i,:]/=np.max(X[i,:])

	idx1  	= [(p,i) for i,p in enumerate(SUM)]
	idx1.sort()
	idx1 		= [u[1] for u in idx1]

	X 			= X[idx1,:]
	F 			= plt.figure(facecolor="white",figsize=(6,10))
	ax1 		= F.add_axes([0.1,0.1,0.5,0.8])
	ax2 		= F.add_axes([0.75,0.1,0.2,0.8])
	center 	= len(idx1)/2
	pos,neg 	= [1 if SUM[i] > 0.1 else 0 for i in idx1[::-1]], [1 if SUM[i] < 0.1 else 0 for i in idx1[::-1]]
	ax2.plot([SUM[i] for i in idx1[::-1]], range(len(idx1)),lw=1,color="black")
	ax2.fill_betweenx( range(len(idx1)),[SUM[i] for i in idx1[::-1]], x2=0.1, where=pos,alpha=0.3,color="blue"  )
	ax2.fill_betweenx( range(len(idx1)),[SUM[i] for i in idx1[::-1]], x2=0.1, where=neg,alpha=0.3,color="red"  )
	ax2.plot([0.1,0.1], [0,641],ls="--",color="green")
	ax2.set_ylim(0,641)
	ax1.imshow(X, cmap=cm.GnBu,interpolation='nearest', aspect='auto',vmin=0, vmax=1)
	for ax in (ax1,):
		ax.set_xticks(np.linspace(0,X.shape[1], 5)  )
		ax.set_xticklabels(["-1500", "-750", "0", "750", "1500"])
	# Hide the right and top spines
	ax1.set_xticklabels(["-1500","-750","0","750","1500"],fontsize=16)
	ax1.set_yticklabels([ "", "0","100","200","300","400", "500", "600"],fontsize=26)
	ax2.set_yticks([])
	ax2.set_xticks([0,0.15,0.3,0.45])
	ax2.set_xticklabels(["0","0.15","0.3","0.45"],rotation=45)
	ax2.set_xlim(0,0.55)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.spines['left'].set_visible(False)
	for ax in (ax1,):
		# Only show ticks on the left and bottom spines
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')

		ax.spines['left'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		# Only show ticks on the left and bottom spines
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')

	plt.savefig('../svg_final/M3A_MDS_K562_all.svg')

	plt.show()

main()