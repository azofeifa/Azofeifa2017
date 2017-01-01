import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.spatial.distance import pdist,squareform
import scipy.cluster.hierarchy as sch
import matplotlib as mpl
from matplotlib import cm
import matplotlib
from scipy.interpolate import spline, UnivariateSpline
import math

def draw(ax1,ax2, motif="HO_GATA1_HUMAN.H10MO.A",BINS=150):
	X,SUM 	= list(),list()
	DIR 		= "/Users/joazofeifa/Lab/new_motif_distances/human/300/"
	SRRS 		= list()
	for i,FILE in enumerate(os.listdir(DIR)):
		FH 		= open(DIR+FILE, "r")
		collect 	= False
		SRRS.append(FILE.split('_')[0])
		for line in FH:
			if  "Emp" in line and collect:
				break
			elif "Binned" in line:
				collect 	= True
			elif collect:
				line_array 	= line.strip("\n").split("\t")
				if line_array[0]==motif:
					x 				= map(float,line_array[1].split(","))
					counts,x1 		= np.histogram(np.linspace(-1500,1500,len(x)), weights=x,bins=BINS)
					sum1,sum2 		= sum([counts[j] for j,x in enumerate(x1[1:]) if abs(x) < 150  ]), sum(counts)
					if sum1 and sum2:
						p 				= sum1/sum2
						XX 					= np.zeros((BINS,2))
						XX[:,0] 			= (x1[1:] + x1[:-1])/2.
						XX[:,1] 			= counts

						X.append(counts)
						SUM.append(p)
	X 		= np.array(X)
	for i in range(X.shape[0]):
		X[i,:]/=np.max(X[i,:])
	A 		= squareform(pdist(X, "cosine"))

	Y 		= sch.linkage(A, method='ward')
	Z2 		= sch.dendrogram(Y, no_plot=True)
	idx1 		= Z2['leaves']
	X 			= X[idx1,:]
	center 	= len(idx1)/2
	pos,neg 	= [1 if SUM[i] > 0.1 else 0 for i in idx1[::-1]], [1 if SUM[i] < 0.1 else 0 for i in idx1[::-1]]
	ax2.plot([SUM[i] for i in idx1[::-1]], range(len(idx1)),lw=1,color="black")
	ax2.fill_betweenx( range(len(idx1)),[SUM[i] for i in idx1[::-1]], x2=0.1, where=pos,alpha=0.3,color="blue"  )
	ax2.fill_betweenx( range(len(idx1)),[SUM[i] for i in idx1[::-1]], x2=0.1, where=neg,alpha=0.3,color="red"  )
	ax2.plot([0.1,0.1], [0,len(X) ],ls="--",color="green")
	ax2.set_ylim(0,len(X))
	ax1.imshow(X, cmap=cm.GnBu,interpolation='nearest', aspect='auto',vmin=np.min(X), vmax=np.max(X))
	ax1.set_title(motif.split("_")[1])


	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.spines['left'].set_visible(False)
	for ax in (ax1,ax2):
		ax.set_xticks([]);ax.set_yticks([])
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

def main():
	F 	= plt.figure(figsize=(15,6),facecolor="white")
	ax1 = F.add_axes([0.1,0.1,0.15,0.8])
	ax2 = F.add_axes([0.26,0.1,0.03,0.8])

	draw(ax1,ax2, motif="HO_NRF1_HUMAN.H10MO.A",BINS=150)


	ax1 = F.add_axes([0.3,0.1,0.15,0.8])
	ax2 = F.add_axes([0.46,0.1,0.03,0.8])

	draw(ax1,ax2, motif="HO_JUND_HUMAN.H10MO.A",BINS=150)


	ax1 = F.add_axes([0.5,0.1,0.15,0.8])
	ax2 = F.add_axes([0.66,0.1,0.03,0.8])

	draw(ax1,ax2, motif="HO_REL_HUMAN.H10MO.C",BINS=150)

	ax1 = F.add_axes([0.7,0.1,0.15,0.8])
	ax2 = F.add_axes([0.86,0.1,0.03,0.8])

	draw(ax1,ax2, motif="HO_NANOG_HUMAN.H10MO.A",BINS=150)
	plt.savefig("../svg_final/M3E_barcodes_all.svg")
	plt.show()

main()
