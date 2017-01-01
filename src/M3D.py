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


def despine(ax):
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')

	ax.spines['left'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	# Only show ticks on the left and bottom spines
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.set_xticks([])
	ax.set_yticks([])
def main():
	SRRS 	= "SRR1105736", "SRR1596501", "SRR1648886", "SRR1552480", "SRR579299","SRR1745519"
	DIR 	= "/Users/joazofeifa/Lab/new_motif_distances/human/300/"
	M 		= {}
	VAR 	= {}
	R 		= {}
	DD 	= {}
	for i,SRR in enumerate(SRRS):
		collect 	= False
		M 			= {}
		for line in open(DIR+ SRR + "_MDS.tsv"):
			if "Binned" in line:
				break
			elif line[0]!="#":
				line_array 	= line.strip("\n").split("\t")
				motif 		= line_array[0]

				M[motif] 	=float(line_array[2].split(",")[2]) -float(line_array[-2].split(",")[0])
		DD[SRR] 	= M
		vals 		= zip(M.values(), M.keys())
		vals.sort()
		for i,(p,m) in enumerate(vals):
			if m not in R:
				R[m] 	= list()
			R[m].append(i)
	F 		= plt.figure(figsize=(7,4),facecolor="white")
	ax 	= plt.gca()
	x_plot 	= np.linspace(-1, 1, 100)
	y_plot 	= 1./(1+np.exp(-x_plot*7))
	x_plot 	= np.linspace(0, 1, 100)
	xs 		= np.linspace(-1,1,100)
	ticks 	= list()
	MS 		= list()
	F 			= ["HO_GATA1_HUMAN.H10MO.A", "HO_NANOG_HUMAN.H10MO.A", 
					"HO_ATF3_HUMAN.H10MO.A", "HO_HOMEZ_HUMAN.H10MO.D" ]
	for m in R:
		X,Y 			= list(),list()		
		for i in range(len(R[m])-1):
			slope 	= (R[m][i+1]-R[m][i])
			base 		= R[m][i]
			y_scaled = y_plot*slope + base
			X+=list(x_plot  + i)
			Y+=list(y_scaled)
		if [1 for f in F if f in m]   :
			ticks.append(Y[-1])
			MS.append(m)
			ax.plot(X,Y,lw=0.75,alpha=1.0,color="blue")
		elif np.random.uniform(0,1) < 0.4:
			ax.plot(X,Y,lw=0.5,alpha=0.05,color="green")

	norm = mpl.colors.Normalize(vmin=-0.2, vmax=0.2)
	cmap = cm.seismic
	for k in range(len(SRRS)):
		D 		= DD[SRRS[k]].values()
		D.sort()

		m 		= cm.ScalarMappable(norm=norm, cmap=cmap)
		cs 	= np.arange(0,len(D), 10)
		for i in range(1,len(cs)) :
			ax.plot([k,k], [cs[i-1], cs[i]]  ,lw=2,color=m.to_rgba(D[cs[i]]))


	ax.set_xlim(-0.1,len(SRRS)-0.9)
	despine(ax)
	ax.yaxis.tick_right()
	ax.set_yticks([])
	ax.set_xticks(range(len(SRRS)))
	ax.set_xticklabels(["Allen\n(2014)", "Andersson\n(2014)", "Puc\n(2015)","Core\n(2014)",
							"Danko\n(2013)","Estaras\n(2015)" ],fontsize=15)
	ax.set_ylabel("MD Score Rank",fontsize=20)
	plt.savefig("../svg_final/M3D_rank.svg")
	plt.show()


main()



