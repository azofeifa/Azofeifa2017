import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats.kde import gaussian_kde
'''
Required Files

1) processed_overlapped_files_with_FCS.tsv
'''
def load(FILE):
	T,t 	= {},0
	X,Y 	= None, None
	with open(FILE) as FH:
		for line in FH:
			if line[0]=="#":
				if X is not None:
					T[TF].append((X,Y,ct1, ct2))
				TF 	= line[1:].strip("\n").split(",")[0]
				ct1, ct2 	= line[1:].strip("\n").split(",")[1:]
				if TF not in T:
					T[TF] 	= list()
				X 	= np.zeros((32,))
				Y 	= [list() for i in range(32)]
				t+=1
			else:
				line_array 	= line.strip("\n").split("\t")
				pos 		= int(line_array[0],2)
				X[pos] 		= float(line_array[1])
				Y[pos] 		= [float(x) for x in line_array[2].split(",") if x]
	return T
def main():
	FILE 			= "../files/processed_overlapped_files_with_FCS.tsv"
	T 				= load(FILE)
	FCS 			= list()
	C 				= 1
	EXAMPLE 		= None
	BEST 			= -np.inf
	EXAMPLE2 	= None
	BEST2 		= -np.inf
	for TF in T:
		for X,Y,ct1, ct2 in T[TF]:
			pos1 		= int("11100" ,2)			
			pos1b 	= int("11101" ,2)			

			pos2 		= int("11010" ,2)
			pos2b 	= int("11011" ,2)

			pos3 		= int("11000" ,2)
			pos4 		= int("11110" ,2)
			pos5 		= int("11111" ,2)
			pos6 		= int("11001" ,2)

			if len(Y[pos1]) > C  and len(Y[pos2]) > C and len(Y[pos3])> C and len(Y[pos4])> C:
				NORM 		= np.mean(Y[pos1] + Y[pos2] + Y[pos3] + Y[pos4] + Y[pos5] + Y[pos6] + Y[pos1b] + Y[pos2b] )

				med1 		= np.mean(Y[pos1])
				med2 		= np.mean(Y[pos2])
				
				med3 		= np.mean(Y[pos3]+ Y[pos6])
				med4 		= np.mean(Y[pos4]+ Y[pos5])

				ARG 		= np.median(Y[pos1]-med3) -np.median(Y[pos2]-med4)

				if ARG > BEST and len(Y[pos1]) > 50:
					BEST 		= ARG
					print TF,ct1,ct2
					EXAMPLE 	= np.random.normal(0,2,len(Y[pos3])),np.random.normal(0,2,len(Y[pos3])),Y[pos1]-med3, Y[pos2b]+med3
				elif ARG > BEST2 and len(Y[pos1])>10:
					BEST2 	= ARG
					EXAMPLE2 = np.random.normal(0,2,len(Y[pos3])),np.random.normal(0,2,len(Y[pos3])),Y[pos1]-med3, Y[pos2b]+med3

				FCS.append(( med4-NORM,med3-NORM, med1-med3  , med2-med3 ))


	TS 		= ("CT_1: BTE, CT_2: No BTE", "CT_1: No BTE, CT_2: BTE","CT_1: BTE, CT_2: BTE","CT_1: No BTE, CT_2: No BTE" )

	sns.set_style("white")
	sns.set_context("notebook", font_scale=1.7, rc={"lines.linewidth": 3.5})
	F 		= plt.figure(figsize=(9,10))
	axes 	= F.add_subplot(4,2, 2), F.add_subplot(4,2, 4),F.add_subplot(4,2, 6),F.add_subplot(4,2, 8)  
	axesl = F.add_subplot(4,2, 1), F.add_subplot(4,2, 3),F.add_subplot(4,2, 5),F.add_subplot(4,2, 7)  
	x 		= np.linspace(-2,2,100)
	
	vs 				=  sns.color_palette("colorblind", n_colors=12)
	cs 				= (vs[5], vs[8])


	for i,ax in enumerate(axes):
		samp 	=[ fc[i] for fc in FCS ]
		my_pdf = gaussian_kde(samp)
		ax.hist(samp,range=(-2,2),bins=50,edgecolor="white",normed=1,alpha=0.6)
		ax.plot(x,my_pdf(x),lw=0.5, ls="--", color="green")
		bp 	= axesl[i].boxplot((EXAMPLE[i],),vert=False,patch_artist=True,notch=True)
		for u,box in enumerate(bp['boxes']):
			# change outline color
			box.set( color=cs[u], linewidth=0.5)
			# change fill color
			box.set( facecolor = cs[u])
		for u,whisker in enumerate(bp['whiskers']):
			if u < 2:
				color 	= cs[0]
			else:
				color 	= cs[1]
			whisker.set(color=color, linewidth=0.5)
		for u,cap in enumerate(bp['caps']):
			if u < 2:
				color 	= cs[0]
			else:
				color 	= cs[1]

			cap.set(color=color, linewidth=0.5)
		axesl[i].set_xlim(-10,10)
		axesl[i].set_ylim(0.7,1.30)
		axesl[i].set_yticklabels(["NR2F2", "MAFF"])
		axesl[i].legend(loc="best",fontsize=7)
		ax.set_xlim(-2,2)
		if i < 3:
			ax.set_xticks([])
			axesl[i].set_xticks([])	
		else:
			ax.set_xticks([int(x) for x in np.linspace(-2,2,5)])

	sns.despine()
	plt.tight_layout()
	plt.savefig("../svg_final/M2E_nearest_hists.svg")
	plt.show()

main()