import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

'''
Required Files

1) TF_mark_jonathan_results.txt
'''

def main():
	results_file 	= "../files/TF_mark_jonathan_results.txt"
	M 				= {}
	with open(results_file) as FH:
		for line in FH:
			if line[0]=="#":
				mark 	= line[1:].strip("\n")
				M[mark] = list()
			else:
				line_array 	= line.strip("\n").split("\t")
				percents 			= [float(x) for x in line_array[1:5] ]
				M[mark].append([line_array[0]] + percents)
	space 			= 3
	label_fontsize 	= 20

	average_diffs 	= [(np.mean([x[1]  for x in M[mark]]), mark)   for mark in M]
	average_diffs.sort(reverse=True)
	marks 			= [y for x,y in average_diffs]
	sns.set_style("white")
	F 				= plt.figure(figsize=(5,10))
	
	ax1 			= plt.gca()
	ax1.set_xlabel("Percent\n Mark Association\n(unannotated)",fontsize=label_fontsize)
	start 			= 0
	sns.despine(left=True,bottom=True)
	vs 				=   sns.color_palette("colorblind", n_colors=12, desat=.5)
	cs 	= (vs[1], vs[0])
	for i,m in enumerate(marks[::-1]):
		TFS 		= [x[0]*100 for x in M[m]]
		tB 			= [x[1]*100 for x in M[m]]
		eB 			= [x[2]*100 for x in M[m]]
		tN 			= [x[3]*100 for x in M[m]]
		eN 			= [x[4]*100 for x in M[m]]
		bp1 		= ax1.boxplot( [tN, tB], positions =[start, start + 1] ,widths=0.6,notch=True, vert=False,patch_artist=True,sym='')
		
		for bp in (bp1, ):
			for i,box in enumerate(bp['boxes']):
				box.set( color=cs[i], linewidth=2)
				box.set( facecolor = cs[i])
			for i,whisker in enumerate(bp['whiskers']):
				if i < 2:
					color 	= cs[0]
				else:
					color 	= cs[1]
				whisker.set(color=color, linewidth=2)
			for i,cap in enumerate(bp['caps']):
				if i < 2:
					color 	= cs[0]
				else:
					color 	= cs[1]
				cap.set(color=color, linewidth=2)
		start+=space
	ax1.plot([-10,-10],[-10,-10], lw=10,label="ChIP site\nwith BTE",color=cs[0])
	ax1.plot([-10,-10],[-10,-10], lw=10,label="ChIP site\nlacking BTE",color=cs[1])
	ax1.set_yticks(np.arange(0.5,start, space))
	ax1.set_yticklabels([x.strip("-human") for x in marks[::-1]],rotation=0,fontsize=label_fontsize)
	for ax in (ax1, ):
		ax.set_ylim(-space,start)
		ax.set_xlim(0,110)
		ax.set_xticks([0,25,50,75, 100])
		ax.set_xticklabels(["0","25","50","75", "100"], fontsize=label_fontsize,rotation=0)
	plt.tight_layout()
	plt.savefig('../svg_final/M2C_boxplot_marks.svg')
	plt.show()

main()