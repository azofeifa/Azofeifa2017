import matplotlib.pyplot as plt
def main():
	FILE 	= "../files/SRR1552480_MDS.tsv"
	L,TYPE 		= list(),0
	X,Y,colors 	= list(),list(),list()
	for line in open(FILE):
		if "Binned" in line:
			break
		elif line[0]!="#":
			line_array 				= line.strip("\n").split("\t")
			(mds,pv1,pv2),null 	= [ map(float, x.split(",")) for x in line_array[2:5]],map(float,line_array[-2].split(","))
			X.append(null[TYPE])
			Y.append(mds[TYPE])

			if pv1[TYPE] < -1 and pv2[TYPE] < -1:
				color = "red"
			elif pv1[TYPE] > -0.00001 and pv2[TYPE] > -0.00001:
				color = "green"
			else:
				color = "blue"
			colors.append(color)

	F 		= plt.figure(facecolor="white")
	ax 	= plt.gca()
	ax.scatter(X,Y,c=colors,edgecolor="")
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlabel("K562 Observed MDS",fontsize=20)
	ax.set_ylabel("Simulated Null MDS",fontsize=20)
	plt.savefig("../svg_final/M3B_obs_vs_null.svg")
	plt.show()
main()