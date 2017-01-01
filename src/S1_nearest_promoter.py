'''
	Requires:
		1) STable_to_1_promoter_distances.csv  (Table from ST_S1_nearest_promoter.py )
	Outputs:
		1) Figure with GMM and Nearest to Promoter distribution
		2) Also generates a table of nearest distances
'''

import scipy.spatial as ss
import numpy as np
import matplotlib.pyplot as plt
import time,math
from sklearn import mixture
import scipy.stats as SS



def load_nearest(FILE):
	header,FH 	= True,open(FILE)
	D 		= list()
	for line in FH:
		if not header:
			line_array 	= line.strip('\n').split(',')

			D.append(math.log(float(line_array[-1])+1,10))
		else:
			header 		= False
	X 			= np.zeros((len(D),2))
	X[:,0] 	= D
	return X
def density_estimation(D,K=2):
	clf = mixture.GaussianMixture(n_components=K, covariance_type='full')
	clf.fit(D)
	norms	= [SS.norm(clf.means_[i,0],math.sqrt(clf.covariances_[i,0,0] ))  for i in range(K)]


	F  	= plt.figure(facecolor="white",figsize=(15,7),tight_layout=True)
	xs 	= np.linspace(min(D[:,0]), max(D[:,0]),1000)
	ys 	= [ sum([ n.pdf(x)*clf.weights_[i] for i,n in enumerate(norms)]) for x in xs]
	ax 	= plt.gca()
	for i,n in enumerate(norms):
		ys 	= [n.pdf(x)*clf.weights_[i] for x in xs]
		ax.plot(xs,ys,label=r'$\mu_' + str(i+1) + '=10^{' + str(clf.means_[i,0])[:5] + "}$ bp" ,lw=2 )
	ax.hist(D[:,0],bins=50,normed=1,alpha=0.7,edgecolor="white",color="blue",label="N=31,613")
	ax.legend(loc="best",fontsize=20)
	ax.set_xlabel("Distance (log) to\nNearest Promoter (RefSeq)",fontsize=20)
	ax.set_ylabel("Normalized Frequency",fontsize=20)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)	
	plt.savefig("../svg_final/S1_NearestPromoter.svg")
	plt.show()
def main():
	TABLE 	= "../Tables/STable_to_1_promoter_distances.csv"
	D 			= load_nearest(TABLE)
	density_estimation(D)
	
if __name__ == "__main__":
	main()

