from math import exp,log
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
class uniform:
	def __init__(self, a,b, w):
		self.a 	= a
		self.b 	= b
		self.l 	= abs(b-a)
		self.w 	= w
		self.type 	= "uniform"
	def __str__(self):
		return "uniform (" + str(self.a) + "," + str(self.b) + "," + str(self.w)+ ")"

	def pdf(self, x):
		return int(self.a<=x<=self.b)*(self.w/self.l)
class laplace:
	def __init__(self, mu, b, w ):
		self.mu 	= float(mu)
		self.b 		= float(b)
		self.w 		= w
		self.type 	= "laplace"
	def __str__(self):
		return "laplace (" + str(self.mu) + "," + str(self.b) + "," + str(self.w) + ")"
	def pdf(self, x):
		return (self.w / (2.0*self.b))*exp(-abs(x-self.mu)/self.b)
class EM:
	def __init__(self, k=1, T=300, ct=pow(10,-4)):
		self.k 	= k
		self.T 	= T
		self.ct = ct
		self.rvs=list()
		self.K 	= {}
		self.N 	= 0.0
	def run_all_three(self, X):
		for k in range(3):
			rvs, ll 	= self.fit(X,k)
			self.K[k] 	= (ll, rvs)
	def get_bic(self):
		val 	= min([ (-2*self.K[k][0] + log(self.N)*pow(k+1,1),k) for k in self.K])
		return val[1]
	def fit(self,X,k):
		a,b 	= X[0,0], X[-1,0]
		self.N 	= sum(X[:,1])
		if k == 0:
			rvs 	= [uniform(a,b,0.5)]
		elif k == 1:
			rvs 	= [uniform(a,b,0.5) , laplace(0,150,0.5)]
		elif k == 2:
			rvs 	= [uniform(a,b,0.33) , laplace(-30,150,0.33),laplace(30,150,0.33)]
		t 			= 0
		self.rvs 	= rvs
		prevll 		= -np.inf
		while t < self.T and k > 0:
			EX 	= np.zeros((k+1,))
			mus = np.zeros((k,))
			bs 	= np.zeros((k,))
			ll 	= 0.0

			for i in range(X.shape[0]):
				norms 	= np.zeros((k+1,))
				for j,rv in enumerate(rvs):
					norms[j] 	= rv.pdf(X[i,0])
				ll+=log(sum(norms))*X[i,1]
				
				norms[:]/=sum(norms)
				for j,rv in enumerate(rvs):
					if rv.type == "laplace":
					 	if k > 1:
							mus[j-1]+=norms[j]*X[i,0]*X[i,1]
						bs[j-1]+=norms[j]*abs(X[i,0]-rv.mu)*X[i,1]
				for j in range(k+1):
					EX[j]+=norms[j]*X[i,1]
			for j,rv in enumerate(rvs):
				rv.w 	= EX[j] / sum(EX)
				if rv.type=="laplace":# and k ==1:
					rv.b 	= bs[j-1] / EX[j]
			if k > 1:
				b 	= (bs[0] / EX[1]) + (bs[1] / EX[2])
				b 	/= 2.0
				v 	= abs(mus[0]/EX[1]) + abs(mus[1]/EX[2])
				v 	/= 2.0
				rvs[1].mu 	=-v #mus[0]/EX[1]
				rvs[2].mu 	= v #mus[1]/EX[2]
				rvs[1].b, rvs[2].b 	= b,b

				W 	= (rvs[1].w + rvs[2].w )/2.0
				rvs[1].w 	= W
				rvs[2].w 	= W
				S 			= sum([rv.w for rv in rvs])
				for rv in rvs:
					rv.w/=S 	
			if abs(ll - prevll) <self.ct:
				break

			prevll 	= ll
			t+=1
		self.rvs 	= rvs
		ll 			= self.get_ll(X)
		return rvs,ll
	def draw(self,X):
		ax 	= plt.gca()
		ax.hist(X[:,0], weights=X[:,1],normed=1.0,edgecolor="white",bins=X.shape[0],alpha=0.75)
		xs 	= np.linspace(X[0,0],X[-1,0], 200)
		ax.plot(xs, map(self.pdf, xs))
		plt.show()
	def pdf(self, x):
		return sum([rv.pdf(x) for rv in self.rvs])
	def get_ll(self, X):


		return sum([log(self.pdf(X[i,0]))*X[i,1]  for i in range(X.shape[0])])
