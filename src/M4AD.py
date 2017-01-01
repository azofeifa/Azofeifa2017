'''
	Requires:
		1) 4 csv tables from S_Table_4AD_differential_example.py scripy
			to be in the ../Tables/ directory
		2) The the motif displacement package!
	Outputs
		1) an svg figure showing the "MA Plot" of the four examples

'''


import sys
sys.path.append("/Users/joazofeifa/Lab/motif_displacement_analysis/MDAP/")
import differential_MD as dMD
import display,matplotlib.pyplot as plt


def main():
	
	F 			= plt.figure(facecolor="white", figsize=(17,5),tight_layout=True) 
	ax1 		= F.add_subplot(1,4,1)
	ax2 		= F.add_subplot(1,4,2)
	ax3 		= F.add_subplot(1,4,3)
	ax4 		= F.add_subplot(1,4,4)

	out 		= "../Tables/STable_to_4A_Allen2014_Differential.csv"
	display.show(out, ax=ax1,FDR=pow(10,-10))
	ax1.set_title('Allen(2014)')
	ax1.set_ylabel('Nutlin-DMSO')


	out 		= "../Tables/STable_to_4B_Luo2014_Differential.csv"
	display.show(out, ax=ax2,FDR=pow(10,-10))
	ax2.set_title('Luo(2014)')
	ax2.set_ylabel('TNFalpha-DMSO')


	
	out 		= "../Tables/STable_to_4C_Hah2013_Differential.csv"
	display.show(out, ax=ax3,FDR=pow(10,-10))
	ax3.set_title('Hah(2013)')
	ax3.set_ylabel('Estradiol-DMSO')

	out 		= "../Tables/STable_to_4D_Puc2015_Differential.csv"
	display.show(out, ax=ax4,FDR=pow(10,-5))
	ax4.set_title('Puc(2015)')
	ax4.set_ylabel('DHT-DMSO')



	for ax in (ax1,ax2,ax3,ax4):
		ax.set_ylim(-0.3,0.3)
		ax.set_xlabel("")
	plt.savefig("../svg_final/M4AD_differential_examples.svg")
	plt.show()
	
	pass
main()