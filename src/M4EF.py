'''
	Requires:
		1) 2 csv tables from S_Table_4EF_differential_example.py scripy
			to be in the ../Tables/ directory
		2) The motif displacement package!
	Outputs
		1) an svg figure showing the delta MD scores of the 2 examples
			across time
'''

import sys
sys.path.append("/Users/joazofeifa/Lab/motif_displacement_analysis/MDAP/")
import differential_MD as dMD
import display,matplotlib.pyplot as plt


def main():
	
	F 			= plt.figure(facecolor="white", figsize=(17,10),tight_layout=True) 
	ax1 		= F.add_subplot(1,2,1)
	ax2 		= F.add_subplot(1,2,2)

	out 		= "../Tables/STable_to_4E_Flavopiridol_Treatment.csv"
	display.show(out, ax=ax1,FDR=pow(10,-7))
	ax1.set_title('Jonkers(2014)')
	ax1.set_ylabel('TNFalpha-DMSO')
	out 		= "../Tables/STable_to_4F_KLA_Treatment.csv"
	display.show(out, ax=ax2,FDR=pow(10,-1))
	ax2.set_title('Kaikkonen(2013)')
	ax2.set_ylabel('Nutlin-DMSO')
	ax2.set_ylim(-0.1,0.1)

	plt.savefig("../svg_final/M4EF_time_series_examples.svg")
	plt.show()
main()