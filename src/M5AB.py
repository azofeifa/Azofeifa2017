'''
	Requires:
		1) 2 csv tables from S_Table_5AB_differential_example.py scripy
			to be in the ../Tables/ directory
		2) The motif displacement package!
	Outputs
		1) an svg figure showing the delta MD scores
'''

import sys
sys.path.append("/Users/joazofeifa/Lab/motif_displacement_analysis/MDAP/")
import differential_MD as dMD
import display,matplotlib.pyplot as plt

def main():
	F 			= plt.figure(facecolor="white", figsize=(17,7),tight_layout=True) 
	ax1 		= F.add_subplot(1,2,1)
	ax2 		= F.add_subplot(1,2,2)

	out 		= "../Tables/STable_to_5A_GM12878_vs_K562.csv"
	display.show(out, ax=ax1,FDR=pow(10,-9))
	ax1.set_title('Core(2014)')
	ax1.set_ylabel('K562-GM12878')

	out 		= "../Tables/STable_to_5B_hESC_to_ENDODERM.csv"
	display.show(out, ax=ax2,FDR=pow(10,-15))
	ax2.set_title('Wang(2014)')
	ax2.set_ylabel('Diff Stage-hESC')

	plt.savefig("../svg_final/M4EF_time_series_examples.svg")
	plt.show()



if __name__ == "__main__":
	main()


