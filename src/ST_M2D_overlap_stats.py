'''
	Requires: (a lot)
		1)	5 directories of bed files (from ENCODE with associated )
			metadata file; named as such
				-> HCT116 
				-> K562
				-> HeLa
				-> GM
				-> MCF7
		2) 5 TFit bidirectional files from GRO-seq datasets in these cell lines
			name as such
				-> HCT116_EMG.bed 
				-> K562_EMG.bed
				-> HeLa_EMG.bed
				-> GM_EMG.bed
				-> MCF7_EMG.bed
	Outputs:
		A table of overlap statistics conditioned on Bidir presence or abscence
		-> THIS DOES NOT FILTER FOR PROMOTER ASSOCIATION
'''
import os
def load_bed_file(FILE,W=None):
	G,FH 	= {},open(FILE,'r')
	for line in FH:
		if line[0]!="#":
			chrom,start,stop 	= line.split("\t")[:3]
			if chrom not in G:
				G[chrom] 		= list()
			if W is None:
				G[chrom].append((int(start), int(stop)))
			else:
				x 		= int((int(stop) + int(start))/2.)
				G[chrom].append((x-W, x+W))

	for chrom in G:
		G[chrom].sort()
	return G
def overlap(A,B,OVER=False, LACK=False):
	L,O 	= {},{}
	for chrom in A:
		if chrom in B:
			a,b 	 	= A[chrom], B[chrom]
			j,N 		= 0, len(b)
			L[chrom], O[chrom] = list(),list()
			for start, stop in a:
				while j < N and b[j][1] < start:
					j+=1
				if j < N and b[j][0] < stop:
					O[chrom].append((start, stop))
				else:
					L[chrom].append((start, stop))

	if not OVER and not LACK:
		return sum([len(O[chrom]) for chrom in O])
	return O,L

def load_meta(FILE):
	FH 		= open(FILE)
	header 	= True
	M 			= {}
	for line in FH:
		if not header:
			line_array 			= line.split("\t")
			if "eGFP" not in line_array[15]:
				M[line_array[0]] 	= line_array[15].split("-")[0]
		else:
			header 		= False
	return M

def load_TF_dir(DIR):
	META 		= load_meta(DIR+"metadata.tsv")
	E 			= {}
	ENC 		= {}
	for f in os.listdir(DIR):
		if f.strip(".bed") in META:
			G 	= load_bed_file(DIR+f, W=750)
			m 			= META[f.strip(".bed")]
			n 	= sum(map(len, G.values()))
			if m not in E or n > sum(map(len, E[m].values())):
				E[m] 	= G
				ENC[m] 	= f.strip(".bed")
	return E,ENC
def make_table(vals, ENCS, bvals, cts, OUT):
	FHW 	= open(OUT, "w")
	FHW.write("TF Name, ENC1, ENC2, cell-type 1, cell-type 2,ct1 total,ct2 total,Overlap Total,#[00],#[01],#[10],#[11]\n")
	for i in range(len(vals)):
		for j in range(i+1, len(vals)):
			ct1,ct2 	= vals[i],vals[j]
			e1,e2 	= ENCS[i],ENCS[j]
			for TF in ct1:
				if TF in ct2:
					m1,m2 			= e1[TF], e2[TF]
					ct1_total 		= sum([len(ct1[TF][chrom]) for chrom in ct1[TF] ])
					ct2_total 		= sum([len(ct2[TF][chrom]) for chrom in ct2[TF] ])
					O,L 				= overlap(ct1[TF], ct2[TF], OVER=True)
					overlap_total 	= sum([len(O[chrom]) for chrom in O ])
					Obct1,Lbct1 	= overlap(O, bvals[i], OVER=True)
					Obct2,Lbct2 	= overlap(O, bvals[j], OVER=True)
					O11 				= overlap(Obct1, Obct2)
					O01 				= overlap(Lbct1, Obct2)
					O10 				= overlap(Obct1, Lbct2)
					O00 				= overlap(Lbct1, Lbct2)
					
					FHW.write(TF+ "," + m1 + "," + m2+ "," + 
						cts[i] + ',' + cts[j] + "," +  str(ct1_total) + "," + str(ct2_total) + "," + str(overlap_total) + "," + 
						str(O00) + "," + str(O01) + "," + str(O10) + "," + str(O11) + "\n" )

def main():




	DIR1 	= "../files/ChIP_files/HCT116/"
	DIR2 	= "../files/ChIP_files/MCF7/"
	DIR3 	= "../files/ChIP_files/K562/"
	DIR4 	= "../files/ChIP_files/GM/"
	DIR5 	= "../files/ChIP_files/HeLa/"
	DIR6  = "../files/"
	bvals = (load_bed_file(DIR6 + "HCT116_EMG.bed", W=750),
				load_bed_file(DIR6 + "MCF7_EMG.bed", W=750),
				load_bed_file(DIR6 + "K562_EMG.bed", W=750),
				load_bed_file(DIR6 + "GM12878_EMG.bed", W=750),
				load_bed_file(DIR6 + "HeLa_EMG.bed", W=750),
				)


	HCT116,H_ENC= load_TF_dir(DIR1)
	print "loaded HCT116"
	MCF7,M_ENC 	= load_TF_dir(DIR2)
	print "loaded MCF7"
	K562,K_ENC 	= load_TF_dir(DIR3)
	print "loaded K562"
	GM,G_ENC 	= load_TF_dir(DIR4)
	print "loaded GM"
	HeLa,He_ENC 	= load_TF_dir(DIR5)
	print "loaded HeLa"

	vals 	= (HCT116,MCF7,GM,K562, HeLa  )
	ENCS 	= (H_ENC,M_ENC,G_ENC,K_ENC,He_ENC )
	cts 	= ("HCT116", "MCF7", "GM12878", "K562", "HeLa")
	OUT 	= "../Tables/STable_to_2D_PreservedTFBindingOverlapStats.csv"
	make_table(vals, ENCS, bvals, cts, OUT)







