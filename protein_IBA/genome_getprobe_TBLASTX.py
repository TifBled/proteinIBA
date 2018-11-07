#!/usr/bin/env python


# Script was written by Jesse Breinholt (jessebreinholt@gamil.com, jbreinholt@rapid-genomics.com) for Espeland et al 2018.
#This script takes fasta file for which there is a single seqeunce per probe loci where each sequnce
#name begins with a "L" and and number followed by an "_" (ex. L1_REFTAXANAME) and a assembled genome
#or transcrptome in fasta formate with sequnces on a single line, the number of threads,  a taxa name
#you want on the resulting sequences, and flank size the number of base pairs before and after the hit
#of the probe sequnces. The script outputs two files a blast table of hits and a fasta file which
#contains sequnces that reasonably resemble the loci in the probe fasta file. Note if the genome or
#transcriptome file does not contain the nucleotide data of each sequence on a single line the output
#will not produce the correct output.


#this script has been modified for (Wolfe et 2018) to use tblastx instead of blastn to screen loci from denovo assemblies of capture data against single copy exon.
#because it used tblastx the same fragment but diffrent coding frame will make it into results you may want to collapse redundant sequnces one way to do this is with USEARCH example: usearch -cluster_fast taxaname_targets.fasta -id 0.98 -sort length -centroids taxaname_targets.fasta_cen

import os, sys


blastdb="makeblastdb -dbtype \"nucl\" -in %s -out %s"
#realblast= "blastn -task blastn -query %s -db %s -out %s -outfmt 6 -perc_identity 60 -qcov_hsp_perc 50 -num_threads %s"
realblast= "tblastx -query %s -db %s -out %s -outfmt 6 -evalue 0.0001 -qcov_hsp_perc 25 -num_threads %s"

blastdict= {}
locidict={}
seq = []
list = []
arguments = sys.argv

if len(arguments) == 1:
	sys.exit("try python genome_getprobe_BLAST.py -h for help")
try:
	hflag = arguments.index("-h")
except:
	hflag = None

if hflag:
	sys.exit("\n\n#################################################################\n\n./genome_getprobe_BLAST.py -inp probe.fasta -ing genome.fast -threads ## -tname taxaname -flanksize ## \n\n\n\t\t-inp: fasta file of probe sequences\n\t\t-ing: genome or transcriptome file with sequnces data for each sequnces on a single line\n\t\t-threads: number of threads to use in blast\n\t\t-tname: taxa name you want to put on the output sequnces\n\t\t-flanksize: number of base pairs you want inclulded before or after the probe hit\n\n#################################################################\n\n")

try:
	probeset = arguments[arguments.index("-inp")+1]
except:
	sys.exit("Error: need probe fasts")

try:
	genome = arguments[arguments.index("-ing")+1]
except:
	sys.exit("Error: need genome fasta")

try:
	taxaname = arguments[arguments.index("-tname")+1]
except:
	sys.exit("Error: need taxa name")

try:
	threads = int(arguments[arguments.index("-threads")+1])
except:
	threads = 1

try:
	flank = int(arguments[arguments.index("-flanksize")+1])
except:
	flank = 0

try:
	probeset = arguments[arguments.index("-inp")+1]
except:
	sys.exit("Error: need probe fasta")


print "making blast database"
os.system(blastdb % (genome, genome+"_db"))
print "Blasting........"
os.system(realblast % (probeset,genome+"_db",taxaname + "_out", threads))
blastfile=open(taxaname + "_out")



print "Making hit dictionary"

for line in blastfile:
	line = line.strip().split()
	if line[1] in blastdict.keys():
		list = blastdict[line[1]]
		locilist = line[0].split("_")
		list.append(locilist[0]+"_"+line[8]+"_"+line[9])
		blastdict[line[1]] = list
	else:
		list = []
		locilist = line[0].split("_")
		list.append(locilist[0]+"_"+line[8]+"_"+line[9])
		blastdict[line[1]] = list
		
print "Getting seqs with flanks"

with open(taxaname + "_targets.fasta", "w") as out:
	with open(genome) as GG:
		line=GG.readline()
		while line:
			if line[0] == ">":
				try:
					line = line.split()
					id = line[0].lstrip(">")
				except:
					id = line.strip().lstrip(">")
			if id in blastdict.keys():
				line=GG.readline()
				for info in blastdict[id]:
					a, b, c = info.split("_")
					if a in locidict:
						locidict[a] += 1
					else:
						locidict[a] = 1
					b = int(b)
					c = int(c)
					if c > b:
						if b > flank:
							start_rec = b - flank
						else:
							start_rec = 0
						end_rec = c + flank
						out.write(">"+a+"_"+taxaname+"__comp"+str(locidict[a])+"\n")
						out.write(line[start_rec:end_rec].strip()+"\n")
					elif b > c:
					#backwards
						if c > flank:
							start_rec = c - flank
						else:
							start_rec = 0
						end_rec = b + flank
						out.write(">"+a+"_"+taxaname+"__comp"+str(locidict[a])+"\n")
						out.write(line[start_rec:end_rec].strip()+"\n")
				line = GG.readline()
			else:
				line=GG.readline()
				line=GG.readline()	

sys.exit()
