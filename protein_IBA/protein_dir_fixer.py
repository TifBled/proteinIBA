#!/usr/bin/env python


#Script written by Jesse Breinholt (jessebreinholt@gmail.com)
#The script requires Biopython
#This script takes the output of protein_IBA (*_finaltable.table and *_finalseqs.fasta) and put sequences in the same direction of the protein reference and trims to blast hit from the reference ie. probe region.
#example: python protein_dir_fixer.py -t taxa_finaltable.table -f taxa_finalseqs.fasta
#note: taxa_finalseqs.fasta must be in single line fasta formate which means dna sequence data is on a single line per sequences 

import sys
import os
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.Alphabet import generic_rna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

arguments = sys.argv





if len(arguments) == 1:
	arguments.append("-h")
try:
	hflag = arguments.index("-h")
except:
	hflag = None

if hflag:
	sys.exit("\n\n#################################################################\n\n\t\tprotein direction fixer: protein_dir_fixer.py\
			 \n\n\n\t-t\ttable from protein_IBA\n\t-f\tfasta file from from protein_IBA\n\tNote:fasta file must be a single line fasta file and file name must only have a single \".\" in name\n\texample: python protein_dir_fixer.py -t taxa_finaltable.table -f taxa_finalseqs.fasta")

try:
	tablefile = arguments[arguments.index("-t")+1]
except:
	sys.exit("Error: need table from protein_IBA")

try:
	fasta_file = arguments[arguments.index("-f")+1]
	lname, crap1 = fasta_file.split(".")
except:
	sys.exit("Error: need fast file from from protein_IBA")

lname, crap1 = fasta_file.split(".")





#set up list, set and counts
coordL=list()
rctaxa=set([])
count=0



#parse fasta file get list of reverse complemented sequences
with open(tablefile, "r") as table:
	makein=table.readlines()
	for i in makein:
		query, target, id, length, mismatch, gaps, qs, qend, ts, tend, evalue, bit, frame, crap1, crap2=i.split()
#Use frame (-,+) to make a set of sequences that need to be reverse complemented
		if "-" in frame:
			coordL.append(query+"|"+str(int(qs)-1)+"|"+qend)
			rctaxa.add(query)
			count +=1		
		if "-" not in frame:
			coordL.append(query+"|"+str(int(qs)-1)+"|"+qend)
			pass

#parse fasta file reverse complemente if in set and  trim to probe region and translate to amino acid
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
end = False

#disabled translate as often the best hit is a short segment of the amino acid sequence enough for ortholog filtering but not ideal to use for final protein sequences. To enable feature comment out line 80 and uncomment line 79, 94,95,97,103,105,107
#with open(lname + "_rightdir.fa", "w") as f, open(lname + "_rightdirPROBE.fa", "w") as g, open(lname + "_rightdirPROBE_PRO.fa",  "w") as h:
with open(lname + "_rightdir.fa", "w") as f, open(lname + "_rightdirPROBE.fa", "w") as g:

	for seq in fasta_sequences:
		if seq.id not in rctaxa:
			SeqIO.write(seq, f, "fasta")
		if seq.id in rctaxa:
			SeqIO.write(seq.reverse_complement(id=seq.id, description=""), f, "fasta")
		for s in coordL:
			name,start,stop=s.split("|")
			if name == seq.id and seq.id in rctaxa:
				seqC=seq[int(start):int(stop)]
				
				seqRC=seqC.reverse_complement()
				nucrecord = SeqRecord(seqRC.seq,id=seqC.id,description="")
#				seqP=seqRC.seq.translate()
#				aarecord = SeqRecord(seqP,id=seqC.id,description="")
				SeqIO.write(nucrecord, g, "fasta")
#				SeqIO.write(aarecord, h, "fasta")
			
			
			
			if name == seq.id and seq.id not in rctaxa:
				seqC=seq[int(start):int(stop)]
#				seqP=seqC.seq.translate()
				
#				aarecord = SeqRecord(seqP,id=seqC.id,description="")
				SeqIO.write(seqC, g, "fasta")
#				SeqIO.write(aarecord, h, "fasta")
				


print str(count) + " seqs were reversed complement"
table.close()
fasta_sequences.close()
sys.exit()