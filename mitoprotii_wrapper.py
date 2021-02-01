# python mitoprotii_wrapper.py $1 $2

from Bio import SeqIO
from os import system, path
import sys


file = open(sys.argv[2] + "mitoprotii1.101_results.txt", "w")
file.write("#\tlocus_tag\tseq_length\tDFM_NOT-MITO\tDFM_MITO(/CHLORO)\tDFMC_NOT-MITO\tDFMC_MITO(/CHLORO)\n")
file.close()

record_iter = SeqIO.parse(open(sys.argv[1] + "thaps_prots_all.faa"),"fasta")

for seq_record in record_iter:
	handle = open(sys.argv[2] + seq_record.id + ".faa", "w")
	handle.write(">%s\n%s\n" %(seq_record.id, seq_record.seq))
	handle.close()
	system("mitoprot -f a -o " + sys.argv[2] + "mitoprotii1.101_results.txt " + sys.argv[2] + seq_record.id + ".faa")
	if path.isfile(sys.argv[2] + seq_record.id + ".faa"):
		system("rm " + sys.argv[2] + seq_record.id + ".faa")
	if path.isfile(seq_record.id + ".mitoprot"):
		system("rm " + seq_record.id + ".mitoprot")

