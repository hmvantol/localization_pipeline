from os import system
from sys import argv
from Bio import SeqIO

def batch_iterator(iterator, batch_size, total_amino_acids):
    entry = True
    while entry:
        batch = []
        aa_total = 0
        while len(batch) < batch_size and aa_total <= total_amino_acids:
            try:
                entry = iterator.next()
                aa_total += len(entry.seq)
            except StopIteration:
                entry = None
            if entry is None:
                break
            if entry.seq.strip("*").count("*") >= 1: # Remove sequences with internal stop codons
                continue
            if len(entry.seq) > 6000: # Don't append sequences longer than 6000 residues
                continue
            elif len(entry.seq) > 4000: # Don't append sequences longer than 4000 residues
                continue
            else:
                entry.id = entry.id.split(",")[0]
                entry.name = ''
                entry.description = ''
                entry.seq = entry.seq.strip("*").upper()
                batch.append(entry)
        if batch:
            aa_total = 0
            yield batch

record_iter = SeqIO.parse(open("thaps_optGeneCatalog.faa"),"fasta")
for i, batch in enumerate(batch_iterator(record_iter, 2000, 198000)):
    filename = "group_%i.fasta" % (i + 1)
    handle = open(argv[1]+"/"+filename, "w")
    count = SeqIO.write(batch, handle, "fasta")
    handle.close()
    print("Wrote %i records to %s" % (count, filename))
    
    
filename1 = "really_long_sequences.faa"
filename2 = "long_sequences.faa"
seq6000 = []
seq4000 = []
record_iter = SeqIO.parse(open("thaps_optGeneCatalog.faa"),"fasta")
for seq_record in record_iter:
    if seq_record.seq.strip("*").count("*") >= 1: # Remove sequences with internal stop codons
        continue
    seq_record.id = seq_record.id.split(",")[0]
    seq_record.name = ''
    seq_record.description = ''
    seq_record.seq = seq_record.seq.strip("*").upper()
    if len(seq_record.seq) > 6000: # Put sequences longer than 6000 residues in a separate file
        seq6000.append(seq_record)
    elif len(seq_record.seq) > 4000: # Put sequences longer than 4000 residues in a separate file
        seq4000.append(seq_record)
handle = open(argv[1]+"/"+filename1, "w")
count = SeqIO.write(seq6000, handle, "fasta")
print("Wrote %i records to %s" % (count, filename1))
handle.close()
handle = open(argv[1]+"/"+filename2, "w")
count = SeqIO.write(seq4000, handle, "fasta")
print("Wrote %i records to %s" % (count, filename2))
handle.close()


output_path = argv[1]

system("cat " + output_path + "/group_* " + output_path + "/long_sequences.faa " + output_path + "/really_long_sequences.faa > " + output_path + "/thaps_prots_all.faa")
system("cat " + output_path + "/group_* " + output_path + "/long_sequences.faa > " + output_path + "/thaps_prots_shorter.faa")
system("cat " + output_path + "/group_* > " + output_path + "/thaps_prots_shortest.faa")


int_stop = []
int_stop_id = []
count = 0
record_iter = SeqIO.parse(open("thaps_optGeneCatalog.faa"),"fasta")
for seq_record in record_iter:
    count += 1
    if seq_record.seq.strip("*").count("*") >= 1:
        int_stop.append(seq_record)
        int_stop_id.append(seq_record.id.split("_")[1])
        

handle = open(argv[1]+"/"+"deleted_genes.faa", "w")
SeqIO.write(int_stop, handle, "fasta")
handle.close()