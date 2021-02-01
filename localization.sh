#!/bin/bash
# Wrapper for protein localization tools

mkdir -p $1
mkdir -p $2

# Split whole proteome fasta file into components for tool restrictions
python tool_restrictions.py $1

NFILES=$(ls -l $1/group_* | wc -l)

# SignalP 4.1
echo "Running SignalP..."
for ((i = 1; i <= ${NFILES}; i += 1)); do
	signalp $1/group_${i}.fasta >> $2/signalp4.1_results.txt
done
signalp $1/long_sequences.faa >> $2/signalp4.1_results.txt

# ASAFind 1.1
echo "Running ASAFind..."
python ASAFind.py -f $1/thaps_prots_shorter.faa -p $2/signalp4.1_results.txt -o $2/ASAFind1.1_results.txt

# MitoprotII 1.101
echo "Running MitoprotII..."
python mitoprotii_wrapper.py $1 $2

# scanProsite (ps_scan)
echo "Running scanProsite..."
ps_scan.pl -d /share/apps/prosite/prosite.dat -p "[SAC]-[KRH]-[LM]>" -p "S-S-L>" -p PS00342 $1/thaps_prots_all.faa >> $2/scanProsite_results_peroxSignal.txt
ps_scan.pl -d /share/apps/prosite/prosite.dat -p "[KD]-[DE]-E-L>" -p PS00014 $1/thaps_prots_all.faa >> $2/scanProsite_results_erSignal.txt

# TargetP 1.1
echo "Running TargetP..."
for ((i = 1; i <= ${NFILES}; i += 1)); do
	targetp -P $1/group_${i}.fasta >> $2/targetp1.1b_results.txt
done

# predictNLS 1.3
echo "Running predictNLS..."
python predictnls.py $1/thaps_prots_all.faa $2/predictnls1.3_results.txt My_NLS_list

# Parse localizations
python parse_localization.py $1 $2