# Parse localization prediction files
from sys import argv
from os import listdir
from os.path import join
from Bio import SeqIO

def parse_signalP(file_path):
    signalP = {}
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                row = line.split(" ")
                row = [e for e in row if e != '']
                if row[9] == "Y":
                    signalP.setdefault("Y",[]).append(row[0])
                else:
                    signalP.setdefault("N",[]).append(row[0])
    return signalP

def parse_predictNLS(file_path):
    predictNLS = {}
    with open(file_path, "r") as f:
        next(f)
        for line in f:
            row = line.split("\t")
            if row[4] == "Experimental":
                predictNLS.setdefault("Experimental",[]).append(row[0])
            else:
                predictNLS.setdefault("Potential",[]).append(row[0])
        return predictNLS

def parse_ASAFind(file_path):
    ASAFind = {}
    with open(file_path, "r") as f:
        next(f)
        for line in f:
            row = line.split("\t")
            if row[12] == "Plastid, high confidence":
                ASAFind.setdefault("Plastid, high confidence",[]).append(row[0])
            elif row[12] == "Plastid, low confidence":
                ASAFind.setdefault("Plastid, low confidence",[]).append(row[0])
            elif row[12] == "Not plastid, SignalP positive":
                ASAFind.setdefault("Not plastid, SignalP positive",[]).append(row[0])
            else:
                ASAFind.setdefault("Not plastid, SignalP negative",[]).append(row[0])
    return ASAFind

def parse_HECTAR(file_path):
    HECTAR = {}
    with open(file_path, "r") as f:
        next(f)
        for line in f:
            row = line.split("\t")
            if row[1] == "chloroplast":
                HECTAR.setdefault("chloroplast",[]).append(row[0].split(" ")[0])
            elif row[1] == "mitochondrion":
                HECTAR.setdefault("mitochondrion",[]).append(row[0].split(" ")[0])
            elif row[1] == "signal peptide":
                HECTAR.setdefault("signal peptide",[]).append(row[0].split(" ")[0])
            elif row[1] == "signal anchor":
                HECTAR.setdefault("signal anchor",[]).append(row[0].split(" ")[0])
            elif row[1] == "other localisation":
                HECTAR.setdefault("other localisation",[]).append(row[0].split(" ")[0])
            else:
                HECTAR.setdefault("no signal peptide or anchor",[]).append(row[0].split(" ")[0])
    return HECTAR

def parse_invalidHECTAR(file_path):
    invalid = []
    with open(file_path, "r") as f:
        next(f)
        for line in f:
            row = line.split("\t")
            invalid.append(row[0])
    return invalid

def parse_erSignal(file_path):
    erSignal = {}
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                row = line.strip(">").split(" ")
                gene = row[0]
            else:
                row = line.strip("\n").split("  ")
                erSignal.setdefault(row[len(row)-1],[]).append(gene)
    return erSignal

def parse_peroxSignal(file_path):
    peroxSignal = {}
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                row = line.strip(">").split(" ")
                gene = row[0]
            else:
                row = line.strip("\n").split("  ")
                peroxSignal.setdefault(row[len(row)-1],[]).append(gene)
    return peroxSignal

def parse_mitoprotII(file_path):
    mitoprotII = {}
    with open(file_path, "r") as f:
        next(f)
        for line in f:
            row = line.split(" ")
            row = [e for e in row if e != '']
            if float(row[3]) > 0.9:
                mitoprotII.setdefault("Score > 0.9",[]).append(row[0])
            elif float(row[3]) > 0.8:
                mitoprotII.setdefault("Score > 0.8",[]).append(row[0])
            else:
                continue
    return mitoprotII

def parse_targetP(file_path):
    targetP = {}
    with open(file_path, "r") as f:
        for line in f:
            row = line.split(" ")
            row = [e for e in row if e != '']
            if len(row) == 8:
                if row[0].startswith("N"):
                    continue
                else:
                    if int(row[7]) == 1:
                        if row[6] == "C":
                            targetP.setdefault("C1",[]).append(row[0])
                        elif row[6] == "M":
                            targetP.setdefault("M1",[]).append(row[0])
                        elif row[6] == "S":
                            targetP.setdefault("S1",[]).append(row[0])
                        else:
                            targetP.setdefault("O1",[]).append(row[0])
                    elif int(row[7]) == 2:
                        if row[6] == "C":
                            targetP.setdefault("C2",[]).append(row[0])
                        elif row[6] == "M":
                            targetP.setdefault("M2",[]).append(row[0])
                        elif row[6] == "S":
                            targetP.setdefault("S2",[]).append(row[0])
                        else:
                            targetP.setdefault("O2",[]).append(row[0])
                    elif int(row[7]) == 3:
                        if row[6] == "C":
                            targetP.setdefault("C3",[]).append(row[0])
                        elif row[6] == "M":
                            targetP.setdefault("M3",[]).append(row[0])
                        elif row[6] == "S":
                            targetP.setdefault("S3",[]).append(row[0])
                        else:
                            targetP.setdefault("O3",[]).append(row[0])
                    elif int(row[7]) == 4:
                        if row[6] == "C":
                            targetP.setdefault("C4",[]).append(row[0])
                        elif row[6] == "M":
                            targetP.setdefault("M4",[]).append(row[0])
                        elif row[6] == "S":
                            targetP.setdefault("S4",[]).append(row[0])
                        else:
                            targetP.setdefault("O4",[]).append(row[0])
                    else:
                        if row[6] == "C":
                            targetP.setdefault("C5",[]).append(row[0])
                        elif row[6] == "M":
                            targetP.setdefault("M5",[]).append(row[0])
                        elif row[6] == "S":
                            targetP.setdefault("S5",[]).append(row[0])
                        else:
                            targetP.setdefault("O5",[]).append(row[0])
            else:
                continue
    return targetP



mypath1 = argv[2]

signalP1 = parse_signalP(mypath1 + "/signalp4.1_results.txt")
predictNLS1 = parse_predictNLS(mypath1 + "/predictnls1.3_results.txt")
ASAFind1 = parse_ASAFind(mypath1 + "/ASAFind1.1_results.txt")
HECTAR1 = parse_HECTAR(mypath1 + "/HECTAR1.3_results.tab")
invalidHECTAR1 = parse_invalidHECTAR(mypath1 + "/HECTAR1.3_invalid_sequences.tab")
erSignal1 = parse_erSignal(mypath1 + "/scanProsite_results_erSignal.txt")
peroxSignal1 = parse_peroxSignal(mypath1 + "/scanProsite_results_peroxSignal.txt")
mitoprotII1 = parse_mitoprotII(mypath1 + "/mitoprotii1.101_results.txt")
targetP1 = parse_targetP(mypath1 + "/targetp1.1b_results.txt")

## Plastid-encoded and mitochondrial-encoded genes           
# Find gbk files in dir

gbk_files = [f for f in listdir('.') if f.endswith(".gbk")]
print gbk_files

# Parse gbk files
plastid_genes = []
mito_genes = []

for f in gbk_files:

    gbk_filename = f
    input_handle  = open(gbk_filename, "r")

    for seq_record in SeqIO.parse(input_handle, "genbank"):
        for seq_feature in seq_record.features:
            if seq_feature.type=="CDS":
                if f == "ThapsPlastid.gbk" or f == "NC_008589.gbk":
                    plastid_genes.append(seq_feature.qualifiers["gene"][0])
                if f == "ThapsMito.gbk" or f == "NC_007405.gbk":
                    mito_genes.append(seq_feature.qualifiers["gene"][0])

    input_handle.close()
        
print str(len(plastid_genes)), 'plastid-encoded genes'
print str(len(mito_genes)), 'mitochondria-encoded genes'

## Nuclear-encoded endoplasmic reticulum targeted genes
ER_target = []

for g in ASAFind1["Not plastid, SignalP positive"] + ASAFind1["Plastid, low confidence"]:
    for k in erSignal1.keys():
        if g in erSignal1[k] and g.startswith("T"):
            ER_target.append(g)
            
print str(len(ER_target)), 'nuclear-encoded endoplasmic reticulum targeted genes'

## Nuclear-encoded chloroplast targeted genes
plastid_target = []

for g in ASAFind1["Plastid, high confidence"] + ASAFind1["Plastid, low confidence"]:
    if g not in ER_target and g.startswith("T"):
        plastid_target.append(g)
        
print str(len(plastid_target)), 'nuclear-encoded chloroplast targeted genes'

## Nuclear-encoded mitochondria targeted genes
mito_target = []

for g in signalP1["N"]:
    if g not in predictNLS1["Experimental"]+predictNLS1["Potential"] and g.startswith("T"):
        if g in HECTAR1["mitochondrion"] and g in targetP1["M1"]+targetP1["M2"]+targetP1["M3"]+targetP1["M4"]+targetP1["M5"]:
            mito_target.append(g)
        elif g in mitoprotII1["Score > 0.9"]:
            mito_target.append(g)
        elif g in mitoprotII1["Score > 0.8"]:
            if g in HECTAR1["mitochondrion"] or g in targetP1["M1"]+targetP1["M2"]+targetP1["M3"]+targetP1["M4"]+targetP1["M5"]:
                mito_target.append(g)


print str(len(mito_target)), 'nuclear-encoded mitochondria targeted genes'


## Nuclear-encoded peroxisome targeted genes
perox_target = []

for g in signalP1["N"]:
    if g not in predictNLS1["Experimental"]+predictNLS1["Potential"]:
        if g not in mito_target and g.startswith("T"):
            for k in peroxSignal1.keys():
                if g in peroxSignal1[k]:
                    perox_target.append(g)
                    
print str(len(perox_target)), 'nuclear-encoded peroxisome targeted genes'


## Unknown localization
unknown = []
del_genes = []

files = ["deleted_genes.faa", "long_sequences.faa", "really_long_sequences.faa"]

for f in files:

    input_handle  = open(argv[1] + '/' + f, "r")

    for seq_record in SeqIO.parse(input_handle, "fasta"):
        del_genes.append(seq_record.id)
    input_handle.close()
        

for g in del_genes:
    if g not in plastid_genes + mito_genes + plastid_target + mito_target + ER_target + perox_target:
        unknown.append(g)
        
for g in invalidHECTAR1:
    if g not in plastid_genes + mito_genes + plastid_target + mito_target + ER_target + perox_target:
        unknown.append(g)
        
print str(len(unknown)), 'genes with unknown localization'

results = ["plastid_genes", "mito_genes", "er_target", "plastid_target", "mito_target", "perox_target", "locUnknown"]
lists = [plastid_genes, mito_genes, ER_target, plastid_target, mito_target, perox_target, unknown]
mypath1

for i in range(0,len(results)):
	with open(mypath1 + "/" + results[i] + ".txt", 'w') as f:
		for x in lists[i]:
			f.write(x+'\n')
