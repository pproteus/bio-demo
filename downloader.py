from Bio import Entrez, SeqIO
import os

try:
    Entrez.email = open("email.txt", "r").read()
except FileNotFoundError as e:
    print("Use your own email!")
    raise e


def save_seq(seqrecord, idnum):
    filepath = f"seqs/{seqrecord.annotations["source"]}__{idnum}__.txt"
    with open(filepath, "w") as f: 
        SeqIO.write(seqrecord, f, "gb")

def open_seq_from_cache(substring):
    substring = substring.lower()
    for filename in os.listdir("seqs"):
        if substring in filename.lower():
            return SeqIO.read(f"seqs/{filename}", "gb")

def export_genes(seqs, features, labels, outfilename):
    genelist = []
    for featurename in features:
        for i, seq in enumerate(seqs):
            #TODO: be graceful if feature does not exist (or isn't labelled)
            # either way is dangerous to assume, so best to leave the whole seq out?
            feature = next(feature for feature in seq.features 
                        if feature.qualifiers.get("gene") is not None 
                        and feature.qualifiers.get("gene")[0] == featurename)
            gene = feature.extract(seq)
            gene.description = ""
            gene.id = labels[i].strip("_() ").replace(" ", "_")
            genelist += gene,
    outfilepath = f"genes/{outfilename}"
    SeqIO.write(genelist, outfilename, "fasta")
    return outfilepath

def download_entrez_idlist(cladename, max_results):
    print(f"Fetching ids for {cladename}")
    stream = Entrez.esearch(db="nucleotide", term=f"{cladename}[Organism] AND mitochondrion[Title] AND refseq[Filter]", retmax=max_results) 
    records = Entrez.read(stream)
    return records["IdList"]

def download_from_id(id):
    print(f"Fetching gene seq {id}")
    seq = SeqIO.read(Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text"), "gb")
    save_seq(seq, id)
    return seq



def build_from_listfile(filepath, outfilepath):
    seqs = []
    labels = []
    with open(filepath, "r") as f:
        for line in f:
            if len(line) > 1:
                line = line.strip("\n")
                seq = open_seq_from_cache(f"{line}__")
                if seq is None:
                    try:
                        id = download_entrez_idlist(line, 1)[0]
                    except IndexError:
                        print(f"Found no results for '{line}'")
                        continue
                    seq = download_from_id(id)
                seqs += seq,
                labels += line,
    export_genes(seqs, DESIRED_FEATURENAMES, labels, outfilepath)


def build_from_clade(cladename, outfilepath, max_results=100):
    records = download_entrez_idlist(cladename, max_results)
    seqs = []
    labels = []
    for id in records:
        seq = open_seq_from_cache(f"__{id}__")
        if seq is None:
            seq = download_from_id(id)
        seqs += seq,
        string = seq.annotations["source"]
        if string.count("(") == 1 and string[string.index("("):].count(")") == 1:
            string = string[string.index("("):]
            string = string[:string.index(")")]
        else:
            string = string.strip(" ")
            string = string[string.find(" ")+1:]
            string = string[:string.find("__")]
        labels += string,
    export_genes(seqs, DESIRED_FEATURENAMES, labels, outfilepath)



if __name__ == "__main__":
    DESIRED_FEATURENAMES = ['COX1'] #['ATP6','ATP8','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND4L','ND5','ND6']
    build_from_listfile("testlist.txt", "testout.txt")

