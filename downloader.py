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
    for teststring in [f"({substring})__", f"__{substring}__", f" {substring} ", substring]:
        for filename in os.listdir("seqs"):
            if teststring in filename.lower():
                return SeqIO.read(f"seqs/{filename}", "gb")

def export_genes(seqs, features, labels, outfolder):
    for featurename in features:
        genelist = []
        for i, seq in enumerate(seqs):
            #TODO: be graceful if feature does not exist (or isn't labelled)
            # either way is dangerous to assume, so best to leave the whole seq out?
            try:
                feature = next(feature for feature in seq.features 
                        if feature.qualifiers.get("gene") is not None 
                        and feature.qualifiers.get("gene")[0] == featurename)
            except Exception as e:
                print(featurename)
                print(seq.description)
                print([f.qualifiers.get("gene") for f in seq.features])
                raise e
            gene = feature.extract(seq)
            gene.description = ""
            gene.id = labels[i].strip("_() ").replace(" ", "_")
            genelist += gene,
        os.makedirs(outfolder, exist_ok=True)
        outfilename = f"{outfolder}/{featurename}.txt"
        SeqIO.write(genelist, outfilename, "fasta")

def download_entrez_idlist(cladename, max_results=40):
    print(f"Fetching ids for {cladename}")
    stream = Entrez.esearch(db="nucleotide", term=f"{cladename}[Organism] AND mitochondrion[Title] AND refseq[Filter]", retmax=max_results) 
    records = Entrez.read(stream)
    if not int(records["Count"]):
        stream = Entrez.esearch(db="nucleotide", term=f"{cladename} AND mitochondrion[Title] AND refseq[Filter]", retmax=max_results) 
        records = Entrez.read(stream)
    if not int(records["Count"]):
        raise Exception(f"Zero results in database for '{cladename}'. Typo?")
    return records["IdList"]

def download_from_id(id, featurenames):
    print(f"Fetching gene seq {id}")
    seq = SeqIO.read(Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text"), "gb")
    found_features = [f.qualifiers.get("gene")[0] for f in seq.features if f.qualifiers.get("gene") is not None]
    #check if we have all the features. if not, don't save!
    if all([f in found_features for f in featurenames]):
        print(f"Found '{seq.description}'")
        save_seq(seq, id)
        return seq



def build_from_listfile(listpath, outfilepath, featurenames):
    seqs = []
    labels = []
    with open(listpath, "r") as f:
        for line in f:
            if len(line) > 1:
                line = line.strip("\n")
                if line[0] == "*":
                    newseqs, newlabels = build_from_clade(line.lstrip("*"), featurenames)
                    seqs += newseqs
                    labels += newlabels
                else:
                    seq = open_seq_from_cache(line)
                    if seq is None:
                        idlist = download_entrez_idlist(line)
                        for id in idlist:
                            seq = download_from_id(id, featurenames)
                            if seq is not None:
                                break  # as soon as we get one valid seq, we're done with this 
                        else:
                            print(f"Found no results for '{line}'")
                            continue
                    seqs += seq,
                    labels += line,
    export_genes(seqs, featurenames, labels, outfilepath)


def build_from_clade(cladename, featurenames):
    records = download_entrez_idlist(cladename)
    seqs = []
    labels = []
    for id in records:
        seq = open_seq_from_cache(id)
        if seq is None:
            seq = download_from_id(id, featurenames)
        seqs += seq,
        string = seq.annotations["source"]
        if string.count("(") == 1 and string[string.index("("):].count(")") == 1:
            string = string[string.index("("):]
            string = string[:string.index(")")].lower()
        else:
            string = string.strip(" ")
            string = string[string.find(" ")+1:]
            string = string[:string.find("__")]
        labels += string,
    return [seqs, labels]




if __name__ == "__main__":
    import sys
    build_from_listfile(sys.argv[1], sys.argv[2], sys.argv[3:])

