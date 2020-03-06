#!/usr/bin/env python3

"""Searches UniProt database for a specific gene, generates a MUSCLE MSA
alignment from 20 entries and produces a consensus sequence. The pickle module will
be used to serialise the alingment object. Employs PSIPRED to make a secondary
structure prediction for a particular protein sequence. The NucleicAcid class is
created that represents nucleic acids as a string seqeunce, and includes a method to
translate the sequence.

Author: Nicole Chai"""

from Bio import SeqIO
import gzip
import pickle
import re
from Bio.Align.Applications import MuscleCommandline
import os
import shutil
from Bio.Align import AlignInfo
from io import StringIO
from Bio import AlignIO
import requests
import json

MAX_OUTPUT = 20
mycgene = []
count = 0

print('20 UniProt Entries for your gene:' + '\n')
# Unzip UniProt Database file
with gzip.open("uniprot_sprot.xml.gz") as f:
    for entry in SeqIO.parse(f, "uniprot-xml"):
        gene_found = False
# Load 20 entries with primary/secondary gene names that match 'myc'. May take a while.
        gene = re.compile(r'^myc$')
        if ("gene_name_primary" in entry.annotations.keys()):
            char = entry.annotations["gene_name_primary"].lower()
            if (gene.match(char)):
                gene_found = char
            print("UniProt Id:{} Primary Gene Name:{}".format(entry.id, char)) if gene_found else ""
        if (not gene_found and "gene_name_synonym" in entry.annotations.keys()):
            chars = [y.lower() for y in entry.annotations["gene_name_synonym"]]
            for char in chars:
                if (gene.match(char)):
                    gene_found = True
                    break
            print("UniProt Id:{}, Secondary Gene Name:{}".format(entry.id, str(char))) if gene_found else ""
        if (gene_found):
            mycgene.append(entry)
            count += 1
        if (count == MAX_OUTPUT):
            break

# Write 20 entries to fasta file
with open("mycgene.fasta", 'w') as f:
    SeqIO.write(mycgene, f, "fasta")

# Generate an MSA from the 20 entries using MUSCLE
fasta_myc = MuscleCommandline(input="mycgene.fasta", out = "mycgene.html", html=True)
fasta_myc()
# Copy the html output into the web server's templates folder
shutil.copy("mycgene.html", "./site_align/align/templates/align/")

print('\n' + 'Your MUSCLE alignment:' + '\n')

# Print alignment in the terminal
cline = MuscleCommandline(input="mycgene.fasta")
stdout, stderr = cline()
alignment = AlignIO.read(StringIO(stdout), "fasta")
print(alignment)

# Pickle and unpickle the alignment object

filename = 'Your Alignment'

outfile = open(filename, 'wb')
pickle.dump(alignment, outfile)
outfile.close()

infile = open(filename, 'rb')
new_alignment = pickle.load(infile)
infile.close()

print('\n' + 'Your consensus sequence:' + '\n')

# Consensus sequence for your alignment
summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus()
print(consensus)

# Create a NucleicAcid class

class NucleicAcid:
    """Nucleic acids as a string sequence."""

    base_complement = {'C': 'G', 'A': 'T', 'G': 'C', 'T': 'A'}

    translations = {"GGT":"G","GCG":"A","CAG":"Q","TAG":"*","AGG":"R","TTT":"F","AGC":"S",
    "AAT":"N","AAG":"K","GTT":"V","TCT":"S","CGT":"R","CTT":"L","AAA":"K","GAT":"D","TTA":"L",
    "GGG":"G","CCT":"P","TGC":"C","CCG":"P","GTA":"V","AAC":"N","GGC":"G","GCT":"A",
    "TCG":"S","ACA":"T","TGA":"*","GCA":"A","CAC":"H","GAC":"D","ATA":"I","CCC":"P",
    "CGC":"R","AGA":"R","TAA":"*","AGT":"S","TCA":"S","ACC":"T","TCC":"S","ATG":"M",
    "GCC":"A","ATT":"I","ACG":"T","GTG":"V","GTC":"V","CAA":"Q","CTC":"L","CTG":"L",
    "TGT":"C","CTA":"L","GGA":"G","CAT":"H","CGA":"R","GAA":"E","GAG":"E","ATC":"I",
    "TAC":"Y","TTC":"F","TTG":"L","ACT":"T","TAT":"Y","TGG":"W","CCA":"P","CGG":"R",}

    def __init__(self, seq):
        self.seq = seq

    def complement(self):
        chars = list(self.seq)
        chars = [self.base_complement[b] for b in chars]
        return ''.join(chars)

    def codon(self):
        seq = self.seq
        length = len(seq)
        tail = length - 1 - (length % 3)
        codon = [seq[xx:xx+3] for xx in range(0, tail, 3)]
        return codon

# Returns translated AA sequence fo a DNA sequence
    def translate(self):
        seq = self.seq
        codon = self.codon()
        aminoacid = [self.translations[aminoacid] for aminoacid in codon]
        return ''.join(aminoacid)

# Open a DNA sequence file and translate the sequence, writing output into a txt file
with open('myc_dna.txt', 'r') as file:
    data = file.read().replace('\n', '')
    seq = NucleicAcid(data)
    with open('mycgene.txt', 'w') as file2:
        file2.write(seq.translate())

print('\n' + 'The UUID for your PSIPRED secondary structure prediction:' + '\n')

# Send protein sequence file to PSIPRED webserver to
# generate secondary structure prediction
url = 'http://bioinf.cs.ucl.ac.uk/psipred/api/submission.json'
files = {'input_data': ('mycgene.txt', open('./mycgene2.txt', 'rb'))}
data = {'job': 'psipred',
        'submission_name': 'myc_gene',
        'email': 'y.chai3@ncl.ac.uk', }
result = requests.post(url, data= data, files= files)
print(result.text)
