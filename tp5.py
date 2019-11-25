"""Travail Pratique 5 - Alexandre Binette."""

import sys
from Bio.Blast import NCBIXML

# Importe le fichier au argv[1] et attrape les erreurs de filepath
try:
    file = sys.argv[1]
except (IOError, IndexError):
    error_message = """Error : Genbank file missing at argv[1],
    please specify filepath"""
    print(error_message, file=sys.stderr)
    exit

with open(file, "r") as file:
    blast_record = NCBIXML.read(file)

possible_hsp = []
for alignement in blast_record.alignments:
    for hsp in alignement.hsps:
        if hsp.identities / hsp.align_length == 1:
            possible_hsp.append(hsp)

longest_align = 0
for hsp in possible_hsp:
    if hsp.align_length > longest_align:
        best_hsp = hsp
