"""Travail Pratique 5 - Alexandre Binette."""

import sys
import json
import pandas as pd
from urllib.request import urlopen
from Bio.Blast import NCBIXML


def get_request(query_url):
    """Returns the result of a URL GET in json format with urlopen"""

    stream = urlopen(query_url)
    result = json.loads(stream.read().decode())
    return result


# Importe le fichier au argv[1] et attrape les erreurs de filepath
try:
    file = sys.argv[1]
except (IOError, IndexError):
    error_message = """Error : Genbank file missing at argv[1],
    please specify filepath"""
    print(error_message, file=sys.stderr)
    exit

# Établit le blast_record avec le fichier
with open(file, "r") as file:
    blast_record = NCBIXML.read(file)

# Trouver les alignements valides et trouver le plus long
longest_align = 0
possible_align = []

for alignement in blast_record.alignments:
    for hsp in alignement.hsps:
        if hsp.identities / hsp.align_length == 1:
            possible_align.append(alignement)

for alignement in possible_align:
    for hsp in alignement.hsps:
        if hsp.align_length > longest_align:
            longest_align = hsp.align_length
            symbol = alignement.title.split("|")[3]

# Construire l'URL du xref
url = "http://rest.ensembl.org/xrefs/symbol/{species}/{symbol}"
opt = ["content-type=application/json", "feature=gene"]
species = "homo_sapiens"

url_xrefs = "{}?{}".format(url, ";".join(opt)).format(
    species=species, symbol=symbol
    )

xref = get_request(url_xrefs)

# Construire l'URL du overlap
id = xref[0].get("id")
url = "http://rest.ensembl.org/overlap/id/{id}"
opt = ["?content-type=application/json;feature=variation"]

url_overlap = "{}{}".format(url, ";".join(opt)).format(
    id=id,
    )

overlap = get_request(url_overlap)

# Établir le DataFrame avec les variables au CSV et l'imprimer au fichier
df = pd.DataFrame.from_dict(overlap)

variables = {
    "name": "id", "chromosome": "seq_region_name",
    "start": "start", "end": "end",
    "alleles": "alleles", "consequence": "consequence_type",
    "clinical": "clinical_significance"
}

df = df[list(variables.values())]
df.columns = list(variables.keys())

df.to_csv("./{}.variations.csv".format(id))
