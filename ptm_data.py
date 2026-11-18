"""
Reference data for the PTM Site Scanner (Group 15).
"""

AMINO_ACID_NAMES = {
    "A": "Alanine",       "C": "Cysteine",     "D": "Aspartic acid",
    "E": "Glutamic acid", "F": "Phenylalanine", "G": "Glycine",
    "H": "Histidine",     "I": "Isoleucine",   "K": "Lysine",
    "L": "Leucine",       "M": "Methionine",   "N": "Asparagine",
    "P": "Proline",       "Q": "Glutamine",    "R": "Arginine",
    "S": "Serine",        "T": "Threonine",    "V": "Valine",
    "W": "Tryptophan",    "Y": "Tyrosine",
    "X": "Unknown / non-standard",
}

# Residue (one-letter code) -> modifications that can occur on it
MODIFICATION_TYPES = {
    "S": ["Phosphorylation", "O-glycosylation"],
    "T": ["Phosphorylation", "O-glycosylation"],
    "Y": ["Phosphorylation"],
    "K": ["Acetylation", "Methylation", "Ubiquitination", "Sumoylation"],
    "R": ["Methylation"],
    "N": ["N-glycosylation"],
    "C": ["S-palmitoylation"],
    "P": ["Hydroxylation"],
}

# Short plain-language descriptions shown in the interface
MOD_INFO = {
    "Phosphorylation": "A phosphate group is added (by kinases) to Ser, Thr or Tyr. Common switch for signalling.",
    "Acetylation": "An acetyl group is added to Lys. Well known on histones; affects charge and protein interactions.",
    "Methylation": "Methyl groups are added to Lys or Arg. Important in histone and RNA-binding proteins.",
    "N-glycosylation": "A sugar chain is attached to Asn, normally inside an N-X-S/T sequon.",
    "O-glycosylation": "A sugar is attached to the hydroxyl group of Ser or Thr.",
    "Ubiquitination": "Ubiquitin is attached to Lys. Often marks a protein for degradation.",
    "Sumoylation": "SUMO protein is attached to Lys, often within a hydrophobic-K-x-D/E motif.",
    "S-palmitoylation": "A fatty acid is attached to Cys, helping anchor proteins to membranes.",
    "Hydroxylation": "A hydroxyl group is added to Pro (best known in collagen).",
}

# Motif rules (motif-based detection)
MOTIFS = [
    {"name": "N-glycosylation sequon (N-X-S/T, X not P)",
     "pattern": "N[^P][ST]", "offset": 0, "modification": "N-glycosylation"},
    {"name": "Proline-directed kinase site (S/T-P)",
     "pattern": "[ST]P", "offset": 0, "modification": "Phosphorylation"},
    {"name": "Basophilic kinase site, PKA-like (R-R/K-x-S/T)",
     "pattern": "R[RK].[ST]", "offset": 3, "modification": "Phosphorylation"},
    {"name": "Acidic kinase site, CK2-like (S/T-x-x-D/E)",
     "pattern": "[ST]..[DE]", "offset": 0, "modification": "Phosphorylation"},
    {"name": "SUMO consensus (hydrophobic-K-x-D/E)",
     "pattern": "[VILMF]K.[DE]", "offset": 1, "modification": "Sumoylation"},
    {"name": "Arginine methylation, RG motif",
     "pattern": "RG", "offset": 0, "modification": "Methylation"},
]