"""
Read a protein file or text input and extract the sequence(s).
Works with PDB files (.pdb / .ent) and FASTA files (.fasta / .fa / .faa).
"""

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def parse_pdb(lines):
    """Reads SEQRES lines. Returns {chain: sequence}"""
    sequences = {}
    for line in lines:
        if line.startswith("SEQRES"):
            parts = line.split()
            if len(parts) < 5:
                continue
            chain = parts[2]
            three_letter_codes = parts[4:]

            if chain not in sequences:
                sequences[chain] = ""

            for code in three_letter_codes:
                sequences[chain] += THREE_TO_ONE.get(code, "X")
    return sequences


def parse_fasta(lines):
    """Reads FASTA text. Returns {name: sequence}"""
    sequences = {}
    name = None
    count = 0

    for line in lines:
        line = line.strip()

        if line.startswith(">"):
            count += 1
            words = line[1:].split()
            name = words[0] if words else "seq" + str(count)
            base, n = name, 2
            while name in sequences:
                name = base + "_" + str(n)
                n += 1
            sequences[name] = ""
        elif line and name is not None:
            sequences[name] += line.upper()

    return sequences


def read_pdb(filename):
    with open(filename, "r") as file:
        return parse_pdb(file)


def read_fasta(filename):
    with open(filename, "r") as file:
        return parse_fasta(file)


PDB_EXTENSIONS = (".pdb", ".ent")
FASTA_EXTENSIONS = (".fasta", ".fa", ".faa")


def read_sequences(filename):
    """Picks the right reader based on the file extension."""
    lower = filename.lower()

    if lower.endswith(PDB_EXTENSIONS):
        return read_pdb(filename)
    elif lower.endswith(FASTA_EXTENSIONS):
        return read_fasta(filename)
    else:
        raise ValueError("Unsupported file type. Use .pdb, .ent, .fasta, .fa or .faa")


def read_sequences_from_text(text, filename=None):
    """
    Reads sequences from in-memory text (uploaded file content or pasted text).
    """
    lines = text.splitlines()

    if filename:
        lower = filename.lower()
        if lower.endswith(PDB_EXTENSIONS):
            return parse_pdb(lines)
        if lower.endswith(FASTA_EXTENSIONS):
            return parse_fasta(lines)
        raise ValueError("Unsupported file type. Use .pdb, .ent, .fasta, .fa or .faa")

    if any(line.startswith("SEQRES") for line in lines):
        return parse_pdb(lines)
    first = next((l.strip() for l in lines if l.strip()), "")
    if first.startswith(">"):
        return parse_fasta(lines)
    return {"pasted_sequence": text} if text.strip() else {}