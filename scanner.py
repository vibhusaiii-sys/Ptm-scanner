"""
Core scanning engine: clean + validate -> look up residues -> motifs -> report.
"""
import re
from ptm_data import AMINO_ACID_NAMES, MODIFICATION_TYPES, MOTIFS, MOD_INFO


class SequenceError(ValueError):
    pass


_WHITESPACE_AND_DIGITS = re.compile(r"[\s\d]+")


def clean_sequence(raw, lenient=False):
    """
    Remove spaces, line breaks and position numbers, convert to capitals,
    drop a trailing stop symbol (*), and reject invalid characters.
    With lenient=True invalid characters become X instead of raising an error.
    Returns (sequence, list_of_invalid_characters_found).
    """
    seq = _WHITESPACE_AND_DIGITS.sub("", raw).upper().rstrip("*")
    invalid = sorted({c for c in seq if c not in AMINO_ACID_NAMES})

    if invalid:
        if not lenient:
            raise SequenceError("Invalid character(s): " + " ".join(invalid))
        seq = "".join(c if c in AMINO_ACID_NAMES else "X" for c in seq)
    return seq, invalid


def _as_list(value):
    return [value] if isinstance(value, str) else list(value)


def scan_sequence(seq):
    """Returns a list of sites, each with its possible modifications."""
    sites = {}

    # Step 1: dictionary lookup for every residue
    for position, residue in enumerate(seq, start=1):
        if residue in MODIFICATION_TYPES:
            sites[position] = {
                "position": position,
                "residue": residue,
                "residue_name": AMINO_ACID_NAMES[residue],
                "modifications": {
                    m: {"name": m, "evidence": "residue", "motifs": []}
                    for m in _as_list(MODIFICATION_TYPES[residue])
                },
            }

    # Step 2: motif detection (lookahead so overlapping matches are found)
    for motif in MOTIFS:
        for match in re.finditer("(?=(" + motif["pattern"] + "))", seq):
            text = match.group(1)
            start = match.start()
            site = sites.get(start + motif["offset"] + 1)
            if site is None:
                continue
            entry = site["modifications"].get(motif["modification"])
            if entry is None:
                continue
            entry["evidence"] = "motif"
            entry["motifs"].append({
                "name": motif["name"], "match": text,
                "start": start + 1, "end": start + len(text),
            })

    result = []
    for position in sorted(sites):
        site = sites[position]
        mods = list(site["modifications"].values())
        site["modifications"] = mods
        site["evidence"] = "motif" if any(m["evidence"] == "motif" for m in mods) else "residue"
        result.append(site)
    return result


def scan_sequences(sequences, lenient=False):
    """sequences: {name: raw_sequence}. Returns a report dictionary."""
    chains = []
    mod_names = []

    for name, raw in sequences.items():
        try:
            seq, invalid = clean_sequence(raw, lenient)
        except SequenceError as err:
            chains.append({"name": name, "error": str(err)})
            continue
        if not seq:
            chains.append({"name": name, "error": "Sequence is empty."})
            continue

        sites = scan_sequence(seq)
        counts = {}
        for site in sites:
            for m in site["modifications"]:
                c = counts.setdefault(m["name"], {"total": 0, "motif": 0})
                c["total"] += 1
                c["motif"] += m["evidence"] == "motif"
                if m["name"] not in mod_names:
                    mod_names.append(m["name"])

        warnings = []
        if invalid:
            warnings.append("Replaced invalid character(s) with X: " + " ".join(invalid))
        x_fraction = seq.count("X") / len(seq)
        if x_fraction > 0.5:
            warnings.append("Over half of this chain is non-standard residues (X); "
                            "it may be DNA/RNA or a ligand rather than a protein.")

        chains.append({
            "name": name, "sequence": seq, "length": len(seq),
            "sites": sites, "counts": counts, "warnings": warnings,
        })

    return {"chains": chains, "modifications": mod_names, "mod_info": MOD_INFO}