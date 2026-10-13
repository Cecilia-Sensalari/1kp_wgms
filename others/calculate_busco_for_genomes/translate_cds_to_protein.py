#!/usr/bin/env python3
"""Translate a CDS FASTA file into a protein FASTA file (standard genetic code), e.g. as input for BUSCO.

Sequence names are kept (first word of the header). The trailing stop codon is removed; a stop codon inside
the sequence is written as X, and so is any codon with an ambiguous base. A final incomplete codon is ignored.
The input is never changed and an existing output file is never overwritten.

Usage:   translate_cds_to_protein.py INPUT_CDS.fa OUTPUT_PROTEINS.fa
"""

import sys
from itertools import product
from pathlib import Path

BASES = "TCAG"
AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {"".join(codon): aa for codon, aa in zip(product(BASES, repeat=3), AMINO_ACIDS)}


def translate(cds: str) -> str:
    cds = cds.upper().replace("U", "T")
    protein = [CODON_TABLE.get(cds[i:i + 3], "X") for i in range(0, len(cds) - len(cds) % 3, 3)]
    if protein and protein[-1] == "*":
        protein.pop()
    return "".join(protein).replace("*", "X")


def read_fasta(path: Path):
    name, chunks = None, []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                fields = line[1:].split(None, 1)
                name, chunks = (fields[0] if fields else ""), []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    source, target = Path(sys.argv[1]), Path(sys.argv[2])
    if not source.is_file():
        sys.exit(f"ERROR: cannot read [{source}]")
    if target.exists():
        sys.exit(f"ERROR: [{target}] already exists, not overwriting it")

    count = 0
    with target.open("w") as out:
        for name, cds in read_fasta(source):
            protein = translate(cds)
            if not protein:
                continue
            out.write(f">{name}\n")
            for start in range(0, len(protein), 60):
                out.write(protein[start:start + 60] + "\n")
            count += 1
    print(f"{source.name}: {count} proteins written to {target}")


if __name__ == "__main__":
    main()
