def rnatoDNA(rna: str) -> str:
    return rna.replace("U", "T")

def DNAtorna(DNA: str) -> str:
    return DNA.replace("T", "U")

def revcomp(DNA: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[nuc] for nuc in reversed(DNA))
