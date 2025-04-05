from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import *
from primers import create
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import MeltingTemp as mt

hydrophobicity_scale = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

def get_all_pairs(list):
    return [(list[i],list[j]) for i in range(len(list)) for j in range(i+1, len(list))]

def get_coord_before_nick(enzyme, seq):
    return enzyme.search(seq)[0] - 2

def get_site_start_coord(enzyme, seq):
    upper_nick = enzyme.search(seq)[0] - 2
    return upper_nick - len(enzyme.elucidate().replace("_", "").split("^")[0]) + 1

def get_site_end_coord(enzyme, seq):
    upper_nick = enzyme.search(seq)[0] - 2
    return upper_nick + len(enzyme.elucidate().replace("_", "").split("^")[1])

def get_aa_properties(aa):
    pa = ProteinAnalysis(aa)
    return {
        'isoelectric_point': pa.isoelectric_point(),
        'molecular_weight': pa.molecular_weight(),
        'hydrophobicity': hydrophobicity_scale.get(aa, 0)  # Default to 0 if not found
    }

def generate_codons(partial_codon):
    nucleotides = ['A', 'T', 'C', 'G']
    if len(partial_codon) == 1:
        return [partial_codon + n1 + n2 for n1 in nucleotides for n2 in nucleotides]
    elif len(partial_codon) == 2:
        return [partial_codon + n for n in nucleotides]
    else:
        raise ValueError("Partial codon must be 1 or 2 nucleotides long")

def select_codon(partial_codon, target_aa, reverse=False):
    target_properties = get_aa_properties(target_aa)
    best_codon = None
    best_match = float('inf')
    
    for codon in generate_codons(partial_codon):
        if reverse:
            codon = codon[::-1]
        aa = str(Seq(codon).translate())
        if aa == '*':
            continue
        aa_properties = get_aa_properties(aa)
        
        difference = (
            abs(aa_properties['isoelectric_point'] - target_properties['isoelectric_point']) * 0.5 +
            abs(aa_properties['molecular_weight'] - target_properties['molecular_weight']) * 0.3 +
            abs(aa_properties['hydrophobicity'] - target_properties['hydrophobicity']) * 0.2
        )
        if difference < best_match:
            best_match = difference
            best_codon = codon
    #print(best_codon, best_match)
    return best_codon