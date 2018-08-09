
from py.utils.BioUtils import *; # Imports utils
from py.utils.GenBankToolbox import *; # Imports utils

#### Constants ####
# Annotation colors
annColors = {
    "LHRColor" : "#56B4E9",
    "RHRColor" : "#56B4E9",
    "gRNAColor" : "#E69F00",
    "recodedRegionColor" : "#D55E00",
    "primerColor" : "#F0E442",
    "targetGeneColor" : "#009E73",
    "gBlockColor" : "#CC79A7",
    "otherAnnColor" : "#BBBBBB"
}

# Restriction enzyme cut sequences used
cut_FseI = "ggccggcc"; # FseI cut site
cut_AsiSI = "gcgatcgc"; # AsiSI cut site
cut_IPpoI = "ctctcttaaggtagc"; # I-PpoI cut site
cut_ISceI = "attaccctgttatcccta"; # I-SceI cut site
cut_AflII = "cttaag"; # AflII cut site
cut_AhdI = "gacaacttgtc"; # AhdI cut site
cut_BsiWI = "cgtacg"; # BsiWI cut site
cut_NheI = "gctagc"; # NheI cut site

# Other sequences
ha_tag = "tacccatacgatgttccagattacgct"; # HA tag sequence

# Plasmids
pSN054_V5_Cas9 = GenBank();
pSN054_V5_Cas9.load("input/plasmids/psn054-updated_v5-tagged-jn_final.gb",loadFromFile=True); # load Cas9 plasmid sequence from GenBank format
pSN054_V5_Cas9.setAllColors(annColors['otherAnnColor'])
pSN054_V5_Cpf1 = GenBank();
pSN054_V5_Cpf1.load("input/plasmids/psn054_v5ha-tags_lbcpf1.gb",loadFromFile=True); # load Cpf1 plasmid sequence from GenBank format
pSN054_V5_Cpf1.setAllColors(annColors['otherAnnColor'])
pSN150_Cas9 = GenBank();
pSN150_Cas9.load("input/plasmids/psn150_cas9.gb",loadFromFile=True); # load Cas9 plasmid sequence from GenBank format
pSN150_Cas9.setAllColors(annColors['otherAnnColor'])
pSN150_Cpf1 = GenBank();
pSN150_Cpf1.load("input/plasmids/psn150_cpf1.gb",loadFromFile=True); # load Cas9 plasmid sequence from GenBank format
pSN150_Cpf1.setAllColors(annColors['otherAnnColor'])

# Codon usage tables
codonUsageTables = {
    'P. falciparum 3D7': codonUsage('input/codonUsageTables/Pfalciparum3D7.txt'),
    'P. vivax': codonUsage('input/codonUsageTables/Pvivax.txt'),
    'T. gondii': codonUsage('input/codonUsageTables/Tgondii.txt'),
    'E. coli K12': codonUsage('input/codonUsageTables/EcoliK12.txt'),
    'S. cerevisiae': codonUsage('input/codonUsageTables/Scerevisiae.txt'),
    'H. sapiens': codonUsage('input/codonUsageTables/Hsapiens.txt'),
    'R. norvegicus': codonUsage('input/codonUsageTables/Rnorvegicus.txt'),
    'scramble': codonUsage()
}
