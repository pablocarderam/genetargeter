
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

# Cut sites to check for each plasmid
cut_sites = {
    'pSN054' : [cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI],
    'pSN150' : [cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI,cut_NheI],
    'pSN150-KO' : [cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI,cut_NheI],
    'pSN054_V5' : [cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI],
    'pSN150-Ter' : [cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI,cut_NheI],
    'pSN150-KO-Ter' : [cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI,cut_NheI]
    }

# Other sequences
ha_tag = "tacccatacgatgttccagattacgct"; # HA tag sequence

# Plasmids
pSN054_Cas9 = GenBank();
pSN054_Cas9.load("input/plasmids/psn054_cas9.gb",loadFromFile=True); # load Cas9 plasmid sequence from GenBank format
pSN054_Cas9.setAllColors(annColors['otherAnnColor'])
pSN054_V5_Cas12 = GenBank();
pSN054_V5_Cas12.load("input/plasmids/pSN054_cas12.gb",loadFromFile=True); # load Cas12 plasmid sequence from GenBank format
pSN054_V5_Cas12.setAllColors(annColors['otherAnnColor'])
pSN054_V5_Cas9 = GenBank();
pSN054_V5_Cas9.load("input/plasmids/psn054-updated_v5-tagged-jn_final.gb",loadFromFile=True); # load Cas9 plasmid sequence from GenBank format
pSN054_V5_Cas9.setAllColors(annColors['otherAnnColor'])
pSN054_V5_Cas12 = GenBank();
pSN054_V5_Cas12.load("input/plasmids/psn054_v5ha-tags_lbcas12.gb",loadFromFile=True); # load Cas12 plasmid sequence from GenBank format
pSN054_V5_Cas12.setAllColors(annColors['otherAnnColor'])
pSN150_Cas9_Ter = GenBank();
pSN150_Cas9_Ter.load("input/plasmids/psn150_cas9_editedpcr.gb",loadFromFile=True); # load Cas9 plasmid sequence from GenBank format
pSN150_Cas9_Ter.setAllColors(annColors['otherAnnColor'])
pSN150_Cas12_Ter = GenBank();
pSN150_Cas12_Ter.load("input/plasmids/psn150_cas12_editedpcr.gb",loadFromFile=True); # load Cas9 plasmid sequence from GenBank format
pSN150_Cas12_Ter.setAllColors(annColors['otherAnnColor'])
pSN150_Cas9 = GenBank();
pSN150_Cas9.load("input/plasmids/psn150_cas9.gb",loadFromFile=True); # load Cas9 plasmid sequence from GenBank format
pSN150_Cas9.setAllColors(annColors['otherAnnColor'])
pSN150_Cas12 = GenBank();
pSN150_Cas12.load("input/plasmids/psn150_cas12.gb",loadFromFile=True); # load Cas9 plasmid sequence from GenBank format
pSN150_Cas12.setAllColors(annColors['otherAnnColor'])

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

# SignalP DB
fOb = open('input/signalP/Pf3D7_GenesWithSignalPeptide_Summary.csv');
signalPDB = fOb.read();
fOb.close();

# Assembly instructions
instructions_pSN054 = '\n*** Assebly instructions for pSN054-type constructs: ***\n\n   1. Obtain primers in Oligo csv file and gene fragments in gBlock fasta file \n      from DNA synthesis\n   2. PCR the RHR and LHR fragments from genomic DNA using the corresponding \n      Gibson Assembly primers in Oligo csv file\n   2.5 If design includes a recoded region and uses a Klenow fragment for its \n      assembly instead of a pre-synthesized gene fragment, run a Klenow reaction with \n      the corresponding oligos to obtain recoded region\n   3. Digest empty parent vector with FseI and AsiSI restriction enzymes\n   4. Gibson Assembly to insert LHR into digestion product; if design includes \n      recoded region, include in assembly as well\n   5. Digest resulting vector with restriction enzyme ISceI\n   6. Gibson Assembly to insert RHR into digestion product\n   6.5 If using a Klenow fragment for sgRNA instead of a pre-synthesized gene \n      fragment, run a Klenow reaction with the corresponding oligos to obtain sgRNA\n   7. Digest resulting vector with restriction enzyme IPpoI\n   8. Gibson Assembly to insert sgRNA cassette or Klenow fragment into digestion \n      product\n   9. Check all steps by sequencing\n   10. Transfect into Cas9 or Cas12-containing cell lines!'

instructions_pSN150 = '\n*** Assebly instructions for pSN150-type constructs for conditional knock-down: ***\n\n   1. Obtain primers in Oligo csv file and gene fragments in gBlock fasta file \n      from DNA synthesis\n   2. PCR the RHR and LHR fragments from genomic DNA using the corresponding \n      Gibson Assembly primers in Oligo csv file\n   3. Digest empty parent vector with restriction enzyme FseI\n   4. Gibson Assembly to insert LHR into digestion product\n   4.5 If design includes a recoded region and uses a Klenow fragment for its \n      assembly instead of a pre-synthesized gene fragment, run a Klenow reaction with \n      the corresponding oligos to obtain recoded region\n   5. Digest resulting vector with restriction enzyme AhdI (add Nhe1 if removing \n      leading HA tag for better expression)\n   6. Gibson Assembly to insert RHR into digestion product; if design includes \n      recoded region, include in assembly as well\n   6.5 If using a Klenow fragment for sgRNA instead of a pre-synthesized gene \n      fragment, run a Klenow reaction with the corresponding oligos to obtain sgRNA\n   7. Digest resulting vector with restriction enzyme IPpoI\n   8. Gibson Assembly to insert sgRNA cassette or Klenow fragment into digestion \n      product\n   9. Check all steps by sequencing\n   10. Transfect into Cas9 or Cas12-containing cell lines!'

instructions_pSN150_KO = '\n*** Assebly instructions for pSN150-type constructs for knock-out: ***\n\n   1. Obtain primers in Oligo csv file and gene fragments in gBlock fasta file \n      from DNA synthesis\n   2. PCR the RHR and LHR fragments from genomic DNA using the corresponding\n      Gibson Assembly primers in Oligo csv file\n   3. Digest empty parent vector with restriction enzyme FseI\n   4. Gibson Assembly to insert LHR into digestion product\n   4.5 If design includes a recoded region and uses a Klenow fragment for its \n      assembly instead of a pre-synthesized gene fragment, run a Klenow reaction with \n      the corresponding oligos to obtain recoded region\n   5. Digest resulting vector with restriction enzyme AsiSI\n   6. Gibson Assembly to insert RHR into digestion product; if design includes \n      recoded region, include in assembly as well\n   6.5 If using a Klenow fragment for sgRNA instead of a pre-synthesized gene \n      fragment, run a Klenow reaction with the corresponding oligos to obtain sgRNA\n   7. Digest resulting vector with restriction enzyme IPpoI\n   8. Gibson Assembly to insert sgRNA cassette or Klenow fragment into digestion \n      product\n   9. Check all steps by sequencing\n   10. Transfect into Cas9 or Cas12-containing cell lines!'

instructions = {
'pSN054':instructions_pSN054,
'pSN150':instructions_pSN150,
'pSN150-KO':instructions_pSN150_KO
}
