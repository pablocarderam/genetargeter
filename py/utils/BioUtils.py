# -*- coding: utf-8 -*-
"""
Created on Sun May 17 17:55:17 2015
Methods useful in bioinformatic algorithms
@author: Pablo the awesome molecular jedi
"""

# Libraries
from __future__ import division; # Needed in order to avoid weird rounding practices in Python 2
from builtins import str
from builtins import range
from math import *; # Used for logarithms
import random; # Needed to generate random color codes
import datetime; # Needed to access dates

from Bio.SeqUtils import MeltingTemp as mt # needed for NN melting temp

# Methods

# File handling
"""
Takes file with multiple sequences in Fasta format,
Returns a dictionary of keys=sequence Fasta names, values=sequences
"""
def loadFastas(pFile):
    txt = open(pFile); # Access given file
    d = txt.read(); # Read file
    txt.close(); # Close file

    d = d.lstrip(">"); # Remove initial >
    d = d.split(">"); # Split into separate sequences

    d[:] = [i.split("\n", 1) for i in d]; # Split each sequence (only first occurrence of \n) into lines: first line contains sequence name. Returns d as a list of lists (sequences) containing strings (seq name, seq).
    strNames = [i[0] for i in d]; # Stores sequence names in separate list.

    d[:] = [i[1].split("\n") for i in d]; # Split each sequence into its lines. Keeps only sequences. D is now a list of lists (sequences) containing strings (lines of sequence).
    d[:] = ["".join(i) for i in d]; # Join each sequence into one string (per sequence).

    seqDic = dict(zip(strNames, d)); # Creates dictionary with seq names as keys and seq strings as values.
    return seqDic; # returns dictionary


"""
Takes file with multiple DNA sequences separated by a newline \n escape sequence,
Returns a list of strings (sequences)
"""
def loadSeqs(pFile):
    txt = open(pFile); # Access given file
    d = txt.read(); # Read file

    d = d.split("\n"); # Split each sequence into its lines. Keeps only sequences. D is now a list containing strings (sequences).
    txt.close(); # close file
    return d; # returns dictionary


"""
Output data to fileName.
"""
def output(out, fileName, wipe=False, maxSizeFragment=10000):
    mode = "a"; # default mode is to append
    if wipe:
        mode = "w+"; # change to overwrite file if wiping

    chunks = [out[i:i+maxSizeFragment] for i in range(0, len(out), maxSizeFragment)]; # divide filestring into chunks
    fOut = open(fileName,mode); # create or open file in appending mode
    c = 0; # counter if keeping track of progress
    for chunk in chunks: # break into smaller chunks and write
        # print c # print progress
        fOut.write(chunk); # write output to file
        c +=1 # advance counter

    fOut.close(); # close file




# Sequence manipulation
"""
Takes two strings,
returns True if they contain matching nucleotides.
Can handle multi-case and "N" nucleotides.
"""
def compareSeq(pSeq1, pSeq2):
    seq1 = pSeq1.upper(); # makes seq all uppercase to enable multi-case comparing
    seq2 = pSeq2.upper(); # makes seq all uppercase to enable multi-case comparing
    r = seq1 == seq2; # compares sequences
    if not r and len(seq1) == len(seq2): # if not equal and same size, pay closer attention
        if seq1.find("N") < 0 or seq2.find("N") < 0: # if seqs contain N character
            r = True; # reset response variable to true
            for i in range(0,len(seq1)):
                if seq1[i] != seq2[i] and seq1[i] != "N" and seq2[i] != "N":
                    r = False;
                    break; # escapes iteration
            pass
    return r; # returns result


"""
Takes two strings,
returns index of first occurrence of second string inside first
"""
def findFirst(pSeq, pMotif):
    seq = pSeq.upper(); # makes seq all uppercase to enable multi-case searching
    motif = pMotif.upper(); # makes motif all uppercase to enable multi-case searching
    index = seq.find(motif); # find index of motif occurrence

    return index;


"""
Takes two strings,
returns list with indexes of occurrences of second string inside first
"""
def findMotif(pSeq, pMotif):
    indexes = []; # Stores indexes of all occurences of the given string
    start = 0; # Starting index for the search
    notDone = True;
    while notDone: # As long as end not reached
        seq = pSeq.upper(); # makes seq all uppercase to enable multi-case searching
        motif = pMotif.upper(); # makes motif all uppercase to enable multi-case searching
        index = seq.find(motif, start); # find index of motif occurrence
        if index == -1: # if not found
            notDone = False; # search done
        else: # else (if found)
            start = index+1; # Move index start number
            indexes.append(index); # store index

    #if len(indexes) == 0: # If motif has not been found,
    #    indexes = "Motif not found."; # say so.

    return indexes;


"""
Counts As, Ts, Cs, and Gs in DNA sequence. Returns dictionary with nucleotides as keys.
"""
def count(seq):
    atcg = {"A":(seq.count("A")+seq.count("a")), "T":(seq.count("T")+seq.count("t")), "G":(seq.count("G")+seq.count("g")), "C":(seq.count("C")+seq.count("c")) };
    return atcg;


"""
Returns GC content of sequence as a fraction.
"""
def gcContent(pSeq):
    seq = pSeq.upper(); # uppercases everything
    GC = 0; # counter
    for bp in seq: # for every letter in sequence
        if bp == "G" or bp == "C": # if letter is G or C
            GC = GC + 1; # add one to GC counter

    GC = GC/len(seq); # Stores GC content for this sequence as decimal. Since using future division, we're fine with the division operator on ints.
    return GC;


"""
Takes dictionary of keys=sequence Fasta names, values=sequences
Returns list with max GC content, name of sequence with max GC content
"""
def maxGCContent(pData):
    max = [0, 0]; # stores maximum GC
    for seq in list(pData.keys()): # for every sequence
        GC = 0; # Counts number of occurrences of G or C
        for j in pData[seq]: # for every letter in sequence
            if j == "G" or j == "C" or j == "g" or j == "c" : # if letter is G or C
                GC += 1; # add one to GC counter

            if GC > max[0]: # if GC content of this sequence is greater than the previous maximum,
                max = [GC, seq]; # store this content and the sequence name as maximum

    max[0] = 100*max[0]/len(pData[max[1]]); # Stores GC content for this sequence as percentage
    return max;


"""
Returns approximate melting temperature of sequence to its complement using
Howley formula. Cation concentration is in moles/L (M). RNA-DNA interactions are
equivalent to RNA-RNA interactions.
http://www.sigmaaldrich.com/technical-documents/articles/biology/oligos-melting-temp.html
"""
def meltingTempSimple(seq, cationConcentration=0.05, interactionType="DNA"):
    Tm = 81.5 + 41*(gcContent(seq)) - 500/len(seq) + 16.6*log(cationConcentration,10);
    if interactionType == "RNA":
        Tm = 79.8 + 58.4*(gcContent(seq)) + 11.8*(gcContent(seq))**2 - 820/len(seq) + 18.5*log(cationConcentration,10);

    return Tm;


"""
Returns approximate melting temperature of sequence to its complement using
nearest neighbor algorithm (Biopython wrapper). Cation concentration is in
moles/L (M).
"""
def meltingTemp(seq, cationConcentration=0.05, interactionType="DNA"):
    if interactionType=="DNA":
        Tm = mt.Tm_NN(seq,Na=cationConcentration*1000, nn_table=mt.DNA_NN3); # Allawi & SantaLucia (1997) (default)
    if interactionType == "RNA":
        Tm = mt.Tm_NN(seq,Na=cationConcentration*1000, nn_table=mt.RNA_NN3); # Chen et al. (2012)

    return Tm;


"""
Takes dictionary of keys=sequence Fasta names, values=sequences
Return dictionary with { consensus: consensus sequence (string), nucleotide: profile (list of nucleotide frequency in each position)
"""
def consensus(data):
    l = len(list(data.values())[0]); # Stores length of DNA sequences analyzed
    C = { "consensus":"", "A":[0]*l, "T":[0]*l, "C":[0]*l, "G":[0]*l } # Dictionary with consensus sequence and lists of size l, stores consensus sequence and occurrences of nucleotide in all sequences for each position

    # Build profile matrix:
    for seq in list(data.values()): # For every sequence
        for b in range(0,l): # For every base in sequence
            C[seq[b]][b] = C[seq[b]][b] + 1; # Add one to count of corresponding nucleotide at corresponding position

    # Build consensus sequence:
    for p in range(0,l): # For every position in sequence
        S = { C["A"][p]:"A", C["T"][p]:"T", C["C"][p]:"C", C["G"][p]:"G" } # Define a dictionary that relates the number of occurrences in all sequences of each nucleotide at the current position to the nucleotide's letter
        C["consensus"] = C["consensus"] + S[max(S.keys())]; # Adds nucleotide letter with most occurrences at this position to the consensus sequence

    return C;


"""
Transcribes DNA to RNA
"""
def transcribe(dna):
    rna = dna.replace("T", "U");
    rna = rna.replace("t", "u");
    return rna;


"""
Reverse transcribes RNA to DNA
"""
def revTranscribe(rna):
    dna = rna.replace("U", "T");
    dna = dna.replace("u", "t");
    return dna;


"""
Returns genetic code dictionary. Stop codons are "*"
"""
def geneticCode():
    code = {'CGA':'R','UUU':'F','UUC':'F',
      'UUA':'L','UUG':'L','UCU':'S',
      'UCC':'S','UCA':'S','UCG':'S',
      'UAU':'Y','UAC':'Y','UAA':'*',
      'UAG':'*','UGU':'C','UGC':'C',
      'UGA':'*','UGG':'W','CUU':'L',
      'CUC':'L','CUA':'L','CUG':'L',
      'CCU':'P','CCC':'P','CCA':'P',
      'CCG':'P','CAU':'H','CAC':'H',
      'CAA':'Q','CAG':'Q','CGU':'R',
      'CGC':'R','AGA':'R','CGG':'R',
      'AUU':'I','AUC':'I','AUA':'I',
      'AUG':'M','ACU':'T','ACC':'T',
      'ACA':'T','ACG':'T','AAU':'N',
      'AAC':'N','AAA':'K','AAG':'K',
      'AGU':'S','AGC':'S','GGA':'G',
      'AGG':'R','GUU':'V','GUC':'V',
      'GUA':'V','GUG':'V','GCU':'A',
      'GCC':'A','GCA':'A','GCG':'A',
      'GAU':'D','GAC':'D','GAA':'E',
      'GAG':'E','GGU':'G','GGC':'G',
      'GGG':'G'
      };
    return code;


"""
Returns dictionary with keys=codons and values=frequencies for the codon table
whose file path was specified as an argument. If no organism is specified,
uniform codon usage is assumed by default.
"""
def codonUsage(tableFilePath=""):
    p = 1/float(64); # calculates probability of a codon according to uniform distribution
    codonFreqs = dict.fromkeys( list(geneticCode().keys()), 1/float(len(geneticCode())) ); # creates dictionary with codons as keys and a uniform probability of each key as values
    if len(tableFilePath) > 0: # if a file path was passed as an argument,
        txt = open(tableFilePath); # Access given file
        d = txt.read(); # Read file

        d = d.lstrip(); # Remove initial whitespace
        d = d.split("\n"); # Split into separate sequences

        totCodons = 0; # stores total number of codons
        for l in d[1::]: # iterates across file lines; starts at 1 to skip table header
            el = l.split(); # splits lines into elements separated by whitespace
            if len(el) > 0: # if line is not empty,
                codon = el[1]; # stores codon (second column)
                freq = float(el[2]); # stores codon's frequency (third column)
                codonFreqs[codon] = freq; # associates codon to frequency
                totCodons += freq; # adds number of codons to total

        txt.close(); # close file
        for c in codonFreqs: # for every codon,
            codonFreqs[c] = codonFreqs[c]/float(totCodons); # normalize over total number of codons

    return codonFreqs;


"""
Translates RNA to protein (takes rna string, returns peptide string).
"""
def translate(pSeq):
    seq = pSeq.upper(); # everything to uppercase
    code = geneticCode(); # Stores genetic code
    pep = ""; # stores peptide sequence
    i = 0; # position on RNA sequence
    while i < len(seq): # for every codon
        codon = seq[i:(i+3)]; # get current codon
        aa = code[codon]; # get aminoacid
        if (aa != "*"): # if not a stop codon
            pep = pep + aa; # add residue to peptide
            i = i+3; # next codon
        else: # if codon is stop
            i = len(seq); # stop loop
    return pep;


"""
Changes codons in dna sequence to synonyms according to given codon frequency
dictionary. Assumes frame starts in 0. Default codon frequency is equal
probabilities, resulting in codon scramble. Will select codons probabilistically
if sampling is True, will choose most likely codon otherwise.
"""
def optimizeCodons(pSeq, codonFreqs=codonUsage(), codonSampling=True):
    seq = transcribe(pSeq.upper()); # everything to uppercase
    code = geneticCode(); # Stores genetic code
    newSeq = ""; # stores new sequence
    i = 0; # position on DNA sequence
    while i < len(seq): # for every codon
        codon = seq[i:(i+3)]; # get current codon
        aa = code[codon]; # get aminoacid
        synonyms = []; # will store all synonyms
        for c in code: # iterates over genetic code
            if code[c] == code[codon]: # if codon is synonym,
                synonyms.append(c); # add to synonyms list

        totalProb = 0; # saves sum of likelihoods of all synonymous codons
        for syn in synonyms: # for every synonymous codon,
            totalProb += codonFreqs[syn]; # add this synonym's likelihood to total

        newCodon = synonyms[0]; # set first codon synonym as default new codon
        if codonSampling: # if not deterministic (choose codons probabilistically),
            r = random.random(); # random number between zero and 1
            j = 0; # stores index of codon selected
            cumProb = 0; # stores cumulative probability of belonging to all codons in indexes 0:i
            while cumProb < r: # will shift selected codon until cumulative probability surpasses random number
                cumProb += codonFreqs[synonyms[j]]/float(totalProb); # adds this codon's probabilty to the cumulative frequency
                j += 1; # advance counter

            newCodon = synonyms[j-1]; # change codon for chosen synonym

        else: # if deterministic mode (choose most likely codons),
            for c in synonyms: # iterate over synonyms
                if codonFreqs[c] > codonFreqs[newCodon]: # if this codon is more likely than current new codon,
                    newCodon = c; # set this codon as the new codon


        newSeq = newSeq + newCodon; # add codon to sequence
        i = i+3; # next codon

    newDNASeq = revTranscribe(newSeq); # reverse transcribe sequence
    return newDNASeq;


"""
Returns the reverse complement of a sequence
"""
def revComp(seq):
    revSeq = seq[::-1]; # reverse seq

    # replace nucleotides, conserve case
    revSeq = revSeq.replace("A", "1");
    revSeq = revSeq.replace("T", "2");
    revSeq = revSeq.replace("G", "3");
    revSeq = revSeq.replace("C", "4");

    revSeq = revSeq.replace("2", "A");
    revSeq = revSeq.replace("1", "T");
    revSeq = revSeq.replace("4", "G");
    revSeq = revSeq.replace("3", "C");

    revSeq = revSeq.replace("a", "1");
    revSeq = revSeq.replace("t", "2");
    revSeq = revSeq.replace("g", "3");
    revSeq = revSeq.replace("c", "4");

    revSeq = revSeq.replace("2", "a");
    revSeq = revSeq.replace("1", "t");
    revSeq = revSeq.replace("4", "g");
    revSeq = revSeq.replace("3", "c");

    return revSeq;

"""
Returns a list with all possible DNA sequences given an ambiguous DNA sequence.
"""
def ambiguousSeqs(seq):
    ambigCodes = {"R":"AG", "Y":"TC", "S":"CG", "W":"AT", "K":"GT", "M":"AC",
                  "B":"CGT", "D":"AGT", "H":"ATC", "V":"ACG", "N":"ATCG"};

    outSeqs = [""];
    for b in seq:
        if b in list(ambigCodes.keys()):
            newOutSeqs = [];
            for possibleB in ambigCodes[b]:
                newNewOutSeqs = [];
                for s in outSeqs:
                    newNewOutSeqs.append(s + possibleB);

                newOutSeqs = newOutSeqs + newNewOutSeqs;

            outSeqs = newOutSeqs;
        else:
            newOutSeqs = [];
            for s in outSeqs:
                newOutSeqs.append(s + b);

            outSeqs = newOutSeqs;


    return outSeqs;

# Auxiliary methods
"""
Checks if a sequence is hard to synthesize.
"""
def isTricky(seq):
    tricky = False;

    if findFirst(seq,"TATATATATATATATATATA") > -1: # if 10 TA repeats found,
        tricky = True; # it's tricky
    elif findFirst(seq,"GCGCGCGCGCGCGC") > -1: # if 7 GC repeats found,
        tricky = True; # it's tricky
    elif findFirst(seq,"AAAAAAAAAA") > -1: # if 10 A repeats found, (used to be 13)
        tricky = True; # it's tricky
    elif findFirst(seq,"TTTTTTTTTT") > -1: # if 10 T repeats found,
        tricky = True; # it's tricky
    elif findFirst(seq,"GGGGGGGGG") > -1: # if 9 G repeats found,
        tricky = True; # it's tricky
    elif findFirst(seq,"CCCCCCCCC") > -1: # if 9 C repeats found,
        tricky = True; # it's tricky

    return tricky;


"""
Returns a random hex color code.
"""
def randHex():
    chars = [0,1,2,3,4,5,6,7,8,9,'a','b','c','d','e','f']; # stores all possible values
    codeList = [str(chars[random.randint(0,len(chars)-1)]) for i in range(6)]; # generates six random characters from list
    hexCode = "#" + "".join(codeList); # generates hex code string
    return hexCode; # returns


"""
Returns a string with today'd date in format "06-JUN-2016"
"""
def todayDateStr():
    months = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]; # store month strings
    today = datetime.date.today(); # get today's date
    r = "-".join([str(today.day), months[today.month-1], str(today.year)]); # creates output string
    return r; # returns

# //
