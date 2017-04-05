# -*- coding: utf-8 -*-
"""
Created on Sat Jan 7 2017
Python functions used to score gRNAs according to Hsu et al. (2013) and Doench
et al. (2016) algorithms. Off-scoring function dependent on C++ binary
executable and text file with off-target gRNA database.
@author: Pablo CR
"""

import subprocess;
import sys;
import math;
import pickle;
import numpy as np;
from py.BioUtils import *;
import gRNAScores.Rule_Set_2_scoring_v1.analysis.model_comparison as RS2; # Import Doench et al. (2016) old on-target scoring module
from gRNAScores.CFD_Scoring.cfd_score_calculator import *; # Import Doench et al. (2016) off-target scoring module
import gRNAScores.azimuth.model_comparison as azimuth; # Import Doench et al. (2016) new on-target scoring module

# Get models used by old Rule Set 2 scores
model_file_1 = './gRNAScores/Rule_Set_2_scoring_v1/saved_models/V3_model_nopos.pickle'
model_file_2 = './gRNAScores/Rule_Set_2_scoring_v1/saved_models/V3_model_full.pickle'
with open(model_file_1, 'rb') as f:
    model1= pickle.load(f)

with open(model_file_2, 'rb') as f:
    model2= pickle.load(f)

# Get models used by CFD scores
mm_scores,pam_scores = get_mm_pam_scores_remote();

cindelCoefficients = { # dictionary with keys=features and values=coefficients used in CINDEL logistic regression (Kim et al. 2017)
    "(Intercept)":-1.830804463, "A":0.101900819, "A9":0.037608778,
    "AA13":0.167340348, "AA25":0.364227683, "AA6":-0.37092459,
    "AC":0.038389043, "AC6":0.083159983, "AC9":0.459109261,
    "AG16":-0.257994467, "AT22":0.469569195, "C26":0.166195809,
    "CA14":0.11154741,"CA5":-0.317283912, "CC":-0.24886262,
    "CC23":0.130179283, "CC4":1.065828073, "CG":0.006835173,
    "CG24":-0.186173348, "CG6":-0.231597251, "CG7":0.275263279,
    "CT13":0.025086518, "CT19":0.164219823, "CT24":0.083624311,
    "CT6":0.120799685, "CT7":-0.336424167, "CT9":-0.094751841,
    "Free_energy":0.11018728, "G19":-0.076228205, "G22":-0.066751042,
    "G23":-0.102884761, "G26":-0.043253074, "G4":0.489354734,
    "GA10":0.059114363, "GA18":0.077528126, "GA21":0.216211044,
    "GA4":0.077064537, "GC20":-0.123286737, "GG21":-0.124324225,
    "GG22":-0.358412408, "GG5":-0.058789262, "GG8":0.126113307,
    "GT23":-0.155677567, "T21":0.048550229, "T4":-0.781196013,
    "T9":-0.028800337, "TA":0.037868724, "TA23":0.444556704,
    "TA24":0.274164613, "TA7":0.391209275, "TC15":0.123925716,
    "TC25":0.194045073, "TG10":-0.190097449, "TG16":0.036768788,
    "TG6":0.12530035, "TT10":-0.173588498, "TT12":-0.035141319,
    "TT16":-0.127599558
};

'''
Receives on-target gRNA string (20-mer gRNA, 3-mer PAM) and string specifying
scoring method to be used ("hsu" for Hsu et al., 2013; "cfd" for Cutting
Frequency Determination, Doench et al., 2016). Calls C++ binary executable to do
the actual scoring. Executable location in binExecPath. Executable depends on
text file with off-target gRNA database, found in gListFilePath. Returns array
containing total score (0-100 scale), score of highest-scoring individual hit
(0-100 scale), number of hits scoring higher than threshold (threshold specified
in C++ methods), and Boolean value specifying if the on-target gRNA was found in
the off-target database or not.
'''
def offTargetScore(gRNA, method, enzyme, pamSeq, pamType, maxNumMisMatches=4, gListFilePath="", binExecPath="gRNAScores/OffTarget_Scoring/offTargetScoringBin"):
    platform = "Linux"; # default to linux platform
    if sys.platform == "darwin": # if running on mac yosemite,
        platform = "OSX"; # set new osx to access alternate binary file

    if gListFilePath == "": # if no argument was passed to gListFilePath
        if enzyme == "Cas9": # if enzyme is Cas9
            gListFilePath = "gRNAScores/OffTarget_Scoring/gListCas9.txt"; # use this list
        elif enzyme == "Cpf1": # if enzyme is Cpf1
            gListFilePath = "gRNAScores/OffTarget_Scoring/gListCpf1.txt"; # use cpf1 list


    args = (binExecPath+platform, gRNA, gListFilePath, method, pamSeq, pamType, str(maxNumMisMatches), "no"); # stores command to be passed to console. "no" specifies output format as just scores.
    popen = subprocess.Popen(args, stdout=subprocess.PIPE); # passes command to console
    popen.wait(); # waits?
    output = popen.stdout.read(); # get console output
    scores = [float(x) for x in output.replace('\n','').split(',')]; # parse outputs as floats
    scores[2] = int(scores[2]); # parse third output as integer
    scores[3] = scores[3] == 1; # and fourth as Boolean
    return scores;

'''
Returns on-target score based on desired method.
'''
def onTargetScore(pSeq,onTargetMethod):
    score = 0; # stores return score (0-100% format)
    if onTargetMethod == "azimuth": # if using Azimuth (Doench et al., 2016)
        score = onTargetScoreAzimuth(pSeq); # use it
    elif onTargetMethod == "ruleset2": # if using Rule Set 2 (Doench et al., 2016)
        score = onTargetScoreRS2(pSeq); # use it
    elif onTargetMethod == "cindel": # if using CINDEL (Kim et al., 2017)
        score = onTargetScoreCINDEL(pSeq[0:27]); # use it

    return score;

'''
Receives on-target gRNA string (4-mer prefix, 20-mer gRNA, 3-mer PAM, 3-mer
suffix). Returns on-target score as per Doench et al. (2016). This method is
adapted from and accesses the authors' original code.
'''
def onTargetScoreRS2(pSeq,aa_cut=-1,per_peptide=-1):
    seq = pSeq.upper()
    if len(seq)!=30:
        print "Please enter a 30mer sequence."

    if (aa_cut == -1) or (per_peptide == -1):
        model = model1
    else:
        model = model2

    if seq[25:27] == 'GG':
        score = RS2.predict(seq, aa_cut, per_peptide, model=model)*100 # Normalized to 0-100 scale
        return score
    else:
        print >> sys.stderr, 'Calculates on-target scores for sgRNAs with NGG PAM only.'

'''
Receives on-target gRNA string (4-mer prefix, 20-mer gRNA, 3-mer PAM, 3-mer
suffix). Returns on-target score using Azimuth software
https://github.com/MicrosoftResearch/Azimuth
based on Doench et al. (2016). This method is adapted from the Doench original
code and accesses the Azimuth source provided by MicrosoftResearch on GitHub.
'''
def onTargetScoreAzimuth(pSeq,aa_cut=-1,per_peptide=-1):
    seq = pSeq.upper();
    if len(seq)!=30:
        print "Please enter a 30mer sequence."

    if (aa_cut == -1) or (per_peptide == -1):
        model = model1
    else:
        model = model2

    if seq[25:27] == 'GG':
        score = azimuth.predict( np.array([pSeq]),np.array([aa_cut]),np.array([per_peptide]) )[0] * 100 # Normalized to 0-100 scale
        return score
    else:
        print >> sys.stderr, 'Calculates on-target scores for sgRNAs with NGG PAM only.'


'''
Receives on-target gRNA string (4-mer PAM, 23-mer gRNA). Returns on-target score
reverse-engineering CINDEL (http://big.hanyang.ac.kr/cindel/). Based off of the
Kim et al. (2017) paper, uses a linear model with coefficients obtained from the
supplementary material.
'''
def onTargetScoreCINDEL(pSeq):
    seq = pSeq.upper();
    linearProduct = 0; # score will be linear combination of coefficients vector and features
    if len(seq)!=27: # if not right size,
        print "Please enter a 27mer sequence." # say so
    else: # if right size,
        for feature in cindelCoefficients: # for every feature coefficient in dictionary
            if feature[len(feature)-1].isdigit(): # if last character in feature is a digit
                i = 0; # will store position of first occurrence of a digit in feature string
                for char in feature: # loops over characters in string
                    if char.isdigit(): # if the character is a digit,
                        break; # escape loop

                    i += 1; # advance indexer if still in loop


                featureSeq = feature[0:i]; # stores nucleotide sequence of feature to be searched
                featureIndexes = [int(feature[i:len(feature)]),int(feature[i:len(feature)]) + len(featureSeq)]; # stores indexes of feature to be searched
                featurePresent = (featureSeq == seq[featureIndexes[0]:featureIndexes[1]]); # Boolean feature presence value
                linearProduct += cindelCoefficients[feature]*featurePresent; # add coefficient to score if feature present

            elif feature == "Free_energy": # if feature is actually intercept
                dG = freeEnergy(seq); # calculates RNA self-folding Gibbs free energy
                linearProduct += cindelCoefficients[feature]*dG; # add product of free energy and coefficient to score
            elif feature == "(Intercept)": # if feature is actually intercept
                linearProduct += cindelCoefficients[feature]*1; # add to score
            else: # if non of the above, assume feature is a nucleotide sequence whose number of occurrences will be counted
                featureCount = seq.count(feature); # counts occurrences of feature within sequence
                linearProduct += cindelCoefficients[feature]*featureCount; # add product of occurrences and coefficient to score




    score = math.exp(linearProduct)/(1.0+math.exp(linearProduct)) * 100; # calculate probability according to logistic regression model. Returned as percent, not decimal.
    return score;

'''
Obtains Gibbs free energy (kcal/mol) for a given sequence using ViennaRNA's
RNAfold binary.
'''
def freeEnergy(seq, binExecPath="./gRNAScores/RNAfold/RNAfoldBin"):
    platform = "Linux"; # default to linux platform
    if sys.platform == "darwin": # if running on mac yosemite,
        platform = "OSX"; # set new osx to access alternate binary file

    args = (binExecPath+platform, "--noPS"); # stores command to be passed to console.

    popen = subprocess.Popen(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash'); # passes command to console
    output = popen.communicate(input=seq)[0]; # passes sequence to console dialog

    numberIndex = output.find("-"); # finds minus sign to search for free energy
    if numberIndex < 0: # if no minus sign found,
        numberIndex = output.find("0.") - 1; # confirms 0.00 kcal/mol, sets index to appropriate position before number

    score = float(output[numberIndex:numberIndex+4]); # extracts free energy (kcal/mol) from output

    return score;


'''
Calculates individual hit score of a given on-target gRNA and a given off-target
gRNA as per Zhang lab procedure (Hsu et al., 2013; http://crispr.mit.edu/about).
The mean pairwise distance d is calculated as specified in the Haeussler et al.
(2016) algorithm, as opposed to the one proposed in this Zhang Lab online tool
forum thread (https://groups.google.com/forum/#!searchin/crispr/algorithm|sort:relevance/crispr/fkhX7FU3r-I/9Nc14v_j3XgJ)
Both methods are numerically equivalent. NAG weight taken from Hsu et al. (2013)
paper ("approximately five times less than NAG").
Both gRNAs must be same-sense 20-mers, with the PAM type and sequence specified
as parameters.
'''
def pairScoreHsu(gOn,gOff,pamSeq,pamType):
    HsuMatrix = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]; # Hsu score matrix
    score = 1; # stores score
    modifier = 1; # stores multiplicative score modifier
    if pamType == "NGG": # if PAM is NGG
        modifier = 1# - (1-0.2)*(pamSeq[1] == 'A'); # stores factor modifying score in case of NAG PAM or multiple mismatches (Hsu et al. 2013)
    elif pamType == "TTTV": # if PAM is TTTV
        modifier = 1# - (1-0.2)*(pamSeq[0] == 'C'); # stores factor modifying score in case of NAG PAM or multiple mismatches (estimated from Kim et al., 2017)

    sumD = 0; # stores sum of distances between mismatches
    prevMM = 0; # stores number of mismatches
    numMM = 0; # stores previous mismatch's index, used for calculating distance between mismatches
    for i in range(20): # iterate through first 20 positions
        if gOn[i] != gOff[i]: # if there is a mismatch at this position,
            score *= 1-HsuMatrix[i]; # multiply score by the mismatch weight at this position from the weight matrix
            sumD = (sumD + i - prevMM) * (prevMM != 0); # add the distance between this mismatch and previous mismatch to the sum of distances. Note that the first mismatch does not contribute to sum.
            prevMM = i; # set this index as the new previous mismatch index.
            numMM += 1; # advance mismatch counter


    if numMM > 1: # if there are multiple mismatches,
        d = sumD/float(numMM-1); # calculate mean pairwise distance based on sum and number of distances
        modifier *= 1.0/((19.0-d)/19 * 4 + 1) * 1.0/(numMM**2); # modify modifier according to Zhang lab algorithm (http://crispr.mit.edu/about)

    score *= modifier * 100; # calculate final hit score, set scale to 0-100

    return score;


'''
Calculates individual hit score of a given on-target gRNA and a given off-target
gRNA as per CFD scores (Doench et al., 2016). Code based off of original Doench
et al. (2016) code in supplementary material.
Both gRNAs must be same-sense 20-mers, with the PAM type and sequence specified
as parameters.
'''
def pairScoreCFD(gOn,gOff,pamSeq,pamType):
    wt_list = list(gOn.upper()); # convert strings to lists
    s_list = list(gOff.upper());
    pam = pamSeq[1:3].upper();
    score = 1;
    for i,sl in enumerate(s_list): # for every base in the off-target gRNA,
        if wt_list[i] == sl: # if bases are equal,
            score *= 1; # don't modify score
        elif i < 20: # if not,
            key = 'r'+transcribe(wt_list[i])+':d'+revcom(sl)+','+str(i+1); # construct a key to access score penalty at this position for this mismatch
            score *= mm_scores[key]; # modify score


    if pamType == "NGG": # if PAM is NGG form
        score *= pam_scores[pam]; # use CFD module to modify score
    elif pamType == "TTTV": # if PAM is TTTV form
        score *= 1 - (1-0.2)*(pamSeq[0] == 'C'); # modify score in case of CTTV PAM (estimated from Kim et al., 2017)

    return (score*100);
