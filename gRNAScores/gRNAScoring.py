# -*- coding: utf-8 -*-
"""
Created on Sat Jan 7 2017
Python functions used to score gRNAs according to Hsu et al. (2013) and Doench
et al. (2016) algorithms. Off-scoring function dependent on C++ binary
executable and text file with off-target gRNA database.
@author: Pablo CR
"""

import subprocess;
from py.BioUtils import *;
from gRNAScores.Rule_Set_2_scoring_v1.analysis.rs2_score_calculator import *; # Import Doench et al. (2016) on-target scoring module

model_file_1 = './gRNAScores/Rule_Set_2_scoring_v1/saved_models/V3_model_nopos.pickle'
model_file_2 = './gRNAScores/Rule_Set_2_scoring_v1/saved_models/V3_model_full.pickle'
with open(model_file_1, 'rb') as f:
    model1= pickle.load(f)

with open(model_file_2, 'rb') as f:
    model2= pickle.load(f)

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
def offTargetScore(gRNA, method, gListFilePath="gRNAScores/OffTarget_Scoring/gList.txt", binExecPath="gRNAScores/OffTarget_Scoring/offTargetScoring"):
    args = (binExecPath, gRNA, gListFilePath, method, "no"); # stores command to be passed to console. "no" specifies output format as just scores.
    popen = subprocess.Popen(args, stdout=subprocess.PIPE); # passes command to console
    popen.wait(); # waits?
    output = popen.stdout.read(); # get console output
    scores = [float(x) for x in output.replace('\n','').split(',')]; # parse outputs as floats
    scores[2] = int(scores[2]); # parse third output as integer
    scores[3] = scores[3] == 1; # and fourth as Boolean
    return scores;

'''
Receives on-target gRNA string (4-mer prefix, 20-mer gRNA, 3-mer PAM, 3-mer
suffix). Returns on-target score as per Doench et al. (2016). This method is
adapted from and accesses the authors' original code.
'''
def onTargetScore(pSeq,aa_cut=-1,per_peptide=-1):
    seq = pSeq.upper()
    if len(seq)!=30:
        print "Please enter a 30mer sequence."

    if (aa_cut == -1) or (per_peptide == -1):
        model = model1
    else:
        model = model2

    if seq[25:27] == 'GG':
        score = model_comparison.predict(seq, aa_cut, per_peptide, model=model)*100 # Normalized to 0-100 scale
        return score
    else:
        print >> sys.stderr, 'Calculates on-target scores for sgRNAs with NGG PAM only.'
