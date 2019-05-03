# -*- coding: utf-8 -*-
"""
Created on Jan 2017
Python functions used to score gRNAs according to Hsu et al. (2013) and Doench
et al. (2016) algorithms. These ended up being too slow and were implemented in
C++ (see gRNAScoring.py). Other functions in this package were used to build the
gRNA database and scoring matrix used by the C++ functions. The code in this
package is not used by GeneTargeterMethods during runtime, but is kept here just
in case.
@author: Pablo CR
"""
from py.utils.BioUtils import *;
from gRNAScores.Rule_Set_2_scoring_v1.analysis.rs2_score_calculator import *; # Import Doench et al. (2016) on-target scoring module
from gRNAScores.CFD_Scoring.cfd_score_calculator import *; # Import Doench et al. (2016) off-target scoring module
import subprocess;
import multiprocessing
import joblib as jl

mm_scores,pam_scores = get_mm_pam_scores_remote();

model_file_1 = './gRNAScores/Rule_Set_2_scoring_v1/saved_models/V3_model_nopos.pickle'
model_file_2 = './gRNAScores/Rule_Set_2_scoring_v1/saved_models/V3_model_full.pickle'
with open(model_file_1, 'rb') as f:
    model1= pickle.load(f)

with open(model_file_2, 'rb') as f:
    model2= pickle.load(f)

def exportCFDMatrix(scoreDict):
    baseCodes = {'A':0,'U':1,'G':2,'C':3};
    bases = baseCodes.keys();
    for i in range(20):
        for b1 in bases:
            for b2 in bases:
                key = 'r'+b1+':d'+revcom(b2)+','+str(i+1);
                if key not in scoreDict:
                    scoreDict[key] = 1.0;





    print len(scoreDict.keys());
    baseCodes['T'] = 1;
    scoreList = [];
    newKeysScores = {};
    for i,k in enumerate(scoreDict):
        position = int(k[6:])-1;
        baseOff = baseCodes[k[4]];
        baseOn = baseCodes[k[1]];
        newK = position*16+baseOn*4+baseOff;
        newKeysScores[newK] = scoreDict[k];

    keyList = sorted(newKeysScores.keys());
    outStr = "";
    for k in keyList:
        outStr = outStr + str(newKeysScores[k]) + "\n";

    output(outStr,"cfdScores.txt");
    return outStr.replace('\n',',')[:-1];


'''
'''
def buildPfalGRNADB(pathToGenome,enzyme="Cas9"):
    PfalGenome = loadFastas(pathToGenome);
    PfalGenome.pop('M76611 | organism=Plasmodium_falciparum_3D7 | version=2013-03-01 | length=5967 | SO=mitochondrial_chromosome',None);
    PfalGenome.pop('PFC10_API_IRAB | organism=Plasmodium_falciparum_3D7 | version=2013-03-01 | length=34242 | SO=apicoplast_chromosome',None);

    PfalGRNAs = [];
    if enzyme == "Cas9":
        PfalGRNAs = buildGRNADB(PfalGenome,"GG",[-21,2]);
        PfalGRNAs = PfalGRNAs + buildGRNADB(PfalGenome,"AG",[-21,2]);
    elif enzyme == "Cas12":
        PfalGRNAs = buildGRNADB(PfalGenome,"TTTA",[0,27]);
        PfalGRNAs = PfalGRNAs + buildGRNADB(PfalGenome,"TTTC",[0,27]);
        PfalGRNAs = PfalGRNAs + buildGRNADB(PfalGenome,"TTTG",[0,27]);
        PfalGRNAs = PfalGRNAs + buildGRNADB(PfalGenome,"CTTA",[0,27]);
        PfalGRNAs = PfalGRNAs + buildGRNADB(PfalGenome,"CTTC",[0,27]);
        PfalGRNAs = PfalGRNAs + buildGRNADB(PfalGenome,"CTTG",[0,27]);

    #guideListList = [];
    #for g in PfalGRNAs:
    #    guideListList.append(list(g.upper().replace('T','U')));

    return PfalGRNAs;#guideListList;

'''
Accepts a dictionary containing a multifasta file (see loadFastas in BioUtils)
Returns a list of all possible gRNAs within that multifasta file.
Each gRNA sequence in return list contains 20mer gRNA and NGG PAM sequence.
'''
def buildGRNADB(genome,PAM="GG",gRNAIndexes=[-21,2]):
    print "Start building";
    pamSeq = PAM.upper(); # should be derived from PAM, but replacing N would involve regex and I'm lazy. The Doench et al. (2016) gRNA scoring technologies only work with NGG PAMs anyway.
    gRNADB = []; # will store all possible gRNAs in genome
    for chrom in genome: # for every separate sequence within multifasta file,
        i = 0; # indexes through sequence
        chromSeq = genome[chrom].upper(); # chromosome sequence to uppercase
        while i < len(chromSeq): # iterate through all chromosome
            #print "Processing " + str(i) + "/" + str(len(genome[chrom])); # print chromosome
            nextI = chromSeq.find(pamSeq,i+1); # find next PAM downstream on plus strand
            altI = chromSeq.find(revComp(pamSeq),i+1); # find next PAM downstream on minus strand
            comp = False; # set sense to plus strand

            if altI < nextI: # if next PAM downstream on minus strand is further upstream than next one on plus strand,
                nextI = altI; # set next PAM downstream on minus strand
                comp = True; # set sense to complementary strand

            i = nextI; # set index to next PAM sequence
            if i == -1: # if no PAM found,
                break; # escape loop

            gRNA = chromSeq[i+gRNAIndexes[0]:i+gRNAIndexes[1]]; # store sequence gRNA-PAM
            if comp: # if on comp strand
                gRNA = revComp(chromSeq[i-gRNAIndexes[1]+len(PAM):i-gRNAIndexes[0]+len(PAM)]); # store sequence gRNA-PAM on comp strand

            gRNA = gRNA.upper(); # make gRNA uppercase
            gChk = gRNA.replace("A",""); # string used to check whether gRNA contains only ACTG characters
            gChk = gChk.replace("T","");
            gChk = gChk.replace("C","");
            gChk = gChk.replace("G","");
            if len(gRNA) == gRNAIndexes[1]-gRNAIndexes[0] and len(gChk) == 0: # if gRNA is right size and contains only accepted characters,
                gRNADB.append(gRNA); # add gRNA to DB


    return gRNADB;



###
'''

'''
def offScoreGRNADB(guideList):
    guideListList = [];
    gList = [];
    for g in guideList:
        guideListList.append(list(g.upper().replace('T','U')));
        gList.append(g.upper().replace('T','U'));

    scoredGRNAs = {};
    c = 0;
    for i in range(len(gList)):
        g = gList[i];
        scoredGRNAs[g] = offScore(g,guideListList);
        print str(float(c*100)/len(gList)) + "%";
        c += 1;

    return scoredGRNAs;


'''

'''
def offScoreHsu(pGRNA,gList,threshold=0):
    g = pGRNA.upper().replace('T','U');
    pam = g[-2:];
    cumScore = 0;
    scores = [];
    exactMatches = 0;
    for i,g2 in enumerate(gList):
        hit_score = calc_Hsu_score(g,g2);
        if hit_score > threshold: # cutoff determined by Haeussler et al. (2016) https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2
            cumScore = cumScore + hit_score;
            scores.append(hit_score);

    specificityScore = 100*(100.0/(100.0+cumScore)); # taken from http://crispr.mit.edu/about
    scores.sort(reverse=True);
    maxScore = scores[0]; # stores max score
    numScores = len(scores); # stores number of scores above threshold
    return [specificityScore,maxScore,numScores]; # return max value, make call based on that and score.


def offScoreCFD(pGRNA,guideListList,threshold=2.3):
    g = pGRNA.upper().replace('T','U');
    cumScore = 0;
    scores = [];
    exactMatches = 0;
    pam = g[-2:];
    for i in range(len(guideListList)):
        sg_list = guideListList[i][:-3].upper().replace('T','U');
        wt_list = list(g);
        hit_score = calc_cfd_basic(wt_list,sg_list,pam);
        if hit_score == 1:
            exactMatches += 1;

        if exactMatches > 1:
            exactMatches = 1;
        elif hit_score > threshold: # cutoff determined by Haeussler et al. (2016) https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2
            cumScore = cumScore + hit_score;
            scores.append(hit_score);

    specificityScore = 100*(100.0/(100.0+cumScore)); # taken from http://crispr.mit.edu/about
    scores.sort(reverse=True);
    maxScore = scores[0]; # stores max score
    numScores = len(scores); # stores number of scores above threshold
    return [specificityScore,maxScore,numScores]; # return max value, make call based on that and score. Explore Cutoffs: single hit: 0.5, cumulative: 0.25 (CFD) (75 percentile). CFD metrics: # above 0.2, 0.2 and 0.5 thresholds, max score.


'''

'''
def calc_cfd_basic(wt_list,s_list,pam):
    score = 1
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score*= mm_scores[key]
    score*=pam_scores[pam]
    return (score*100)


'''

'''
def calc_Hsu_score(gOn,gOff): # taken from http://crispr.mit.edu/about
    HsuMatrix = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583];
    score = 1;
    modifier = 1 - 0.8*(gOff[21] == 'A');

    sumD = 0;
    prevMM = 0;
    numMM = 0;
    for i in range(20):
        if gOn[i] != gOff[i]:
            score *= 1-HsuMatrix[i];
            sumD = (sumD + i - prevMM) * (prevMM != 0);
            prevMM = i;
            numMM += 1;


    if numMM > 1:
        d = sumD/float(numMM-1);
        modifier *= 1.0/((19.0-d)/19 * 4 + 1) * 1.0/(numMM**2);

    score *= modifier * 100;

    return score;

'''

'''
def onScore(pSeq,aa_cut=-1,per_peptide=-1):
    seq = pSeq.upper()
    if len(seq)!=30:
        print "Please enter a 30mer sequence."

    if (aa_cut == -1) or (per_peptide == -1):
        model = model1
    else:
        model = model2

    if seq[25:27] == 'GG':
        score = model_comparison.predict(seq, aa_cut, per_peptide, model=model)
        return score
    else:
        print >> sys.stderr, 'Calculates on-target scores for sgRNAs with NGG PAM only.'
'''
def scoreGRNA(g):
    gOffScoreHsu = offTargetScore(g,'hsu');
    gOffScoreCFD = offTargetScore(g,'cfd');
    gGCContent = gcContent(g[:len(g)-3]);
    return g +","+ str(gGCContent) +","+ str(gOffScoreCFD[0]) +","+ str(gOffScoreCFD[1]) +","+ str(gOffScoreCFD[2]) +","+ str(gOffScoreHsu[0]) +","+ str(gOffScoreHsu[1]) +","+ str(gOffScoreHsu[2]) +"\n";

def scoreGRNAs(gRNAs):
    outStr = "gRNA,GC Content,Total off-target score (Doench et al. 2016), Max off-target score (Doench et al. 2016), Number of off-target scores (Doench et al. 2016), Total off-target score (Hsu et al. 2013), Max off-target score (Hsu et al. 2013), Number of off-target scores (Hsu et al. 2013)\n";
    n_cores = multiprocessing.cpu_count()
    res = jl.Parallel(n_jobs=n_cores,verbose=8)(jl.delayed(scoreGRNA)(g) for g in gRNAs)
    outStr = outStr + "".join(res)

    output(outStr,"gRNAOutput.csv");
'''

def scoreGRNAs(gRNAs):
    outStr = "gRNA,GC Content,Total off-target score (Doench et al. 2016), Max off-target score (Doench et al. 2016), Number of off-target scores (Doench et al. 2016), Total off-target score (Hsu et al. 2013), Max off-target score (Hsu et al. 2013), Number of off-target scores (Hsu et al. 2013)\n";
    for i,g in enumerate(gRNAs):
        gOffScoreHsu = offTargetScore(g,'hsu');
        gOffScoreCFD = offTargetScore(g,'cfd');
        gGCContent = gcContent(g[:len(g)-3]);
        outStr = outStr + g +","+ str(gGCContent) +","+ str(gOffScoreCFD[0]) +","+ str(gOffScoreCFD[1]) +","+ str(gOffScoreCFD[2]) +","+ str(gOffScoreHsu[0]) +","+ str(gOffScoreHsu[1]) +","+ str(gOffScoreHsu[2]) +"\n";
        print i;

    output(outStr,"gRNAOutput.csv");

def timeFunc(gList):
    import sys
    import time

    default_timer = time.time;
    start = default_timer();
    offTargetScore("TTGTGTTCTCCATATATCGATGG","hsu");
    finish = default_timer();
    elapsed = (finish - start);

    return elapsed;



def offTargetScore(gRNA, method, gListFilePath="gRNAScores/OffTarget_Scoring/gListCas9.txt",pamType="NGG",mm='4'):
    if pamType=="NGG":
        pam = gRNA[20:23]

    elif pamType=="TTTV":
        pam = gRNA[0:4]

    args = ("gRNAScores/OffTarget_Scoring/offTargetScoringBinLinux", gRNA, gListFilePath, method, pam, pamType, mm, "no");
    popen = subprocess.Popen(args, stdout=subprocess.PIPE);
    popen.wait();
    output = popen.stdout.read();
    scores = [float(x) for x in output.replace('\n','').split(',')];
    scores[2] = int(scores[2]);
    scores[3] = scores[3] == 1;
    return scores;

'''
from gRNAScores.gRNADBBuilder import *
gList = buildPfalGRNADB("./input/genomes/PlasmoDB-28_Pfalciparum3D7_Genome.fasta");
g="ATAACATCTGCTTTAAATTCTGG"
offScoreHsu(g,gList);
n,s=offScore(g,guideListList,calc_cfd_basic);

outStr = "\n".join(gList);
output(outStr,'nuclearGListCas9.txt');

from gRNAScores.gRNADBBuilder import *
guideListList = buildPfalGRNADB("./input/genomes/PlasmoDB-28_Pfalciparum3D7_Genome.fasta");
PfalGenome = loadFastas("./input/genomes/PlasmoDB-28_Pfalciparum3D7_Genome.fasta");
PfalGRNAs = buildGRNADB(PfalGenome,"GG");
import numpy as np
subset = np.random.choice(PfalGRNAs,100);
HsuScores = []
CFDScores = []
for i,g in enumerate(subset):
    #nH,sH=offScore(g,guideListList,calc_Hsu_score);
    nC,sC=offScore(g,guideListList,calc_cfd_basic);
    #HsuScores.append(nH)
    CFDScores.append(nC)
    print i

s.sort(reverse=True)
s[0]
x = list("ATAACATCTGCTTTAAATTC".replace("T","U"))

print offScore(g,guideListList);
print onScore("tgaaaatgatttaatgacacccacaggaaa");

off = offScoreGRNADB(PfalGRNAs)

print offScore("acatccatagggaaaaaaga",PfalGRNAs);
scoreGRNADB(PfalGRNAs);
output = open('PfalGRNAs.pkl', 'wb');
pickle.dump(PfalGRNAs, output);
with open('PfalGRNAs.csv', 'wb') as f:  # Just use 'w' mode in 3.x
    w = DictWriter(f, PfalGRNAs.keys())
    w.writeheader()
    w.writerow(PfalGRNAs)


seq="atatgtagctgtagaaaaagtagaggcatcatataaacaatctttatttattagccaaatttttgtatttggatattcctatgaatctgttctcgtttgtgttatttgtccatcgacagattccatcgatatatggagaacacaaaagaaaatcaaagcaactgatgaagaagtaattaaattaccagaatttaaagcagatgttattaatgatttaacatccatagggaaaaaagacggactcaaaggatttgagcaaattaaggatattcatttcactcttgaggcatttactattgaaaatgatttaatgacacccacaggaaaaattaaaagacatgaagctaagaagagatttaaaaaggaaattgatgagatgtatgaaaagttgaagcaatagaacataaaggataatatgcacaaatatatatatgtatatatgtatatatatatatatatatatatatatatatatatatatatatatgtatatgtttgtatttatatatttatagaagaaacaaaaagaagtgaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaccgtgattcctgtttttttttttttctctttttca"
from gRNAScores.gRNADBBuilder import *
guideListList = buildPfalGRNADB("./input/genomes/PlasmoDB-28_Pfalciparum3D7_Genome.fasta");
dSeq = {1:seq}
seqGRNAs=buildGRNADB30mers(dSeq)
scoreGRNAs(seqGRNAs,guideListList)
#from gRNAScores.gRNADBBuilder import *
#exportCFDMatrix(mm_scores)
#offTargetScore("TTGTGTTCTCCATATATCGATGG","cfd");
outStr = "\n".join(gListCpf1);
output(outStr, "gListCpf1.txt");

from gRNAScores.gRNADBBuilder import *
gs = buildPfalGRNADB("all_genes.fasta")
gs = list(dict.fromkeys(gs))
scoreGRNAs(gs)

'''
