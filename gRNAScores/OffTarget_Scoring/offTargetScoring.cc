
/*/ --- gRNA Off-target Scoring --- /*/

// C++ code to calculate off-target specificity scores for a given gRNA, based
// on a given list of off-target gRNAs. The on-target gRNA string and off-target
// gRNA list file path string must be passed as parameters, as well as a third
// string specifying the method used for scoring off-targets. Off-targets can
// have NGG or NAG PAM sequences.
//
// Individual off-target scoring is set either with parameter "hsu" for the Hsu
// et al. (2013) method, which takes into account mismatch position but not
// nature, or "cfd" for Cutting Frequency Determination as defined by Doench et
// al. (2016), which takes into account the nucleotide identities in each
// mismatch. Individual hits are scored on a 0-100 basis, with 100% being equal
// frequency to the on-target gRNA.
//
// The scoring parameters were developed from Cas9 data, however, they can be
// used for Cpf1 while data on Cpf1 is published. It's recommended to use Hsu
// scores when evaluating Cpf1 cut sites.
//
// Aggregate score is calculated as 100 * 100/(100 + sum of individual hits),
// resulting in a 0-100 scale as well. If the on-target gRNA is found in the
// off-target list, it is omitted from the calculation.
//
// Code was based on the information provided by the Zhang group (http://crispr.mit.edu/about and https://groups.google.com/forum/#!searchin/crispr/algorithm|sort:relevance/crispr/fkhX7FU3r-I/9Nc14v_j3XgJ)
// for Hsu calculations, and from the source provided by Doench et al. (2016)
// in the case of the CFD method. The algorithm was changed to optimize speed,
// besides being implemented in C++ as opposed to Python.
//
// --Pablo CR, pablocarderam@gmail.com
//
// Compilation: g++ -o offTargetScoringBin offTargetScoring.cc
//
// Usage: ./offTargetScoringBin "TTGTGTTCTCCATATATCGA" "gListCas9.txt" "cfd" "Cas9" "TGG" "NGG" "yes"
//

/*/ ---  --- /*/

// Imports
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

// Constants
const double hsuMatrix[20] = {0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583}; // obtained from
const double cfdMatrix[320] = {1.0,1.0,0.857142857,1.0,1.0,1.0,0.857142857,0.956521739,1.0,0.9,0.714285714,1.0,1.0,1.0,1.0,0.913043478,0.727272727,1.0,0.785714286,0.8,1.0,0.846153846,0.857142857,0.84,0.636363636,0.846153846,0.692307692,1.0,0.909090909,0.727272727,1.0,0.695652174,0.705882353,1.0,0.428571429,0.611111111,1.0,0.714285714,0.428571429,0.5,0.5,0.75,0.384615385,1.0,0.6875,0.866666667,1.0,0.5,0.636363636,1.0,0.352941176,0.625,1.0,0.476190476,0.647058824,0.625,0.363636364,0.9,0.529411765,1.0,0.8,0.842105263,1.0,0.5,0.363636364,1.0,0.5,0.72,1.0,0.5,1.0,0.64,0.3,0.866666667,0.785714286,1.0,0.636363636,0.571428571,1.0,0.6,0.714285714,1.0,0.454545455,0.714285714,1.0,0.866666667,0.909090909,0.571428571,0.666666667,1.0,0.681818182,1.0,0.928571429,0.928571429,1.0,0.5,0.4375,1.0,0.4375,0.705882353,1.0,0.875,0.6875,0.588235294,0.571428571,1.0,0.6875,1.0,0.8125,0.75,1.0,0.470588235,0.428571429,1.0,0.428571429,0.733333333,1.0,0.8,1.0,0.733333333,0.625,1.0,0.615384615,1.0,0.875,0.65,1.0,0.642857143,0.6,1.0,0.571428571,0.666666667,1.0,0.928571429,0.923076923,0.619047619,0.533333333,0.642857143,0.538461538,1.0,0.875,0.857142857,1.0,0.619047619,0.882352941,1.0,0.333333333,0.555555556,1.0,0.857142857,0.533333333,0.5,0.8125,0.933333333,0.4,1.0,0.941176471,0.866666667,1.0,0.388888889,0.307692308,1.0,0.4,0.65,1.0,0.75,0.666666667,0.4,0.384615385,1.0,0.428571429,1.0,0.307692308,0.75,1.0,0.25,0.333333333,1.0,0.263157895,0.722222222,1.0,0.8,0.947368421,0.5,0.384615385,0.933333333,0.529411765,1.0,0.538461538,0.714285714,1.0,0.444444444,0.3,1.0,0.210526316,0.652173913,1.0,0.692307692,0.789473684,0.260869565,0.3,0.923076923,0.421052632,1.0,0.7,0.384615385,1.0,0.136363636,0.533333333,1.0,0.214285714,0.466666667,1.0,0.619047619,0.285714286,0.0,0.266666667,0.75,0.428571429,1.0,0.733333333,0.35,1.0,0.0,0.2,1.0,0.272727273,0.65,1.0,0.578947368,0.272727273,0.05,0.142857143,0.941176471,0.272727273,1.0,0.066666667,0.222222222,1.0,0.05,0.0,1.0,0.0,0.192307692,1.0,0.909090909,0.666666667,0.346153846,0.0,1.0,0.0,1.0,0.307692308,1.0,1.0,0.153846154,0.133333333,1.0,0.176470588,0.176470588,1.0,0.533333333,0.705882353,0.117647059,0.25,0.933333333,0.235294118,1.0,0.466666667,0.466666667,1.0,0.058823529,0.5,1.0,0.19047619,0.4,1.0,0.666666667,0.428571429,0.333333333,0.666666667,0.692307692,0.476190476,1.0,0.642857143,0.538461538,1.0,0.133333333,0.538461538,1.0,0.206896552,0.375,1.0,0.285714286,0.275862069,0.25,0.666666667,0.714285714,0.448275862,1.0,0.461538462,0.428571429,1.0,0.125,0.6,1.0,0.227272727,0.764705882,1.0,0.5625,0.090909091,0.176470588,0.7,0.9375,0.428571429,1.0,0.3,0.5,1.0,0.058823529};

// Structs
struct Score { // struct to store multiple return values of offScore function
    double totalScore; // total score of on-target gRNA
    double maxScore; // score of highest-scoring off-target hit for this on-target gRNA
    int numHits; // number of hits scoring above threshold for this on-target gRNA
    bool foundInList; // whether or not the on-target sequence was found in the off-target list
};

// Function declarations
double pairScoreHsu(string gOn, string gOff, string pamSeq, string pamType);
double pairScoreCFD(string gOn, string gOff, string pamSeq, string pamType);
Score offScore(string gRNA, string gListFilePath, string pamSeq, string method, string pamType);

// Main
int main (int argc, char* argv[]) { // receive arguments
    string gRNA = argv[1]; // on-target gRNA
    int i = 0; // index through gRNA
    while (gRNA[i]) { // loop through gRNA making uppercase
        gRNA[i] = toupper(gRNA[i]);
        i++;
    }
    string gListFilePath = argv[2]; // file path to off-target gRNA list
    string method = argv[3]; // string defining method
    string pamSeq = argv[4]; // PAM sequence of gRNA
    string pamType = argv[5]; // type of PAM sequence being used
    string explainOutput = "yes"; // stores whether or not output will be printed with format explanation. Default is it will.
    if (argv[6]) { // if there is a sixth argument in call,
        explainOutput = argv[6]; // set the output parameter to it.
    }

    Score score = offScore(gRNA,gListFilePath,pamSeq,method,pamType); // scoring function call

    if (explainOutput == "yes") {
        printf( "Total Score: %f, Max. Hit Score: %f, Num. Hits: %d, gRNA found in DB: %d\n", score.totalScore, score.maxScore, score.numHits, score.foundInList ); // output format
    }
    else {
        printf( "%f, %f, %d, %d\n", score.totalScore, score.maxScore, score.numHits, score.foundInList ); // output score
    }


    return 0; // main function requirement
}

/*
Receives on-target gRNA as a string, filepath string for off-target gRNA list.
Returns off-target score.
Will read off-target list every time function is called because function reads
file line by line and performs comparison without storing whole list in memory,
since the list file is too large to load.
*/
Score offScore(string gRNA, string gListFilePath, string pamSeq, string method, string pamType) {
    double threshold = 0; // sets threshold under which scores will be ignored.
    double (*pairScore) (string gOn, string gOff, string pamSeq, string pamType); // function pointer for individual hit scoring
    if (method == "hsu") { // if Hsu method is chosen,
        pairScore = &pairScoreHsu; // set scoring function accordingly.
        threshold = 0.05; // threshold based on Benchling's precision in CRISPR guide designing process, or at least what it shows us
    }
    else if (method == "cfd") {// if CFD method is chosen,
        pairScore = &pairScoreCFD; // set scoring function accordingly.
        threshold = 2.3; // cutoff determined by Haeussler et al. (2016) https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2
    }

    double singleScore = 0; // will store individual hit score
    double cumScore = 0; // will store sum of individual hit scores
    double maxScore = 0; // stores highest score
    int numHits = 0; // stores number of individual hit scores that are above threshold
    bool foundMatch = false; // keeps track of whether orr not an exact match has been found

    string offG; // will store off-target gRNA (line in file)
    ifstream gListFile (gListFilePath.c_str()); // file stream to off-target gRNA list
    if (gListFile.is_open()) { // if file is open,
        while (getline (gListFile,offG)) { // loop while there is a new gRNA to compare (new line in file). Stores next line in file (next off-target gRNA) in variable offG.
            singleScore = pairScore(gRNA,offG,pamSeq,pamType); // score on-target gRNA with this particular off-target gRNA
            if (singleScore == 100 && !foundMatch) { // if an exact match is found and is the first one,
                foundMatch = true;
            }
            else if (singleScore > threshold) { // if not first exact match and this score is higher than the threshold,
                cumScore = cumScore + singleScore; // add it to the cumulative score
                if (singleScore > maxScore) { // if this score is greater than current max,
                    maxScore = singleScore; // set as new max
                }
                ++numHits; // advance single hit counter
            }
        }
        gListFile.close(); // close file
    }
    else cout << "Unable to open file.\n"; // if file couldn't be opened, say so.

    //std::cout << "Max score: " << maxScore << '\n';
    double totalScore = 100.0 * 100.0/(100.0 + cumScore); // calculate aggregate score based on cumulative score

    Score output = {totalScore,maxScore,numHits,foundMatch}; // stores output variables for returning
    return output; // return aggregate score
}

/*
Calculates individual hit score of a given on-target gRNA and a given off-target
gRNA as per Zhang lab procedure (Hsu et al., 2013; http://crispr.mit.edu/about).
The mean pairwise distance d is calculated as specified in the Haeussler et al.
(2016) algorithm, as opposed to the one proposed in this Zhang Lab online tool
forum thread (https://groups.google.com/forum/#!searchin/crispr/algorithm|sort:relevance/crispr/fkhX7FU3r-I/9Nc14v_j3XgJ)
Both methods are numerically equivalent. NAG weight taken from Hsu et al. (2013)
paper ("approximately five times less than NAG").
Both gRNAs must be same-sense 23-mers, with the 3-mer (NGG or NAG) PAM sequence at the end.
*/
double pairScoreHsu(string gOn, string gOff, string pamSeq, string pamType) {
		double score = 1; // stores hit score
    double modifier = 1; // modifier due to PAM sequence
    if (pamType == "NGG") { // if PAM is NGG
        modifier = 1 - (1-0.2)*(pamSeq[1] == 'A'); // stores factor modifying score in case of NAG PAM or multiple mismatches (Hsu et al. 2013)
    }
    else if (pamType == "TTTV") { // if PAM is TTTV
        modifier = 1 - (1-0.2)*(pamSeq[0] == 'C'); // stores factor modifying score in case of NAG PAM or multiple mismatches (estimated from Kim et al., 2017)
    }

		int sumD = 0; // stores sum of distances between mismatches
	  int numMM = 0; // stores number of mismatches
    int prevMM = 0; // stores previous mismatch's index, used for calculating distance between mismatches
		for (size_t i = 0; i < 20; i++) { // iterate through first 20 positions
				if (gOn[i] != gOff[i]) { // if there is a mismatch at this position,
						score = score * (1.0-hsuMatrix[i]); // multiply score by the mismatch weight at this position from the weight matrix
            sumD = (sumD + i - prevMM) * (prevMM != 0); // add the distance between this mismatch and previous mismatch to the sum of distances. Note that the first mismatch does not contribute to sum.
            prevMM = i; // set this index as the new previous mismatch index.
			      ++numMM; // advance mismatch counter
				}
		}

		if (numMM > 1) { // if there are multiple mismatches,
        double d = (double)(sumD)/(double)(numMM-1); // calculate mean pairwise distance based on sum and number of distances
        modifier = modifier * 1.0/((19.0-d)/19.0 * 4.0 + 1.0) * 1.0/pow((double)numMM, 2.0); // modify modifier according to Zhang lab algorithm (http://crispr.mit.edu/about)
		}

    score = score * modifier * 100.0; // calculate final hit score, set scale to 0-100

		return score;
}

/*
Calculates individual hit score of a given on-target gRNA and a given off-target
gRNA as per Doench et al. (2016) Cutting Frequency Determination scale.
The score weights matrix for each position-mismatch type was obtained from the
original Doench et al. source code, and was completed to include all possible
bp matches (perfect matches have a weight of 1), resulting in 320 weights (20
positions * 4 possible bp in on-target position * 4 possible bp in off-target
position). Weights were sorted by position, on-target base (A>T>G>C), and
off-target base (A>T>G>C), in that order. In the new algorithm below, weights
are accesed from this list based on their order, instead of the Python
dictionary originally used. Off-target bases are reverse-complemented, as in
original algorithm! Only NAG PAM weights are included (other than NGG,
weighted as 1), to keep database size down. Score is on a 0-100 scale instead of
the original 0-1 scale used by the authors.
Both gRNAs must be same-sense 23-mers, with the 3-mer PAM sequence at the end.
*/
double pairScoreCFD(string gOn, string gOff, string pamSeq, string pamType) {
    double score = 1; // hit score
    double modifier = 1; // modifier for PAM mismatches
    if (pamType == "NGG") { // if PAM is NGG
        modifier = 1 - (1-0.25925925899999996)*(pamSeq[1] == 'A'); // modifier for NAG PAM (Doench et al., 2016)
    }
    else if (pamType == "TTTV") { // if PAM is TTTV
        modifier = 1 - (1-0.2)*(pamSeq[0] == 'C'); // stores factor modifying score in case of NAG PAM or multiple mismatches (estimated from Kim et al., 2017)
    }
    int onBase = -1; // stores integer representation of on-target base
    int offBase = -1; // stores integer representation of off-target base

		for (size_t i = 0; i < 20; i++) { // for the first 20 positions in gRNA,
        switch (gOn[i]) { // determine integer representation of on-target base
          case 'A': onBase = 0;
              break;
          case 'T': onBase = 1;
              break;
          case 'G': onBase = 2;
              break;
          case 'C': onBase = 3;
              break;
        }
        switch (gOff[i]) { // determine integer representation of off-target base. REVERSE-COMPLEMENTED!
          case 'T': offBase = 0;
              break;
          case 'A': offBase = 1;
              break;
          case 'C': offBase = 2;
              break;
          case 'G': offBase = 3;
              break;
        }
        score = score * cfdMatrix[16*i + 4*onBase + offBase]; // multiply score by the weight of this particular mismatch (or match) at this specific position. cfdMatrix is 320 (20*4*4) weights long, and is sorted by position first, on-target base identity second and reverse-complemented off-target base identity last.
		}

    score = score * modifier * 100.0; // scale score to 0-100. Note this is different than the original 0-1 scale used by Doench et al. (2016).

		return score;
}
