from py.utils.BioUtils import *;
from py.utils.GenBankToolbox import *;
from os import walk;
import os;

#TODO: make function that processes string input of message files separated by ::: EDIT not anymore? just CLI?
#TODO: make function that processes list of message file strings for each message file by splitting by lines,
# searching for each design step's key words and then filling in success, error, warning or not achieved,
# and returns string output of csv file containing results in tabular form. DONE, I think. Link to designs?

'''
USE JUST THIS: put all files for both enzymes into output/ directory, then:
python
from py.auxiliary.ProcessOutput import *
processOutput()
'''

def processOutput(dirPath='output/'):
    o=createOutputSummaryTable(dirPath);
    o=createOutputSummaryTableNonProtein(dirPath);
    #gbs=generateConcatGenBanks(dirPath);
    #createOutputDirStructure(dirPath);
    print 'Done postprocessing'


def createOutputSummaryTable(dirPath):
    geneResults = {};
    for (dirpath, dirnames, filenames) in walk(dirPath):
        for f in filenames:
            if f.find("Message_File") > -1:
                geneName = f[13:f.find('_pSN')];
                msg = "\n".join(loadSeqs(dirPath+f));
                results = [""]*12;
                if geneName in geneResults:
                    results = geneResults[geneName];

                enzymeMod = -1;
                if msg.find("Cas9") > -1:
                    enzymeMod = 0;
                elif msg.find("Cas12") > -1:
                    enzymeMod = 6;


                status = "Success";
                if msg.find("ERROR") > 0:
                    status = "Error";
                elif msg.find("Warning") > 0:
                    status = "Warning";

                results[enzymeMod] = status;
                for i in range(4):
                    results[i+enzymeMod+1] = "Success";
                    results[enzymeMod+5] = "0";

                lines = msg.split("\n");
                for l in lines:
                    if l.find("Warning") > -1:
                        if l.find("LHR") > -1:
                            results[enzymeMod+1] = "Warning";
                        elif l.find("gBlock") > -1 or l.find("Recoded") > -1:
                            results[enzymeMod+2] = "Warning";
                        elif l.find("RHR") > -1:
                            results[enzymeMod+3] = "Warning";
                        elif l.find("gRNA") > -1:
                            results[enzymeMod+4] = "Warning";


                    if l.find("ERROR") > -1:
                        if l.find("LHR") > -1:
                            results[enzymeMod+1] = "Error";
                        elif l.find("gBlock") > -1 or l.find("Recoded") > -1:
                            results[enzymeMod+2] = "Error";
                        elif l.find("RHR") > -1:
                            results[enzymeMod+3] = "Error";
                        elif l.find("gRNA") > -1:
                            results[enzymeMod+4] = "Error";
                            results[enzymeMod+1] = "-";
                            results[enzymeMod+2] = "-";
                            results[enzymeMod+3] = "-";
                            results[enzymeMod+5] = "0";
                        elif l.find("No gene annotations") > -1:
                            results[enzymeMod+1] = "Error";
                            results[enzymeMod+2] = "Error";
                            results[enzymeMod+3] = "Error";
                            results[enzymeMod+4] = "Error";
                            results[enzymeMod+5] = "0";
                        elif l.find("protein-coding sequence") > -1:
                            results[enzymeMod+1] = "NA";
                            results[enzymeMod+2] = "NA";
                            results[enzymeMod+3] = "NA";
                            results[enzymeMod+4] = "NA";
                            results[enzymeMod+5] = "NA";


                    if l.find("Recoded region with size ") > -1:
                        results[enzymeMod+5] = l[l.find("Recoded region with size ")+len("Recoded region with size "):l.find(" for gene ")];
                    elif l.find("Recoded region not deemed necessary") > -1:
                        results[enzymeMod+5] = "0";


                geneResults[geneName] = results;



    tabStr = "Gene_ID,Cas9_Status,Cas9_LHR,Cas9_gBlock,Cas9_RHR,Cas9_gRNA,Cas9_gBlock_size,Cas12_Status,Cas12_LHR,Cas12_gBlock,Cas12_RHR,Cas12_gRNA,Cas12_gBlock_size,Enzyme_recommended\n";
    for g in geneResults:
        bestOption = "Cas9";
        if geneResults[g][0:5].count("Error") > 0 and geneResults[g][6:11].count("Error") > 0:
            bestOption = "None";
        elif geneResults[g][0:5].count("Success") < geneResults[g][6:11].count("Success") and geneResults[g][6:11].count("Error") == 0:
            bestOption = "Cas12";
        elif geneResults[g][0:5].count("Success") == geneResults[g][6:11].count("Success"):
            if int(geneResults[g][5]) < int(geneResults[g][11]):
                bestOption = "Cas9";
            elif int(geneResults[g][5]) > int(geneResults[g][11]):
                bestOption = "Cas12";
            else:
                bestOption = "Either";

        geneResults[g].append(bestOption);
        tabStr = tabStr + g + "," + ",".join(geneResults[g]) + "\n";

    output(tabStr,"resultsTable.csv", wipe=True);
    return tabStr;


'''
Same method as one above, except outputs only non protein-coding genes
'''
def createOutputSummaryTableNonProtein(dirPath):
    geneResults = {};
    for (dirpath, dirnames, filenames) in walk(dirPath):
        for f in filenames:
            if f.find("Message_File") > -1:
                geneName = f[13:28];
                msg = "\n".join(loadSeqs(dirPath+f));
                if msg.find("not to be a protein-coding") > -1:
                    enzymeMod = 0;

                    results = [""]*12;
                    if geneName in geneResults:
                        results = geneResults[geneName];

                    enzymeMod = -1;
                    if msg.find("Cas9") > -1:
                        enzymeMod = 0;
                    elif msg.find("Cas12") > -1:
                        enzymeMod = 6;


                    status = "Success";
                    if msg.find("ERROR") > 0:
                        status = "Error";
                    elif msg.find("Warning") > 0:
                        status = "Warning";

                    results[enzymeMod] = status;
                    for i in range(4):
                        results[i+enzymeMod+1] = "Success";
                        results[enzymeMod+5] = "0";

                    lines = msg.split("\n");
                    for l in lines:
                        if l.find("Warning") > -1:
                            if l.find("LHR") > -1:
                                results[enzymeMod+1] = "Warning";
                            elif l.find("gBlock") > -1 or l.find("Recoded") > -1 or l.find("gene fragment") > -1:
                                results[enzymeMod+2] = "Warning";
                            elif l.find("RHR") > -1:
                                results[enzymeMod+3] = "Warning";
                            elif l.find("gRNA") > -1:
                                results[enzymeMod+4] = "Warning";


                        if l.find("ERROR") > -1:
                            if l.find("LHR") > -1:
                                results[enzymeMod+1] = "Error";
                            elif l.find("gBlock") > -1 or l.find("Recoded") > -1 or l.find("gene fragment") > -1:
                                results[enzymeMod+2] = "Error";
                            elif l.find("RHR") > -1:
                                results[enzymeMod+3] = "Error";
                            elif l.find("gRNA") > -1:
                                results[enzymeMod+4] = "Error";
                                results[enzymeMod+1] = "-";
                                results[enzymeMod+2] = "-";
                                results[enzymeMod+3] = "-";
                                results[enzymeMod+5] = "0";
                            elif l.find("No gene annotations") > -1:
                                results[enzymeMod+1] = "Error";
                                results[enzymeMod+2] = "Error";
                                results[enzymeMod+3] = "Error";
                                results[enzymeMod+4] = "Error";
                                results[enzymeMod+5] = "0";
                            elif l.find("protein-coding sequence") > -1:
                                results[enzymeMod+1] = "NA";
                                results[enzymeMod+2] = "NA";
                                results[enzymeMod+3] = "NA";
                                results[enzymeMod+4] = "NA";
                                results[enzymeMod+5] = "NA";


                        if l.find("Recoded region with size ") > -1:
                            results[enzymeMod+5] = l[l.find("Recoded region with size ")+len("Recoded region with size "):l.find(" for gene ")];
                        elif l.find("Recoded region not deemed necessary") > -1:
                            results[enzymeMod+5] = "0";


                    geneResults[geneName] = results;



    tabStr = "Gene_ID,Cas9_Status,Cas9_LHR,Cas9_gBlock,Cas9_RHR,Cas9_gRNA,Cas9_gBlock_size,Cas12_Status,Cas12_LHR,Cas12_gBlock,Cas12_RHR,Cas12_gRNA,Cas12_gBlock_size,Enzyme_recommended\n";
    for g in geneResults:
        bestOption = "Cas9";
        if geneResults[g][0:5].count("Error") > 0 and geneResults[g][6:11].count("Error") > 0:
            bestOption = "None";
        elif geneResults[g][0:5].count("Success") < geneResults[g][6:11].count("Success") and geneResults[g][6:11].count("Error") == 0:
            bestOption = "Cas12";
        elif geneResults[g][0:5].count("Success") == geneResults[g][6:11].count("Success"):
            if int(geneResults[g][5]) < int(geneResults[g][11]):
                bestOption = "Cas9";
            elif int(geneResults[g][5]) > int(geneResults[g][11]):
                bestOption = "Cas12";
            else:
                bestOption = "Either";

        geneResults[g].append(bestOption);
        tabStr = tabStr + g + "," + ",".join(geneResults[g]) + "\n";

    output(tabStr,"resultsTableNonProtein.csv", wipe=True);
    return tabStr;


'''
Used to create tree structure containing one folder per gene, each with one
folder per enzyme, each with the six output files.

RUN AFTER GENERATING OUTPUTS/COPYING MSG FILES SOMEWHERE ELSE

Put all output files in a single directory, then run this:

from py.ProcessOutput import *
createOutputDirStructure('output/nameOfDir/')

Will print an error at end if trying to process DS_Store, but should be fine
'''
def createOutputDirStructure(fileDir):
    if not os.path.exists(fileDir+"/output_by_gene/"):
        os.makedirs(fileDir+"/output_by_gene/");

    for (dirpath, dirnames, filenames) in walk(fileDir):
        for f in filenames:
            fileInfo = f.split(".txt")[0].split(".csv")[0].split(".gb")[0].split("_");
            enzyme = fileInfo[-1];
            geneName = "";
            if fileInfo[0] == "Oligos":
                geneName = "_".join(fileInfo[1:-1])
            else:
                geneName = "_".join(fileInfo[2:-1])

            if not os.path.exists(fileDir+"/output_by_gene/"+geneName):
                os.makedirs(fileDir+"/output_by_gene/"+geneName);

            if not os.path.exists(fileDir+"/output_by_gene/"+geneName+"/"+enzyme):
                os.makedirs(fileDir+"/output_by_gene/"+geneName+"/"+enzyme);

            os.rename(fileDir+"/"+f,fileDir+"/output_by_gene/"+geneName+"/"+enzyme+"/"+f)



def generateConcatGenBanks(fileDir,gbToGenerate="ALL"):
    if not os.path.exists(fileDir+"/gb_concats/"):
        os.makedirs(fileDir+"/gb_concats/");

    chromosomes = range(1,15)+["MIT","API"];
    for chromosome in chromosomes:
        if not os.path.exists(fileDir+"/gb_concats/"+"Chromosome_"+str(chromosome)):
            os.makedirs(fileDir+"/gb_concats/"+"Chromosome_"+str(chromosome));

    gbs = {"Locus_Pre-editing_Cas9":GenBank(),"Locus_Post-editing_Cas9":GenBank(),"pSN054_V5_Cas9":GenBank(),
           "Locus_Pre-editing_Cas12":GenBank(),"Locus_Post-editing_Cas12":GenBank(),"pSN054_V5_Cas12":GenBank()};
    genome = dict( (chrom,deepcopy(gbs)) for chrom in chromosomes );

    lenSep = int(1e3);
    sep = "N"*lenSep;

    for (dirpath, dirnames, filenames) in walk(fileDir):
        pctSize = len(filenames)/100;
        fileCt = 0;
        pct = 0;

        for f in filenames:
            if f.find(".gb") > -1:
                gb = GenBank();
                gb.load(os.path.join(dirpath,f),loadFromFile=True);
                if gb.origin > 0:
                    correctedFileName = f.replace("Pre-editing _","Pre-editing_").replace("(1).gb","").replace(".gb","").strip();
                    splitName = correctedFileName.split("PF3D7_")
                    chromosome = "MIT";
                    if len(splitName) > 1:
                        chromosome = splitName[1][:2];
                        if chromosome == "AP":
                            chromosome = "API";
                        else:
                            chromosome = int(chromosome);



                    gbToAdd = '_'.join( correctedFileName.split('_')[:2] + [ correctedFileName[-4:] ] )
                    if chromosome == gbToGenerate or gbToGenerate == "ALL":
                        lenSeq = len(genome[chromosome][gbToAdd].origin);

                        genome[chromosome][gbToAdd].insertSeq(sep,lenSeq);
                        genome[chromosome][gbToAdd].insertGenBank(gb,lenSeq+lenSep);



            fileCt += 1;
            if fileCt > pctSize:
                pct += 1;
                fileCt = 0;
                print str(pct) + '%';

        for chromosome in chromosomes:
            for concat in genome[chromosome].keys():
                if len(genome[chromosome][concat].origin) > 0:
                    genome[chromosome][concat].removeSeq([0,lenSep]);
                    print "Removed initial separator";
                    fName = "Chromosome_"+str(chromosome)+'_'+concat;
                    genome[chromosome][concat].name = fName;
                    genome[chromosome][concat].save(fileDir+"/gb_concats/Chromosome_"+str(chromosome)+"/"+fName+".gb",saveToFile=True);
                    print "Saved " + fName;

        return gbs;







###
### NOT NEEDED ANY MORE:
###
def resultsDict(filePath):
    lines = loadSeqs(filePath);
    results = {};
    for l in lines:
        codeStart = l.find(". ");
        code = l[codeStart+2:codeStart+15];
        resultStart = l.find("|| ");
        result = l[resultStart+3:resultStart+4];
        results[code] = result;

    return results;

def outputResultsTable(rDict1,rDict2):
    outStr = "Code,Cas9,Cas12";
    for c in rDict1:
        outStr = outStr + c + ", " + rDict1[c] + ", " + rDict2[c] + "\n";

    output(outStr, "outResultsTable.csv", wipe=True);


def processErrorFiles(dirPath):
    outStr = "Code,Status\n";
    for (dirpath, dirnames, filenames) in walk(dirPath):
        for f in filenames:
            outStr = outStr + f[13:26];
            ftxt = "\n".join(loadSeqs(dirPath+f));
            if ftxt.find("ERROR") > 0:
                outStr = outStr + ',' + 'E' + '\n';
            elif ftxt.find("Warning") > 0:
                outStr = outStr + ',' + 'W' + '\n';
            else:
                outStr = outStr + ',' + 'S' + '\n';



    output(outStr,'outResultsTable.csv', wipe=True)

'''
NOT NEEDED ANY MORE:
for file in *.txt; do echo ':::' >> "$file"; done
cat Mess* >> msgs.txt
###
for f in *.txt; do mv -- "$f" "${f%.txt}_Cas9.txt"; done
###
'''

def processMsgFileStrings(fileString):
    files = fileString.split(":::");
    return files;
