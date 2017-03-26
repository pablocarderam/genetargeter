from py.BioUtils import *;
from os import walk;

#TODO: make function that processes string input of message files separated by :::
#TODO: make function that processes list of message file strings for each message file by splitting by lines,
# searching for each design step's key words and then filling in success, error, warning or not achieved,
# and returns string output of csv file containing results in tabular form. Link to designs?

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
    outStr = "Code,Cas9,Cpf1";
    for c in rDict1:
        outStr = outStr + c + ", " + rDict1[c] + ", " + rDict2[c] + "\n";

    output(outStr, "outResultsTable.csv");


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



    output(outStr,'outResultsTable.csv')

'''
for file in *.txt; do echo ':::' >> "$file"; done
cat Mess* >> msgs.txt
###
from py.ProcessOutput import *
txt = open('msgs.txt');
d = txt.read(); # Read file
txt.close()
d=processMsgFileStrings(d)
o=createOutputSummaryTable(d)
'''

def processMsgFileStrings(fileString):
    files = fileString.split(":::");
    return files;

def createOutputSummaryTable(msgFiles):
    geneResults = {};
    for msg in msgFiles:
        if msg.find("PF3D7") > -1:
            geneName = msg[msg.find("PF3D7"):msg.find("PF3D7")+13];
            results = [""]*8;
            if geneName in geneResults:
                results = geneResults[geneName];

            enzymeMod = -1;
            if msg.find("Cas9") > -1:
                enzymeMod = 0;
            elif msg.find("Cpf1") > -1:
                enzymeMod = 4;

            for i in range(4):
                results[i+enzymeMod] = "Success";

            lines = msg.split("\n");
            for l in lines:
                if l.find("Warning") > -1:
                    if l.find("LHR") > -1:
                        results[enzymeMod] = "Warning";
                    elif l.find("gBlock") > -1:
                        results[enzymeMod+1] = "Warning";
                    elif l.find("RHR") > -1:
                        results[enzymeMod+2] = "Warning";
                    elif l.find("gRNA") > -1:
                        results[enzymeMod+3] = "Warning";


                if l.find("ERROR") > -1:
                    if l.find("LHR") > -1:
                        results[enzymeMod] = "Error";
                    elif l.find("gBlock") > -1:
                        results[enzymeMod+1] = "Error";
                    elif l.find("RHR") > -1:
                        results[enzymeMod+2] = "Error";
                    elif l.find("gRNA") > -1:
                        results[enzymeMod+3] = "Error";


            geneResults[geneName] = results;

    tabStr = "Gene_ID,Cas9_LHR,Cas9_gBlock,Cas9_RHR,Cas9_gRNA,Cpf1_LHR,Cpf1_gBlock,Cpf1_RHR,Cpf1_gRNA,Enzyme_recommended\n";
    for g in geneResults:
        bestOption = "Cas9";
        if geneResults[g][0:4].count("Success") < geneResults[g][4:8].count("Success"):
            bestOption = "Cpf1";

        geneResults[g].append(bestOption);
        tabStr = tabStr + g + "," + ",".join(geneResults[g]) + "\n";

    return tabStr;
