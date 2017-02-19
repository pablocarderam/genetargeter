from py.BioUtils import *;

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
