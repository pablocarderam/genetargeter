#!/usr/bin/env python

"""
Calls main method to build targeted vector with arguments from the command line.

# Usage:
# To run GeneTargeter as a command line application, use the following command:
#
# genetargeter PATH/TO/FILE.gb ParameterFile.txt PATH/TO/OUTPUTFOLDER NAME_OF_GENE
#
# The last parameter is optional; if empty, the name will be taken from the name
# of the file.
# Additionally, instead of the name of a single file, a directory may be listed.
# In this case, all files in the folder will be processed and each gene's name
# will be taken from the name of its file.
"""
from __future__ import print_function

from builtins import str
from builtins import range
import sys; # Needed for receiving user input
import os; # for joining file paths
import multiprocessing # for parallelizing genes
import joblib as jl # for parallelizing genes
from py.genetargeter.constants import *; # main python library in py folder
from py.genetargeter.GeneTargeterMethods import *; # main python library in py folder
from py.genetargeter.inputProcessing import *; # input handling

'''
Gets parameters from file passed in function call
'''
def getParams(pFile):
    txt = open(pFile); # Access given file
    d = txt.read(); # Read file
    txt.close(); # Close file

    d = d.replace(':','=') # replace to get Python syntax, remove all whitespace
    d = d.split("\n"); # Split into separate lines
    d = [ l[0:l.find('#')] for l in d ] # remove comments
    while("" in d):
        d.remove("")  # remove empty lines

    return d

"""
Takes file
Returns a string with file contents
"""
def getGeneFileStr(pFile):
    if pFile[-3:] == '.gb':
        txt = open(pFile); # Access given file
        d = txt.read(); # Read file
        txt.close(); # close file
        return d; # returns dictionary
    else:
        print('ERROR in '+pFile+' : only GenBank files processed.')

def runGene(geneName, geneFileStr, params, outputDir, geneNum=0):
    params_dict = {}
    for l in params: # for every line in parameter file,
        if l.find('=') > -1:
            params_dict[l.split()[0]] = eval(l.split('=')[-1].strip())

    if params_dict['prefix'] == "*None*":
        params_dict['prefix'] = "";
    elif params_dict['prefixNum'] != "" and params_dict['prefixNum'] > -1:
        params_dict['prefix'] += str(params_dict['prefixNum']+geneNum).zfill(3)+'_';
    else:
        params_dict['prefix'] += '_';

    params_dict['minGRNAGCContent']=params_dict['minGRNAGCContent']/100.0 # convert to decimal
    outMsg = ""; # will store output message
    outMsgHA = "HA_Tag_Design-" # will store output message for HA tag design, if any

    geneGBs = preprocessInputFile(geneName, geneFileStr, useFileStrs=True, doingAltSplices=True); # obtain gb file(s) to be processed

    for gbName in geneGBs: # for every gb
        sigPep = chkSignalPeptide5Prime(gbName,signalPDB); # signalP info
        params_dict['haTag'] = False; # don't use HA tags by default
        if params_dict['plasmidType'] == "pSN150" and ( params_dict['haTag'] == True or (params_dict['haTag'] == "Auto" and not sigPep) ): # if forcing HA tags or 5' end does not contain signal peptide and auto HA-tagging
            params_dict['haTag'] = True; # use HA tags

        if params_dict['filterCutSites'] == "Auto": # if filterCutSites is automatic
            params_dict['filterCutSites'] = cut_sites[params_dict['plasmidType']] # set to its plasmid

        output = targetGene(gbName, geneGBs[gbName], codonOptimize=params_dict['codonOptimize'], useFileStrs=False, outputDir=outputDir, HRannotated=params_dict['HRannotated'],lengthLHR=params_dict['lengthLHR'], lengthRHR=params_dict['lengthRHR'], gibsonHomRange=params_dict['gibsonHomRange'], optimRangeLHR=params_dict['optimRangeLHR'], optimRangeRHR=params_dict['optimRangeRHR'], endSizeLHR=params_dict['endSizeLHR'], endSizeRHR=params_dict['endSizeRHR'], endTempLHR=params_dict['endTempLHR'], endTempRHR=params_dict['endTempRHR'], gibTemp=params_dict['gibTemp'], gibTDif=params_dict['gibTDif'], maxDistLHR=params_dict['maxDistLHR'], maxDistRHR=params_dict['maxDistRHR'], minGBlockSize=params_dict['minGBlockSize'], codonSampling=params_dict['codonSampling'], minGRNAGCContent=params_dict['minGRNAGCContent'], onTargetMethod=params_dict['onTargetMethod'], minOnTargetScore=params_dict['minOnTargetScore'], offTargetMethod=params_dict['offTargetMethod'], minOffTargetScore=params_dict['minOffTargetScore'], maxOffTargetHitScore=params_dict['maxOffTargetHitScore'], enzyme=params_dict['enzyme'], PAM=params_dict['PAM'], gBlockDefault=params_dict['gBlockDefault'], plasmidType=params_dict['plasmidType'], haTag=params_dict['haTag'], sigPep=sigPep, setCoding=params_dict['setCoding'], bulkFile=params_dict['bulkFile'], prefix=params_dict['prefix'], filterCutSites=params_dict['filterCutSites']); # call result


def parallelRun(file,params,outputDir,geneNum):
    if file[-3:] == '.gb':
        geneFileStr = getGeneFileStr(file);
        geneName = file[file.rfind('/')+1:file.find('.')] # get gene name from file
        runGene(geneName, geneFileStr, params, outputDir, geneNum); # call result


def runAll(args=None):
    if args is None:
        args = sys.argv

    if len(args) == 5:
        script, geneFile, paramsFile, outputDir, geneName = args[0:5]; # From console, args returns script name, arguments
    elif len(args) == 4:
        script, geneFile, paramsFile, outputDir = args[0:4]; # From console, args returns script name, arguments
        geneName = geneFile[geneFile.rfind('/')+1:geneFile.find('.')] # get gene name from file
    else:
        print("Wrong number of arguments! Pass GENEBANK_FILE.gb ParameterFile.txt outputDirectory NAME_OF_GENE (last one is optional)")
        exit()


    params = getParams(paramsFile);

    if os.path.isdir(geneFile): # if path given is a directory,
        files = os.listdir(geneFile) # get full list
        c = 0; # counts files to process
        geneNums = [] # will contain indexes
        geneFiles = [] # will contain files
        for file in files: # add file number to each file
            if file[-3:] == '.gb': # if a genbank file,
                geneNums.append(c); # add to number list
                geneFiles.append(file); # add to file list


        n_jobs = min(len(files),multiprocessing.cpu_count()) # jobs is min of num cores and files
        jl.Parallel(n_jobs=n_jobs,verbose=10) (jl.delayed(parallelRun)(os.path.join(geneFile,geneFiles[i]),params,outputDir,geneNums[i]) for i in range(len(geneNums)))

    else:
        geneFileStr = getGeneFileStr(geneFile);
        runGene(geneName, geneFileStr, params, outputDir); # call result

    print(geneName + " targeted. Results in output folder.")


if __name__ == '__main__':
    runAll() # run all
