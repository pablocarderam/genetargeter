# -*- coding: utf-8 -*-
"""
Created on Sat Jun 5 2017
Python functions used to get gRNA info according to Hsu et al. (2013) and Doench
et al. (2016) algorithms from a pre-calculated database.
@author: Pablo CR
"""

def searchDB(searchSeq,dbFilePath):
    returnInfo = [];
    with open(dbFilePath) as fp:
        for i,row in enumerate(fp):
            r = row.split(','); # split row into elements
            #print r[0].upper().replace('"',"")
            if r[0].upper().replace('"',"") == searchSeq: # if gRNA found,
                r = row.replace("NA","0").replace('"',"").split(',');
                r[0] = r[0].upper(); # extended seq to uppercase
                r[1] = float(r[1]); # Azimuth as float,
                r[2] = float(r[2]); # CINDEL as float,
                r[3] = float(r[3]); # Hsu_Normal as float,
                r[4] = float(r[4]); # Hsu_Normal_Max as float,
                r[5] = int(r[5]); # Hsu_Normal_NumHits as int,
                r[6] = float(r[6]); # Hsu_Agnostic as float,
                r[7] = float(r[7]); # Hsu_Agnostic_Max as float,
                r[8] = int(r[8]); # Hsu_Agnostic_NumHits as int,
                r[9] = float(r[9]); # CFD as float,
                r[10] = float(r[10]); # CFD_Max as float,
                r[11] = int(r[11]); # CFD_NumHits as float,
                returnInfo = r;
                break;



    return returnInfo;

def getGRNAInfoFromDB(pExtSeq,enzyme):
    searchSeq = pExtSeq.upper(); # to uppercase
    returnInfo = [];
    dbFilePaths = ['./gRNAScores/data/cas9Small.csv','./gRNAScores/data/cas12Small.csv']
    fp='';
    if enzyme=='Cas9':
        fp='./gRNAScores/data/cas9Small.csv';
    elif enzyme == 'Cas12':
        fp='./gRNAScores/data/cas12Small.csv';

    returnInfo = searchDB(searchSeq,fp);
    return returnInfo;



'''grnaDBFile1 = open('./gRNAScores/grnasSmallDB_1.csv'); # Access grna DB file
grnaDBFile2 = open('./gRNAScores/grnasSmallDB_2.csv'); # Access grna DB file
d = grnaDBFile1.read() + grnaDBFile2.read(); # Read file
d = d.replace("NA","0");
d = d.replace('"',"");
d = d.split("\n"); # split into lines
grnaDB = [r.split(",") for r in d]; # split each line into elements
grnaDB = grnaDB[2:len(grnaDB)-1]; # remove first line (headers)
if grnaDB[len(grnaDB)-1] < 13: # if last line not complete (empty),
    grnaDB = grnaDB[1:len(grnaDB)-2]; # remove last line

newGRNADB =[]; # will store new list
for r in grnaDB: # loop over all gRNAs,
    if len(r) == 12 and r[0].find("Extended_Seq") < 0: # if last line not complete (empty),
        r[0] = r[0].upper(); # extended seq to uppercase
        r[1] = float(r[1]); # Azimuth as float,
        r[2] = float(r[2]); # CINDEL as float,
        r[3] = float(r[3]); # Hsu_Normal as float,
        r[4] = float(r[4]); # Hsu_Normal_Max as float,
        r[5] = int(r[5]); # Hsu_Normal_NumHits as int,
        r[6] = float(r[6]); # Hsu_Agnostic as float,
        r[7] = float(r[7]); # Hsu_Agnostic_Max as float,
        r[8] = int(r[8]); # Hsu_Agnostic_NumHits as int,
        r[9] = float(r[9]); # CFD as float,
        r[10] = float(r[10]); # CFD_Max as float,
        r[11] = int(r[11]); # CFD_NumHits as float,
        newGRNADB.append(r); # add to new array

grnaDB = sorted(newGRNADB, key=lambda r: r[0]); # sort rows by their extended sequences

def getGRNAInfoFromDB(pExtSeq):
    searchSeq = pExtSeq.upper(); # to uppercase
    grnaRow = []; # will store this gRNA's information
    upSearchLim = len(grnaDB)-1; # upper search limit for binary search
    loSearchLim = 0; # lower search limit for binary search
    midPoint = int((upSearchLim-loSearchLim)/2)+loSearchLim; # midpoint
    while len(grnaRow) == 0 and midPoint != loSearchLim: # while grna not found and still room for finding it,
        if searchSeq > grnaDB[midPoint][0]: # if extended seq greater than midpoint row's extended seq,
            loSearchLim = midPoint; # reset lower search limit for binary search
            midPoint = int((upSearchLim-loSearchLim)/2)+loSearchLim; # reset midpoint
        elif searchSeq < grnaDB[midPoint][0]: # if extended seq less than midpoint row's extended seq,
            upSearchLim = midPoint; # reset lower search limit for binary search
            midPoint = int((upSearchLim-loSearchLim)/2)+loSearchLim; # reset midpoint
        else: # if extended sequence found,
            grnaRow = grnaDB[midPoint]; # set out variable

    return grnaRow;'''
