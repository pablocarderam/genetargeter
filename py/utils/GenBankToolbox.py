# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 19:24:17 2016
Classes useful in rosalind bioinformatic algorithms
@author: Pablo the awesome molecular jedi
"""
from __future__ import print_function
from __future__ import division

from builtins import str
from builtins import range
from builtins import object
from py.utils.BioUtils import *; # Imports utils
from copy import deepcopy; # Import object copying methods for deep copies, used in creating revComp of GenBank objects

class GenBankAnn(object):
    """Stores some GenBank format information for an annotation."""
    def __init__(self, pLab="", pType="misc", pSeq="", pComp=False, pIndex=[], pCol=""): # class constructor
        self.label = pLab # annotation label
        self.type = pType; # type of annotation
        self.seq = pSeq; # sequence of annotation
        self.comp = pComp; # false if on positive strand, true if on complementary strand
        self.index = pIndex; # list with region start and end indexes IN PYTHON INDEXING (starting on 0)
        if pCol == "": # if no color specified,
            pCol = randHex(); # set a random coloring

        self.color = pCol # annotation color


class GenBank(object):
    """
    Stores some GenBank format information for a file.
    Attributes:
    name = name of sequence
    additionalInfo = string containing molecule information ("ds-DNA"), topology information (ie "circular", "linear"), GenBank division, whatever else
    date = date of creation
    definition: a string containing the first line of the file
    features: a list with objects with the following attributes:
        label = annotation labels
        type = type of annotation
        seq = sequence of annotation on the strand it actually is in
        comp = false if on positive strand, true if on complementary strand
        index = list with region start and end indexes IN PYTHON INDEXING (starting on 0)

    origin: a string containing the full sequence.

    Methods:
    GenBank() = Constructor
    reset() = resets attributes
    load(pFile) = Loads a GenBank class object from GenBank file.
    save(pFile) = Saves a GenBank class object to a GenBank file.
    insertSeq(seq, index) = inserts sequence into GenBank object, modifies annotation indexes accordingly
    removeSeq(indexes) = removes sequence from GenBank object, modifies annotation indexes accordingly
    insertGenBank(gb, index) = inserts another GenBank object into GenBank object, modifies annotation indexes accordingly
    """

    def __init__(self): # class constructor
        self.reset(); # resets attributes

    def reset(self): # resets attributes
        self.name = ""; # initializes name string
        self.additionalInfo = ""; # initializes additionalInfo string
        self.date = ""; # initializes date string
        self.definition = ""; # initializes definition string
        self.origin = ""; # initializes full sequence string
        self.features = []; # initializes a new list


    """
    Loads a GenBank-class object from raw file contents passed as a string in
    pFile or from a .gb file saved in the filepath pFile (set loadFromFile =
    True for this option).
    """
    def load(self, pFile, loadFromFile=False):
        self.reset(); # resets attributes
        if loadFromFile: # if a file path has been passed in pFile,
            txt = open(pFile); # Access given file
            d = txt.read(); # Read file.
        else: # if a file string has been passed in pFile,
            d = pFile; # save the read file string as raw data in d

        if len(d) == 0: # if data is empty,
            print("ERROR: empty file or filestring passed for GenBank loading.")
        else: # if there's data,
            d = d.split("\n"); # Split into lines. d is now a list of strings corresponding to the lines of the file.

            i = 0; # indexes through lines
            w = d[i].split(); # splits line into words

            while len(w) == 0 or w[0] != "LOCUS" and i < len(d): # iterates until LOCUS found or end of document reached
                i = i+1; # advance counter
                w = d[i].split(); # splits line into words

            if w[0] == "LOCUS": # if locus found...
                self.name = w[1]; # saves name
                self.additionalInfo = " ".join(w[4:(len(w)-1)]); # saves additional information
                self.date = w[len(w)-1]; # date is always last in line

            while w[0] != "DEFINITION" and i < len(d): # iterates until DEFINITION found or end of document reached
                i = i+1; # advance counter
                w = d[i].split(); # splits line into words

            if w[0] == "DEFINITION": # if definition found...
                self.definition = d[i].lstrip()[12::]; # saves definition information on this line after "DEFINITION" string and ensuing whitespace.
                while w[len(w)-1].find(".") < 0 and i < len(d): # iterates until last word in line contains period, or end of document reached
                    i = i+1; # advance counter
                    w = d[i].split(); # splits line into words
                    self.definition = self.definition + d[i]; # save this line as definition information

            # We need to load the sequence before we load the annotations
            while w[0] != "ORIGIN" and i < len(d): # iterates until ORIGIN found or end of document reached
                i = i+1; # advance counter
                w = d[i].split(); # splits line into words

            if w[0] == "ORIGIN": # if origin found...
                i = i+1; # advance counter
                w = d[i].split(); # splits line into words
                while w[0] != "//" or i == len(d): # iterates until end of document reached
                    self.origin = self.origin + "".join(w[1:]); # joins all sequence information in this line into the full sequence. The first word is omitted due to the fact it is supposed to be a position number.
                    i = i+1; # advance counter
                    w = d[i].split(); # splits line into words

            # Now the annotations
            i = 0; # reset the line counter to search for the features section
            while w[0] != "FEATURES" and i < len(d): # iterates until FEATURES found and end of document reached
                i = i+1; # advance counter
                w = d[i].split(); # splits line into words

            if w[0] == "FEATURES": # if features found...
                newAnn = GenBankAnn("","","",False,[]); # creates new empty GenBankAnn class object
                i = i+1; # advance counter
                w = d[i].split(); # splits line into words
                while w[0] != "ORIGIN" and i < len(d): # iterates until ORIGIN found or end of document reached, assuming ORIGIN comes after features...
                    if d[i].lstrip()[0] != "/" and len(w) == 2: # if the line does not start with "/" and contains two strings separated by spaces, it means it might be a new annotation. lstrip() removes left whitespace.
                        if len(newAnn.label)*len(newAnn.type)*len(newAnn.seq)*len(newAnn.index) != 0: # if all properties in newAnn have been filled in,
                            self.features.append(newAnn); # adds the new annotation
                            newAnn = GenBankAnn("","","",False,[]); # resets newAnn variable

                        posRawString = w[1].split("("); # splits second word in w by "(" to see if on complementary strand
                        posString = posRawString[0]; # assigns posString to the indexes string in the case that the annotation is on the positive strand
                        if len(posRawString) > 1: # if split into more than one piece, it means it's on the complementary strand
                            posString = posRawString[1][:(len(posRawString[1])-1)]; # extracts index information from raw string

                        if posString.split("..")[0].isdigit(): # if index information is actually a number (avoids cases of continuing description lines)
                            newAnn = GenBankAnn("","","",False,[]); # resets newAnn variable
                            newAnn.type = w[0]; # stores annotation type

                            if len(posRawString) > 1: # if split into more than one piece, it means it's on the complementary strand
                                newAnn.comp = True; # set comp to true

                            newAnn.index = [int(index) for index in posString.split("..")]; # takes the index string, splits it into two strings separated by ".." and converts them both to ints.
                            if len(newAnn.index) == 1: # if only one index found (just one bp annotated),
                                newAnn.index.append(newAnn.index[0]); # duplicate the index

                            newAnn.index[0] = newAnn.index[0] - 1; # changes to Python indexing. Upper limit does not have to be changed since Python does not include it in indexing, while GenBank does.
                            newAnn.label = ''

                    elif d[i].find('=') > 0: # if line contains = character, it means it's a new property
                        p = d[i].lstrip().split("="); # splits line into a property's name and value
                        if p[0] == "/gene" or p[0] == "/label" or p[0] == "/ID": # if the property is a gene or label...
                            newAnn.label = (newAnn.label + p[1].strip('"').strip("'")).strip(); # assigns the label (extracting the label from string p)
                            if newAnn.comp: # if on complementary strand,
                                newAnn.seq = revComp(self.origin[(newAnn.index[0]):newAnn.index[1]]); # assigns the sequence of the annotated region to the complement. To find the sequence of the annotated region, indexes in the full sequence according to the position previously extracted.
                            else: # if on positive strand,
                                newAnn.seq = self.origin[(newAnn.index[0]):newAnn.index[1]]; # assigns the sequence of the annotated region. To find the sequence of the annotated region, indexes in the full sequence according to the position previously extracted.

                        elif p[0] == "/ApEinfo_fwdcolor" or p[0] == "/ApEinfo_revcolor": # if property is a color,
                            newAnn.color = p[1]; # save color


                    else: # if line does not contain = character, it's a continuation of the previous line.
                        pass; # do nothing since we're not loading descriptions anyway. Keep your labels short. Sorry.

                    i = i+1; # advance counter
                    w = d[i].split(); # splits line into words

                if len(newAnn.label)*len(newAnn.type)*len(newAnn.seq)*len(newAnn.index) != 0: # if all properties in newAnn have been filled in,
                    self.features.append(newAnn); # adds the new annotation
                    newAnn = GenBankAnn("","","",False,[]); # resets newAnn variable


    """
    Makes GenBank format file with annotations based on GenBank object. Returns
    a string containing the raw file string for a .gb file, or saves it to file
    path pFile (if saveToFile = True). Only encodes locus, definition, features,
    origin information.
    """
    def save(self, pFile, saveToFile=False):
        outStr = ""; # variable contains output string

        outStr = outStr + "LOCUS       " + self.name + " " + str(len(self.origin)) + " bp " + self.additionalInfo + "     " + todayDateStr() + "\n"; # writes LOCUS line

        d = "DEFINITION  " + self.definition; # Stores string with definition and start string
        defOut = d; #"\n".join([d[i:i+80] for i in range(0, len(d), 80)]); # introduces line separators into string ##WILL NOW PRINT STRING WITHOUT SEPARATORS

        outStr = outStr + defOut + "\nFEATURES             Location/Qualifiers\n"; # write locus and definition information to file, writes a new line to start features list

        for ann in self.features:
            annSeq = ann.seq; # default case: positive strand sequence
            if ann.comp: # if on complementary strand,
                annSeq = revComp(ann.seq); # save as rev comp

            if self.origin[ann.index[0]:ann.index[1]].upper() == annSeq.upper(): # if the annotation sequence was found in the region where it should be,
                startIndex = ann.index[0] + 1; # advances index (GenBank indexing starts at 1)
                endIndex = ann.index[1]; # stores end of sequence (GenBank starts at 1, but we don't shift this index because GenBank indexing includes final index)
                indexStr = str(startIndex) + ".." + str(endIndex); # save annotation indexes as string
                if ann.comp:
                    indexStr = "complement(" + indexStr + ")"; # if on complementary strand, rewrite indexStr with this information

                spacing = "".join([" " for i in range(16-len(ann.type))]); # creates enough spacing to complete 16 characters between tab and indexes
                outStr = outStr + "     " + ann.type + spacing + indexStr + "\n"; # write first line of annotation format
                outStr = outStr + '                     /label="' + ann.label + '"\n'; # write label line
                col = ann.color; # generates a random color for annotation
                outStr = outStr + '                     /ApEinfo_revcolor=' + col + '\n'; # write reverse color line
                outStr = outStr + '                     /ApEinfo_fwdcolor=' + col + '\n'; # write forward color line
            else: # if it wasn't found where it should be,
                print("ERROR: Annotation sequence " + ann.label + " not found in main sequence."); # warn of error.

        outStr = outStr + "ORIGIN\n"; # starts origin lines
        orSpaced = " ".join([self.origin[i:i+10] for i in range(0, len(self.origin), 10)]); # introduces spaces into string
        orOut = ""; # will store final sequence string to be saved in file
        for i in range(0, len(orSpaced), 66): # every 66 characters,
            indexStr = str(int(i - i/11 + 1)); # index in actual sequence (not counting space characters), starting from 1
            spacing = "".join([" " for j in range(9-len(indexStr))]); # creates enough spacing to complete 9 characters before sequence starts in this line
            orOut = orOut + spacing + indexStr + " " + orSpaced[i:i+66] + "\n"; # introduces line separators and index numbers into string

        orOut = orOut + "//\n"; # adds file ending
        outStr = outStr + orOut; # writes formatted sequence to file
        if saveToFile: # if user wants file saved,
            output(outStr,pFile,wipe=True); # write to file
        else: # if user wants file string,
            return outStr; # give them file string.


    """
    Inserts new sequence into GenBank object, modifies annotation indexes
    accordingly. The new sequence should be in 5' to 3', positive sense.
    Returns list of annotations whose sequence was altered by the elimination if
    it occurs accross one or more annotated regions, returns empty variable
    otherwise. The returned annotation contains the new sequence, not the
    original one, so it can be removed from the features list if necessary.
    """
    def insertSeq(self, seq, index):
        r = []; # initialize return variable
        self.origin = self.origin[:index] + seq + self.origin[index:]; # inserts into main sequence
        for ann in self.features: # loops through all annotations editing index numbers if necessary
            if len(ann.seq) > 0: # checks whether annotation actually contains information
                if ann.index[0] >= index: # if this annotation's start index is greater or equal to the insertion site,
                    ann.index[0] = ann.index[0] + len(seq); # shift the annotation's start index
                    ann.index[1] = ann.index[1] + len(seq); # shift the annotation's end index
                elif ann.index[1] > index: # if the annotation's start index is less than the insertion site, but the end index is strictly greater, the insertion happens inside the annotation site
                    ann.index[1] = ann.index[1] + len(seq); # shift the annotation's end index
                    ann.seq = self.origin[ann.index[0]:ann.index[1]]; # updates annotation sequence
                    r.append(ann); # appends annotation to return variable in case of annotation alteration
                # the final case is that the whole annotated region is upstream of the insertion site, where nothing changes.

        return r; # returns variable


    """
    Removes sequence between indexes given as list from GenBank object, modifies
    annotation indexes accordingly. The sequence removed is contained between
    indexes[0]:indexes[1], as per Python indexing. The indexes of the sequence
    to be removed should be in 5' to 3', positive sense. Returns list of
    annotations whose sequence was altered by the elimination if it occurs
    accross one or more annotated regions, returns empty variable otherwise. The
    returned annotation contains the new sequence, not the original one, so it
    can be removed from the features list if necessary.
    """
    def removeSeq(self, indexes):
        r = []; # initialize return variable
        lenRem = indexes[1] - indexes[0]; # stores length of removed region
        self.origin = self.origin[:indexes[0]] + self.origin[indexes[1]:]; # removes from main sequence
        i = 0; # variable iterating over features list.
        maxIter = len(self.features); # sets max iterations. May change if annotations are deleted, which is why we use this weird while numeric loop instead of a for loop with the elements in the iterating variable.
        while i < maxIter: # loops through all annotations editing index numbers if necessary. Iterates numbers and not elements in array to allow deletion of elements without using remove() method (which is apparently glitchy).
            ann = self.features[i]; # gets this iteration's annotation
            if len(ann.seq) > 0: # checks whether annotation actually contains information
                if ann.index[0] >= indexes[1]: # if this annotation's start index is greater or equal to the end of the removed region,
                    ann.index[0] = ann.index[0] - lenRem; # shift the annotation's start index
                    ann.index[1] = ann.index[1] - lenRem; # shift the annotation's end index
                elif ann.index[1] > indexes[1]: # else if the annotation's end index is strictly greater than the end of the removed portion, it means a part of the annotated region is removed.
                    if ann.index[0] >= indexes[0]: # if the start of the annotated region is inside the eliminated sequence, it means the part removed from the annotated region includes the starting point (trimmed from left)
                        lenTrim = indexes[1] - ann.index[0]; # stores length of sequence removed from this region (overlap between removal region and annotated region)
                        ann.index[0] = indexes[0]; # shift the annotation's start index to the excision point
                        ann.index[1] = ann.index[1] - lenRem; # shift the annotation's end index
                        ann.seq = ann.seq[lenTrim:]; # updates annotation sequence
                        r.append(ann); # add to return list
                    else: # if not, it means a middle part of the region was removed
                        ann.index[1] = ann.index[1] - lenRem; # shift the annotation's end index
                        indexInSeq = indexes[0]-ann.index[0]; # index of start of removed portion inside annotation sequence is stored
                        ann.seq = ann.seq[:indexInSeq] + ann.seq[(indexInSeq+lenRem):]; # updates annotation sequence
                        r.append(ann); # add to return list
                elif ann.index[0] >= indexes[0]: # else if the annotation's start index is greater than or equal to the start of the removed sequence, then the totality of the annotated region was removed.
                    ann.index = []; # set indexes to none
                    ann.seq = ""; # set sequence to empty
                    r.append(ann); # add to return list
                    self.features = self.features[:(i)] + self.features[(i+1):]; # remove it from list of annotations
                    i = i - 1; # adjusts iterator accordingly
                    maxIter = maxIter - 1; # adjusts maxIter accordingly
                elif ann.index[1] > indexes[0]: # if the end of the annotated region is inside the eliminated sequence, it means the part removed from the annotated region includes the ending point (trimmed from right)
                    lenTrim = ann.index[1] - indexes[0]; # stores length of sequence removed from this region (overlap between removal region and annotated region)
                    ann.index[1] = indexes[0]; # shift the annotation's end index to the start of the excision point
                    ann.seq = ann.seq[:(len(ann.seq)-lenTrim)]; # updates annotation sequence
                    r.append(ann); # add to return list
                # the final case is that the whole annotated region is upstream of the removed sequence, where nothing changes.

            i = i + 1; # advances iterator


        return r; # returns variable


    """
    Inserts new GenBank object into GenBank object, modifies annotation indexes
    accordingly. The new sequence should be in 5' to 3', positive sense.
    Returns list of annotations whose sequence was altered by the elimination if
    it occurs accross one or more annotated regions, returns empty variable
    otherwise. The returned annotation contains the new sequence, not the
    original one, so it can be removed from the features list if necessary.
    """
    def insertGenBank(self, gb, index):
        gb_copy = deepcopy(gb); # copy of sequence to be inserted
        r = self.insertSeq(gb_copy.origin,index); # inserts new sequence without annotations into this GenBank object, modifies this object's annotations
        gb_copy.insertSeq(self.origin[:index],0); # inserts this object's sequence up until the insertion point into other file
        self.features = self.features + gb_copy.features; # now that all annotations have the same coordinate system, save the annotations onto this object's feature lists
        gb_copy = None; # delete copy

        return r; # returns list of annotations whose sequence was altered by the elimination if it occurs accross one or more annotated regions


    """
    Returns reverse complement of GenBank object. Main sequence is reverse
    complement of original, sequences of annotations are the same as the
    original ones, but indexes are reverse complemented complement(X..Y). Adds
    "Reverse complement of " to sequence definition.
    """
    def revComp(self):
        r = deepcopy(self); # initialize return variable
        r.origin = revComp(r.origin); # revComps the main sequence
        r.definition = "Reverse complement of " + r.definition; # Adds mark to sequence definition
        r.definition.replace("Reverse complement of Reverse complement of ",""); # in case the object had been flipped previously
        for ann in r.features: # iterates over features list
            ann.comp = not ann.comp; # switch orientation
            ann.index = [len(r.origin)-ann.index[1], len(r.origin)-ann.index[0]]; # adjusts annotation indexes according to new orientation

        return r; # returns variable


    """
    Returns list of annotations in features list whose label contains the
    given search term.
    """
    def findAnnsLabel(self, searchTerm, suppressError=False):
        annList = []; # stores found gene annotations
        for ann in self.features: # iterates until end of list reached
            if ann.label.find(searchTerm) > -1: # if this annotation contains the term being searched for,
                annList.append(ann); # add annotation to list

        if len(annList) == 0 and not suppressError: # if annotation still not found and not suppressing error messages,
            print("ERROR: Annotations with '" + searchTerm + "' in label not found in sequence " + self.name); # Report error

        return annList; # returns list


    """
    Returns list of annotations in features list whose type contains the
    given search term.
    """
    def findAnnsType(self, searchTerm, suppressError=False):
        annList = []; # stores found gene annotations
        for ann in self.features: # iterates until end of list reached
            if ann.type.find(searchTerm) > -1: # if this annotation contains the term being searched for,
                annList.append(ann); # add annotation to list

        if len(annList) == 0 and not suppressError: # if annotation still not found and not suppressing error messages,
            print("ERROR: Annotations with '" + searchTerm + "' in type not found in sequence " + self.name); # Report error

        return annList; # returns list

    """
    Returns true if given index is inside a gene exon (assuming either exons or
    introns are annotated)
    """
    def checkInExon(self,pIndex): # checks if pIndex is inside an exon in gene
        insideGene = False; # Boolean stores whether pIndex is inside a gene
        insideExon = False; # Boolean stores whether pIndex is inside an exon
        exonsAnnotated = False; # stores whether introns are annotated in this GB file
        for ann in self.features: # loop through all annotations
            if  ann.type == "gene" and ann.index[0] <= pIndex <= ann.index[1]: # if this annotation is a gene and pIndex is inside this gene,
                insideGene = True; # set inside this gene
            elif ann.type == "exon": # if this annotation is an exon,
                exonsAnnotated = True; # exons are annotated
                if ann.index[0] <= pIndex <= ann.index[1]: # if pIndex is inside this exon,
                    insideExon = True; # sets insideExon



        if insideGene and not exonsAnnotated: # if inside gene with no exons annotated,
            insideExon = True; # Assume inside exon

        return insideExon; # return Boolean

    """
    Returns true if given index is inside a CDS with a given label (assuming
    CDS are annotated)
    """
    def checkInCDS(self,pIndex,label): # checks if pIndex is inside an exon in gene
        insideGene = False; # Boolean stores whether pIndex is inside a gene
        insideExon = False; # Boolean stores whether pIndex is inside an exon
        exonsAnnotated = False; # stores whether introns are annotated in this GB file
        for ann in self.features: # loop through all annotations
            if  ann.type == "gene" and ann.index[0] <= pIndex <= ann.index[1]: # if this annotation is a gene and pIndex is inside this gene,
                insideGene = True; # set inside this gene
            elif ann.type == "CDS": # if this annotation is a CDS,
                exonsAnnotated = True; # exons are annotated
                if ann.index[0] <= pIndex <= ann.index[1] and ann.label.find(label)>-1: # if pIndex is inside this CDS and CDS contains right label,
                    insideExon = True; # sets insideExon



        if insideGene and not exonsAnnotated: # if inside gene with no exons annotated,
            insideExon = True; # Assume inside exon

        return insideExon; # return Boolean


    """
    Loops over all annotations changing the color to the one passed as an
    argument. If pCol is random, randomizes all colors.
    """
    def setAllColors(self,pCol):
        if pCol == "random": # if setting random,
            for ann in self.features: # loop through all annotations
                ann.color = randHex(); # set random color
        else: # if using a set color,
            for ann in self.features: # loop through all annotations
                ann.color = pCol; # set color
