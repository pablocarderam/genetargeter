# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 19:24:17 2016
Classes useful in rosalind bioinformatic algorithms
@author: Pablo the awesome molecular jedi
"""

from BioUtils import *; # Imports utils
from copy import deepcopy; # Import object copying methods for deep copies, used in creating revComp of GenBank objects

class GenBankAnn(object):
    """Stores some GenBank format information for an annotation."""
    def __init__(self, pLab="", pType="", pSeq="", pComp=False, pIndex=[]): # class constructor
        self.label = pLab # annotation label
        self.type = pType; # type of annotation
        self.seq = pSeq; # sequence of annotation
        self.comp = pComp; # false if on positive strand, true if on complementary strand
        self.index = pIndex; # list with region start and end indexes IN PYTHON INDEXING (starting on 0)


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

        d = d.split("\n"); # Split into lines. d is now a list of strings corresponding to the lines of the file.

        i = 0; # indexes through lines
        w = d[i].split(); # splits line into words

        while w[0] != "LOCUS" and i < len(d): # iterates until LOCUS found or end of document reached
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
            self.definition = d[i][10].lstrip(); # saves definition information on this line after "DEFINITION" string and ensuing whitespace.
            while w[len(w)-1].find(".") > -1 and i < len(d): # iterates until last word in line contains period, or end of document reached
                self.definition = d[i]; # save this line as definition information
                i = i+1; # advance counter
                w = d[i].split(); # splits line into words

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
                if d[i].lstrip()[0] != "/": # if the line does not start with "/", it means it's a new annotation. lstrip() removes left whitespace.
                    newAnn = GenBankAnn("","","",False,[]); # resets newAnn variable
                    newAnn.type = w[0]; # stores annotation type
                    posRawString = w[1].split("("); # splits second word in w by "(" to see if on complementary strand
                    posString = posRawString[0]; # assigns posString to the indexes string in the case that the annotation is on the positive strand
                    if len(posRawString) > 1: # if split into more than one piece, it means it's on the complementary strand
                        newAnn.comp = True; # set comp to true
                        posString = posRawString[1][:(len(posRawString[1])-1)]; # extracts index information from raw string

                    newAnn.index = [int(index) for index in posString.split("..")]; # takes the index string, splits it into two strings separated by ".." and converts them both to ints.
                    newAnn.index[0] = newAnn.index[0] - 1; # changes to Python indexing. Upper limit does not have to be changed since Python does not include it in indexing, while GenBank does.
                else:
                    p = d[i].lstrip().split("="); # splits line into a property's name and value
                    if p[0] == "/gene" or p[0] == "/label": # if the property is a gene or label...
                        newAnn.label = p[1][1:(len(p[1])-1)]; # assigns the label (extracting the label from string p)
                        if newAnn.comp: # if on complementary strand,
                            newAnn.seq = revComp(self.origin[(newAnn.index[0]):newAnn.index[1]]); # assigns the sequence of the annotated region to the complement. To find the sequence of the annotated region, indexes in the full sequence according to the position previously extracted.
                        else: # if on positive strand,
                            newAnn.seq = self.origin[(newAnn.index[0]):newAnn.index[1]]; # assigns the sequence of the annotated region. To find the sequence of the annotated region, indexes in the full sequence according to the position previously extracted.



                if len(newAnn.label)*len(newAnn.type)*len(newAnn.seq)*len(newAnn.index) != 0: # if all propoerties in newAnn have been filled in,
                    self.features.append(newAnn); # adds the new annotation
                    newAnn = GenBankAnn("","","",False,[]); # resets newAnn variable

                i = i+1; # advance counter
                w = d[i].split(); # splits line into words


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
        defOut = "\n".join([d[i:i+80] for i in range(0, len(d), 80)]); # introduces line separators into string

        outStr = outStr + defOut + "\nFEATURES             Location/Qualifiers\n"; # write locus and definition information to file, writes a new line to start features list

        for ann in self.features:
            annSeq = ann.seq; # default case: positive strand sequence
            if ann.comp: # if on complementary strand,
                annSeq = revComp(ann.seq); # save as rev comp

            if self.origin[ann.index[0]:ann.index[1]] == annSeq: # if the annotation sequence was found in the region where it should be,
                startIndex = ann.index[0] + 1; # advances index (GenBank indexing starts at 1)
                endIndex = ann.index[1]; # stores end of sequence (GenBank starts at 1, but we don't shift this index because GenBank indexing includes final index)
                indexStr = str(startIndex) + ".." + str(endIndex); # save annotation indexes as string
                if ann.comp:
                    indexStr = "complement(" + indexStr + ")"; # if on complementary strand, rewrite indexStr with this information

                spacing = "".join([" " for i in range(16-len(ann.type))]); # creates enough spacing to complete 16 characters between tab and indexes
                outStr = outStr + "     " + ann.type + spacing + indexStr + "\n"; # write first line of annotation format
                outStr = outStr + '                     /label="' + ann.label + '"\n'; # write label line
                col = randHex(); # generates a random color for annotation
                outStr = outStr + '                     /ApEinfo_revcolor=' + col + '\n'; # write reverse color line
                outStr = outStr + '                     /ApEinfo_fwdcolor=' + col + '\n'; # write forward color line
            else: # if it wasn't found where it should be,
                print "ERROR: Annotation sequence " + ann.label + " not found in main sequence."; # warn of error.

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
                        indexInSeq = ann.indexes[0]-ann.index[0]; # stores index of start of removed portion inside annotation sequence
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
                    ann.index[1] = indexes[0]; # shift the annotation's end index to the start of the excision point
                    ann.seq = ann.seq[:(len(ann.seq)-lenTrim)]; # updates annotation sequence
                    r.append(ann); # add to return list
                # the final case is that the whole annotated region is upstream of the removed sequence, where nothing changes.

            i = i + 1; # advances iterator


        return r; # returns variable


    """
    Returns reverse complement of GenBank object. Main sequence is reverse
    complement of original, sequences of annotations are the same as the
    original ones, but indexes are reverse complemented complement(X..Y). Adds
    "Reverse complement of " to sequence definition.
    """
    def revComp(self):
        r = deepcopy(self); # initialize return variable
        r.comp = not r.comp; # switch orientation
        r.origin = revComp(r.origin); # revComps the main sequence
        r.definition = "Reverse complement of " + r.definition; # Adds mark to sequence definition
        for ann in r.features: # iterates over features list
            ann.index = [len(r.origin)-r.index[0]-1, len(r.origin)-r.index[1]-1]; # adjusts annotation indexes according to new orientation

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
            print "ERROR: Annotations with '" + searchTerm + "' in label not found in sequence " + self.name; # Report error

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
            print "ERROR: Annotations with '" + searchTerm + "' in type not found in sequence " + self.name; # Report error

        return annList; # returns list
