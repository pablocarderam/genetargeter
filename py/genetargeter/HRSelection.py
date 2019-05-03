

from py.utils.BioUtils import *; # Imports utils
from py.utils.GenBankToolbox import *; # Imports utils
from py.genetargeter.constants import *; # Imports constants

"""
Selects an appropriate HR for the given gene. GenBank sequence object given as
argument (geneGB) must have an annotation with label=gene, at least one
annotation with gRNA in label (or some similar indicator), and the length of the
gene must be at least lengthLHR[2]. lengthLHR is [min, preferred, max] length in
bp of LHR. minTmEnds is min melting temp in extremes of LHR. endsLength is
length of extremes of the LHR contemplated in analysis. optimizeRange is a list
with the index positions relative to the start of the LHR over which the start
of the LHR is optimized once it is chosen. maxDistanceFromGRNA is maximum
distance in bp between end of LHR and start of gRNA.
LHR: Left Homologous Region used for chromosomal integration by homologous
recombination during repair.
"""
def chooseHR(geneGB, gene, doingHR='LHR', targetExtreme='end', lengthHR=[450,500,750], minTmEnds=59, endsLength=40, codingGene=True, gBlockDefault=True, minGBlockSize=125, optimizeRange=20, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI]):
    log = "" # init log
    gRNAs = [] # List of all gRNAs
    gRNAExt = GenBankAnn() # init var to hold gRNA
    if len(gRNAs) == 0: # if no gRNAs found
        gRNAs = geneGB.findAnnsLabel("gRNA") # List of all gRNAs
        gRNAExt = gRNAs[0] # will store gRNA most extreme
        for g in gRNAs: # loops across gRNAs
            if ((doingHR=='LHR')*2-1) * g.index[0] < ((doingHR=='LHR')*2-1) * gRNAExt.index[0]: # if more extreme (upstream for end targeting, downsteam for beginning targeting),
                gRNAExt = g # replace as most extreme



    genes = geneGB.findAnnsType("gene") # list of all genes, used to verify LHR doesn't start inside any genes (truncating them)

    # Declare variables and check assertions
    seqBeg = 0
    seqEnd = len(geneGB.origin)
    lenMin = lengthHR[0]
    lenMax = lengthHR[2]
    genBeg = gene.index[0]
    genEnd = gene.index[1]
    nxtGen = min( [ seqEnd ] + [ (genEnd>g.index[0]) * seqEnd + g.index[0] for g in genes ] ) # start of next gene downstream
    prvGen = max( [ seqBeg ] + [ (genBeg>g.index[1]) * g.index[1] for g in genes ] ) # end of previous gene upstream
    gRNAEx = gRNAExt.index[doingHR == 'RHR'] # position of most extreme gRNA to be avoided

    if not ( (seqBeg<=lenMin) and (lenMin<=lenMax) and (lenMax<=genBeg) and (genBeg<=genEnd) and (genEnd<=seqEnd) and (seqEnd-genEnd>=lenMax) ) : # If assertions don't hold,
        log = log + "\nERROR: Not enough space on either side of gene for maximum HR size, or you mixed up min and max HR lengths. Aborting." + "\n" # give a warning
        return {"out":None, "log":log}

    # Initialize region

    regBeg = max( (targetExtreme == 'end') * (genBeg-lenMax-1), prvGen ) if doingHR == 'LHR' else max( (targetExtreme == 'end') * genEnd, gRNAEx, min(gBlockDefault*(genBeg+minGBlockSize),nxtGen-lenMin) ) # beginning of possible LHR or RHR region
    regEnd = min( genEnd-3*codingGene, gRNAEx, max(gBlockDefault*(gRNAEx<genEnd)*(targetExtreme == 'end')*(genEnd-minGBlockSize-3), (targetExtreme == 'end')*seqEnd + genBeg) ) if doingHR == 'LHR' else min( nxtGen, seqEnd ) # end of possible LHR or RHR region

    if regBeg > seqEnd-lenMin or regEnd < seqBeg+lenMin: # if not enough space for the HR on the sequence file,
        log = log + "\nERROR: Not enough space on either side of gene for HR chosen. Aborting." + "\n" # give an error
        return {"out":None, "log":log}

    # Partitioning
    cutSites = [] # will store cut site indeces
    for site in filterCutSites: # for every cut site,
        cutSites = cutSites + findMotif( geneGB.origin[regBeg:regEnd], site ) # store all its occurrences

    cutSites = [ s+regBeg for s in sorted(cutSites) ] # sort by ascending index, add back beginning of region index to find global coordinates
    regBegArr = [regBeg] + [ s+6 for s in cutSites ] # start of possible HR regions, get rid of 6 bp from cutsite and assume 0.02% probability of reforming site doesn't happen :P
    regEndArr = cutSites + [regEnd] # Ends of possible HR regions

    regIdxArr = [] # will store indexes of possible HR regions
    for i in range( len(regBegArr) ): # for every partition,
        if abs( regBegArr[i]-regEndArr[i] ) >= lenMin and (doingHR == 'LHR')*genBeg <= regEndArr[i]: # if enough space for HR in partition and LHR end still within bounds,
            regIdxArr.append( [ regBegArr[i],regEndArr[i] ] ) # add to valid partition index array


    if len(regIdxArr) < 1: # if no valid partitions,
        regIdxArr = [ [ regBeg, regEnd ] ] # default to region
        log = log + "\nWarning: No "+doingHR+" without restriction enzyme cut sites inside it \nand without excluding a downstream gene found, check output manually!" + "\n" # give a warning


    # Find beginning and ending indeces
    begIntArr = [ [ regIdxArr[i][0], regIdxArr[i][1]-lenMin ] for i in range(len(regIdxArr)) ] # intervals to search for beginnings
    endIntArr = [ [ regIdxArr[i][0]+lenMin, regIdxArr[i][1] ] for i in range(len(regIdxArr)) ] # intervals to search for ends
    terIdxArr = [ begIntArr, endIntArr ] # contains array of possible beginning search intervals, and array of possible ending search intervals

    # Iteration search
    step = -1 if doingHR == 'LHR' else 1 # step direction for LHR is -1, 1 for RHR
    ss = step * (step<0) # search start position in partitions, start from end (-1) if doing LHR or from start (0) if RHR
    a = terIdxArr[ss+1][ss][ss] # for LHR, start anchor terminal search as the last ending of the final partition; for RHR, start with the first beginning of the first partition
    b = terIdxArr[ss][ss][ss] # for LHR, start oppos. terminal search as the last beginning of the final partition; for RHR, start with the first ending of the first partition
    #a = a if a>=0 else a + len( geneGB.origin ) # if wrapping around back of seq, set to positive indexing
    #b = b if b>=0 else b + len( geneGB.origin ) # if wrapping around back of seq, set to positive indexing
    p = ss if ss>=0 else ss + len( begIntArr )  # will keep track of partition if wrapping around back of partitions, set to positive indexing

    i = [ min(a,b), max(a,b) ] # search each terminus start position in sequence terminus, start from end if doing LHR or from start if RHR, flipped for loop
    bestHR = [ i,
        [ isTricky( geneGB.origin[ i[0]:i[0]+endsLength ] ),
            isTricky( geneGB.origin[ i[1]-endsLength:i[1] ] ) ],
        [ meltingTemp( geneGB.origin[ i[0]:i[0]+endsLength ] ),
            meltingTemp( geneGB.origin[ i[1]-endsLength:i[1] ] ) ] ] # save indeces, synthesis problems, and melting temperatures of best HR so far, initial guess as default

    satisfiedHR = not ( ( bestHR[1][0] + bestHR[1][1] + ( bestHR[2][0] < minTmEnds ) + ( bestHR[2][1] < minTmEnds ) ) > 0 ) # keep track of whether or not adequate HR has been found

    while 0 <= p and p < len(regIdxArr) and not satisfiedHR: # while partitions not exhausted and no good HR found,
        bestInPart = bestHR # will store best anchors in this partition
        ter = (step<0) # which terminus we're adjusting, 0 is beginning and 1 is ending, start with ending if LHR and beginning if RHR (will get switched at start of next loop)
        a = terIdxArr[ss+1][p][ss] # for LHR, start anchor terminal search as the last ending of the final partition; for RHR, start with the first beginning of the first partition
        b = terIdxArr[ss][p][ss] # for LHR, start oppos. terminal search as the last beginning of the final partition; for RHR, start with the first ending of the first partition
        i = [ min(a,b), max(a,b) ] # search each terminus start position in sequence terminus, start from end if doing LHR or from start if RHR, flipped for loop

        bestInTer = [ i, [True,True], [0,0] ] # will keep track of last two terminus searches, initialize as default
        while ( terIdxArr[0][p][0] <= i[0] and i[0] <= terIdxArr[0][p][1] ) and ( terIdxArr[1][p][0] <= i[1] and i[1] <= terIdxArr[1][p][1] ) and i[0] < nxtGen and not satisfiedHR: # while not having exhausted this terminus and not having found a good HR
            if lenMin > i[1]-i[0] or i[1]-i[0] > lenMax: # if size is wrong,
                while ( lenMin > i[1]-i[0] or i[1]-i[0] > lenMax ) and ( terIdxArr[ter][p][0] <= i[ter] and i[ter] <= terIdxArr[ter][p][1] ) and i[0] < nxtGen: # while not within size bounds and still within search range,
                    i[ter] += step # move terminus until within size constraints


            bestInTer[0][ter], bestInTer[1][ter], bestInTer[2][ter] = i[ter], True, 0 # will store best in this terminus search, initialize as default
            satisfiedTer = not ( ( bestInTer[1][ter] + ( bestInTer[2][ter] < minTmEnds ) ) ) # keep track of whether or not adequate HR has been found
            while ( terIdxArr[ter][p][0] <= i[ter] and i[ter] <= terIdxArr[ter][p][1] ) and i[0] < nxtGen and not satisfiedTer: # while within search range and no good terminus found,
                extReg = geneGB.origin[ min(i[ter],i[ter]+(2*(not ter)-1)*endsLength):max(i[ter],i[ter]+(2*(not ter)-1)*endsLength) ] # extreme region to be analyzed
                tricky = isTricky( extReg ) # true if this terminus contains homopolymers or AT repeats
                tm = meltingTemp( extReg ) # Tm for this terminus
                lenHR = i[1] - i[0] # length of HR as it stands
                inIntronWhenNotSupposedTo = ((not geneGB.checkInExon(i[0]) and targetExtreme!='end' and doingHR=='RHR') or (not geneGB.checkInExon(i[1]) and targetExtreme=='end' and doingHR=='LHR')) # shouldn't be in intron if LHR, targeting 3', and end terminus or if RHR, targeting 5', and beginning terminus
                if (not tricky) and tm >= minTmEnds and not inIntronWhenNotSupposedTo: # if sufficiently good and not in intron when not supposed to,
                    # optimize ends
                    ti = i[ter] # test more indeces to optimize
                    bi = i[ter] # best i
                    while ( max(terIdxArr[ter][p][0],i[ter]-optimizeRange) <= ti and ti <= min(terIdxArr[ter][p][1],nxtGen,i[ter]+optimizeRange) ) : # while we can still optimize within bounds,
                        ti += step # advance searcher
                        t_extReg = geneGB.origin[ min(ti,ti+(2*(not ter)-1)*endsLength):max(ti,ti+(2*(not ter)-1)*endsLength) ] # extreme region to be analyzed
                        t_tm = meltingTemp( t_extReg ) # Tm for this terminus
                        if t_tm > tm and not isTricky( t_extReg ) and lenMax >= abs(i[not ter] - ti) >= lenMin: # if this end point has a better Tm, and is still within bounds
                            bi = ti; # make this the ending position
                            tm = t_tm # record this tm as best


                    bestInTer[0][ter], bestInTer[1][ter], bestInTer[2][ter] = bi, tricky, tm # save as best in this terminus search
                    satisfiedTer = True # satisfied with this terminus
                    if ( lenMin <= lenHR and lenHR <= lenMax ):
                        bestInPart[0][ter], bestInPart[1][ter], bestInPart[2][ter] = i[ter], tricky, tm
                        if ( ( bestInTer[1][0]+bestInTer[1][1] <= bestHR[1][0]+bestHR[1][1] ) and ( ( ( ( bestInTer[2][0]<minTmEnds ) + ( bestInTer[2][0]<minTmEnds ) < ( bestHR[2][0]<minTmEnds ) + ( bestHR[2][0]<minTmEnds ) ) or ( bestInTer[2][0]+bestInTer[2][1] >= bestHR[2][0]+bestHR[2][1] ) ) and ( bestInTer[1][0]+bestInTer[1][1] == bestHR[1][0]+bestHR[1][1] ) ) ) : # if size is right and right length and better than best overall
                            bestHR = bestInPart # save as best overall
                            satisfiedHR = not ( ( bestHR[1][0] + bestHR[1][1] + ( bestHR[2][0] < minTmEnds ) + ( bestHR[2][1] < minTmEnds ) ) ) # keep track of whether or not adequate HR has been found


                else: # if not sufficiently good
                    if (tricky < bestInTer[1] or ( tricky == bestInTer[1] and tm > bestInTer[2] )) and not inIntronWhenNotSupposedTo: # if at least better than best of this terminus and not in exon when not supposed to be,
                        bestInTer[0][ter], bestInTer[1][ter], bestInTer[2][ter] = i[ter], tricky, tm # save as best in this terminus search
                        if ( lenMin <= lenHR and lenHR <= lenMax ):
                            bestInPart[0][ter], bestInPart[1][ter], bestInPart[2][ter] = i[ter], tricky, tm
                            if ( ( bestInTer[1][0]+bestInTer[1][1] <= bestHR[1][0]+bestHR[1][1] ) and ( ( ( ( bestInTer[2][0]<minTmEnds ) + ( bestInTer[2][0]<minTmEnds ) < ( bestHR[2][0]<minTmEnds ) + ( bestHR[2][0]<minTmEnds ) ) or ( bestInTer[2][0]+bestInTer[2][1] >= bestHR[2][0]+bestHR[2][1] ) ) and ( bestInTer[1][0]+bestInTer[1][1] == bestHR[1][0]+bestHR[1][1] ) ) ) : # if size is right and right length and better than best overall
                                bestHR = bestInPart # save as best overall



                    i[ter] += step # advance or decrease base index according to search direction


            ter = not ter # switch which terminus we're adjusting, 0 is beginning and 1 is ending

        p += step # advance or decrease partition according to search direction

    if not satisfiedHR:
        log = log + "Warning: No "+doingHR+" without restriction enzyme cut sites inside it, \nwithout excluding a downstream gene found, with adequate Tm at ends and \nwithout excessive homopolymers at ends found, check output manually!" + "\n" # give a warning

    log = log + doingHR + " for gene " + gene.label + " selected.\n\n"; # logs this process finished
    HR = GenBankAnn( gene.label + " " + doingHR, "misc_feature", geneGB.origin[ bestHR[0][0]:bestHR[0][1] ], False, [ bestHR[0][0],bestHR[0][1] ], annColors[doingHR+'Color'] ) # creates GenBankAnn object to hold LHR

    return {"out":HR, "log":log} # returns LHR GenBankAnn object
