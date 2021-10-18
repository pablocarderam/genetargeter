
from builtins import str
import operator; # Needed for sorting
import copy; # Import object copying methods for deep copies

from py.utils.BioUtils import *; # Imports utils
from py.utils.GenBankToolbox import *; # Imports utils
from py.genetargeter.constants import *; # Imports constants
from gRNAScores.gRNAScoring import * # import methods for gRNA scoring
from gRNAScores.gRNAScoreDB import * # import methods for gRNA scoring from precalculated DB

"""
Selects an appropriate gRNA for the given gene. GenBank gene sequence given as
argument must have an annotation with "gRNA 1" as a part of its label, and the
gene must be at least searchRange[0] bp long. The file must include at least
searchRange[1] bp of 3' UTR. searchRange is a list indicating search start and
end indexes, counted with the last bp in the gene's stop codon as index 0.
target3Prime is false if PAM sequence is at the 5' end of gRNA, true if at 3'.
gRNA: guide RNA used by CRISPR enzyme.
Note: argument searchRange is in format [inside_gene,outside_gene] whether
targeting 3' or 5' end.
"""
def chooseGRNA(geneGB, gene, searchRange=[-700,125], searchRangeNonCoding=550, PAM="NGG", minGCContent=0.3, minOnTargetScore=25, minOffTargetScore=75, maxOffTargetHitScore=35, onTargetMethod="azimuth", offTargetMethod="hsu", gLength=20, maxDistanceBetweenGRNAS=50, enzyme="Cas9", gBlockDefault=True, maxTier1GBlockSize=500, gBlockOverlapSize=40, codingGene=True, closestGene=-1, target3Prime=True, targetCenter=False, targetRegionOverride=False, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI,cut_NheI]): # could've been useful at some point: http://grna.ctegd.uga.edu/ http://www.broadinstitute.org/rnai/public/software/sgrna-scoring-help http://crispr.mit.edu/about
    if closestGene < 0: # if closestGene parameter is default,
        closestGene = len(geneGB.origin); # set to total length of gene file as default

    log = "Choosing gRNA with PAM sequence " + PAM + " for use with enzyme " + enzyme + '\n\n'; # init log
    if not codingGene and searchRange[0] < 0: # if gene is non protein-coding and part of the search region is inside the gene,
        log += "Since this looks like a non protein-coding gene, the gRNA search range will be shifted to be entirely outside of the end of the gene.\n\n"; # Notate change in search range
        searchRange = [0,searchRangeNonCoding]; # move search range outside gene

    gRNATable = []; # will store information on each gRNA evaluated. Format: Label, Status, Enzyme, Position, Strand, GC_content, On-target_score, On-target_method, Aggregated_off-target_score, Max_pairwise_off-target_score, Off-target_method, >9_consecutive_A/T, 4-Homopolymer, Triple_T, Sequence, Recoded_sequence, Recoded_sequence_pairwise_off-target_score
    geneList = geneGB.findAnnsLabel(gene.label); # stores gene annotation


    for g in geneList: # search for right gene
        if g.label == gene.label and g.type == "gene" or g.type == "mRNA": # if found,
            gene = g; # save it

    # This line of code increases the search range if there are introns. We decided to roll this back for the time being.
    # searchRange[0] = adjustExonicRange(geneGB, gene, target3Prime, searchRange) # expand search range in cases where there are lots of introns or something
    gRNAs = []; # list of GenBankAnn objects for gRNAs.
    backupGRNAs = []; # stores gRNAs that fail off-target score as possible backups
    gRNAExtreme = GenBankAnn(); # will store gRNA most upstream

    if targetCenter: # if knocking out and targeting center of gene,
        center = int((gene.index[1]-gene.index[0])/2) # find center of gene
        searchRangeLength = int((searchRange[1]-searchRange[0])/2) # find length of search range
        searchRange = [-center-searchRangeLength,-center+searchRangeLength] # center search range on center; flip order and sign since it'll be flipped due to 3' targeting when KO

    searchRange[0] = max(searchRange[0],-1*len(geneGB.origin)); # adjusts search range in case it goes beyond the total length of the gene file (to the left) (SHOULD WORK IN 5' CASE TOO)
    searchRange[1] = min(abs(closestGene-gene.index[target3Prime]),searchRange[1]); # set end of search range as start of next downstream gene if that happens before the current end of search range (SHOULD WORK IN 5' CASE TOO)

    pamSeqs = ambiguousSeqs(PAM); # store actual PAM sequences
    extGRNASeqIndexes = []; # stores indexes of extended gRNA sequence (NNNN-gRNA 20mer-PAM 3mer-NNN for Cas9, PAM 4mer-gRNA 23mer-NNN for Cas12) relative to PAM start point.
    realGRNASeqIndexes = []; # stores indexes of actual gRNA sequence (gRNA 20mer for Cas9, gRNA 23mer for Cas12) relative to PAM start point.
    pamIndexes = []; # stores indexes of PAM sequences within extended gRNA

    if enzyme == "Cas9": # if enzyme is Cas9,
        extGRNASeqIndexes = [-24,6,-3,27]; # Start and end indexes for NNNN-gRNA 20mer-PAM 3mer-NNN extended sequence relative to PAM start position. Second two numbers are for rev comp strand.
        realGRNASeqIndexes = [4,24]; # Start and end indexes for gRNA 20mer sequence within extended sequence.
        pamIndexes = [24,24+len(PAM)]; # Stores PAM seq indexes, including Ns and Vs
    elif enzyme == "Cas12": # if enzyme is Cas12,
        extGRNASeqIndexes = [0,30,-26,4]; # Start and end indexes for PAM 4mer-gRNA 23mer extended sequence relative to PAM start position. Second two numbers are for rev comp strand.
        realGRNASeqIndexes = [len(PAM),23+len(PAM)]; # Start and end indexes for PAM gRNA 23mer sequence within extended sequence.
        pamIndexes = [0,len(PAM)]; # Stores PAM seq indexes, including Ns and Vs
        # assign values

    if len(gene.label) > 0: # if gene found,
        searchStart = gene.index[target3Prime] + (-1+(2*target3Prime)) * searchRange[(not target3Prime)]; # Start searching for gRNAs in a position relative to gene start or end point
        searchEnd = gene.index[target3Prime] + (-1+(2*target3Prime)) * searchRange[target3Prime]; # Finish searching for gRNAs in a position relative to gene start or end point

        # Override for annotated target regions
        if targetRegionOverride:
            targetRegion = geneGB.findAnnsLabel("Target Region");
            if len(targetRegion) > 0:
                searchStart = targetRegion[0].index[0]; # Start searching for gRNAs in target region
                searchEnd = targetRegion[0].index[1]; # Finish searching for gRNAs in target region

        searchSeq = geneGB.origin[searchStart:searchEnd].upper(); # get sequence where gRNAs will be searched for, centered around start or end of gene
        i = len(PAM) + 1; # indexer used to search for PAMs. Start searching for gRNAs by searching for PAM in searchSeq. Start at upstream end, go downstream
        if target3Prime: # if going for 3',
            i = len(searchSeq)-extGRNASeqIndexes[1]+1; # start searching for gRNAs by searching for PAM in searchSeq. Start at downstream end (allow five bases downstream to accomodate on-target 30-mer gRNA sequence), go upstream

        while i > len(PAM) and i < len(searchSeq)-extGRNASeqIndexes[1]+2 and len(gRNAs) < 3: # iterate through all searchSeq until done with searchSeq (allow enough bases upstream of search end point to accomodate on-target 30-mer minus-strand gRNA sequence) or three candidates found.
            comp = False; # set sense to plus strand
            extGRNASeq = ""; # stores extended sequence sequence (NNNN-gRNA 20mer-PAM 3mer-NNN for Cas9, PAM 4mer-gRNA 23mer for Cas12)
            gRNAIndexes = []; # store gRNA indexes
            if searchSeq[i:i+len(PAM)] in pamSeqs: # if PAM found on plus strand,
                extGRNASeq = searchSeq[i+extGRNASeqIndexes[0]:i+extGRNASeqIndexes[1]]; # store extended sequence
                gRNAIndexes = [searchStart+i+extGRNASeqIndexes[0]+realGRNASeqIndexes[0],searchStart+i+extGRNASeqIndexes[0]+realGRNASeqIndexes[1]]; # store gRNA indexes
            elif searchSeq[i:i+len(PAM)] in [revComp(p) for p in pamSeqs]: # if PAM found on minus strand, +extGRNASeqIndexes[0]+realGRNASeqIndexes[0]
                extGRNASeq = revComp(searchSeq[i+extGRNASeqIndexes[2]:i+extGRNASeqIndexes[3]]); # store extended sequence on comp strand
                gRNAIndexes = [searchStart+i+extGRNASeqIndexes[3]-realGRNASeqIndexes[1],searchStart+i+extGRNASeqIndexes[3]-realGRNASeqIndexes[0]]; # store gRNA indexes on comp strand
                comp = True; # set sense to complementary strand

            if len(extGRNASeq) == extGRNASeqIndexes[1]-extGRNASeqIndexes[0] and not (geneGB.checkInExon(gRNAIndexes[0]) and (gRNAIndexes[1] < gene.index[0] or gene.index[1] < gRNAIndexes[0]) ): # if extended gRNA is right size (doesn't overstep boundaries) and isn't inside a gene further upstream or downstream,
                pamSeq = extGRNASeq[pamIndexes[0]:pamIndexes[1]]; # stores PAM sequence
                gRNASeq = extGRNASeq[realGRNASeqIndexes[0]:realGRNASeqIndexes[1]]; # store actual gRNA seq without PAM
                gc = gcContent(gRNASeq); # store gc content
                strandString = '+'; # stores string denoting strand orientation
                if comp: # if on opposite strand,
                    strandString = '-'; # save that info

                gRNATable.append(['Unlabeled','Rejected (low GC content)',enzyme,str(gRNAIndexes).replace(',',' to'),strandString,str(gc),'Not evaluated',onTargetMethod,'Not evaluated','Not evaluated',offTargetMethod,str(findFirst(gRNASeq.replace('T','A'),"AAAAAAAAAA") > -1),'Not evaluated','Not evaluated',gRNASeq,'Not recoded','-']); # starts storing info
                if gc >= minGCContent*0.8 and findFirst(gRNASeq.replace('T','A'),"AAAAAAAAAA") < 0: # if gc content is acceptable and does not contain 10 or more consecutive As or Ts,
                    gRNAInfo = getGRNAInfoFromDB(extGRNASeq,enzyme); # access gRNA scores from DB
                    onTarget = 0;
                    if len(gRNAInfo) == 0: # if not found in DB
                        onTarget = onTargetScore(extGRNASeq,onTargetMethod); # store on-target score
                    else: # if found in DB,
                        if onTargetMethod == "azimuth": # if using Azimuth (Doench et al., 2016)
                            onTarget = gRNAInfo[1]; # use it
                        # elif onTargetMethod == "ruleset2": # if using Rule Set 2 (Doench et al., 2016)
                        #     onTarget = onTargetScore(extGRNASeq,onTargetMethod); # store on-target score
                        elif onTargetMethod == "cindel": # if using CINDEL (Kim et al., 2017)
                            onTarget = gRNAInfo[2]; # use it


                    gRNATable[len(gRNATable)-1][1] = 'Rejected (low on-target score)'; # Edit this gRNA's status
                    gRNATable[len(gRNATable)-1][6] = str(onTarget); # Edit this gRNA's on-target score
                    if onTarget > minOnTargetScore: # if on-target score is passable,
                        offTargetScores = [];
                        if len(gRNAInfo) == 0: # if not found in DB
                            offTargetScores = offTargetScore(gRNASeq,offTargetMethod,enzyme,pamSeq,PAM); # store off-target score
                        else: # if found in DB,
                            if offTargetMethod == "hsu": # if using Hsu
                                offTargetScores = gRNAInfo[3:6]; # use it
                            elif offTargetMethod == "cfd": # if using CFD
                                offTargetScores = gRNAInfo[9:12]; # use it


                        newGRNA = GenBankAnn(gene.label + " " + enzyme + " gRNA ","misc",gRNASeq,comp,gRNAIndexes,annColors['gRNAColor']); # create annotation object with gRNA information
                        newGRNA.onTarget = onTarget; # Add on-target score as attribute
                        newGRNA.offTarget = offTargetScores; # Add off-target scores as attribute
                        newGRNA.gc = gc; # Add gc content as attribute
                        newGRNA.homopolymer = ( findFirst(gRNASeq,"AAAA") > -1 or findFirst(gRNASeq,"TTTT") > -1 or findFirst(gRNASeq,"CCCC") > -1 or findFirst(gRNASeq,"GGGG") > -1 ); # 4-homopolymers are bad https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1581-4
                        newGRNA.tripleT = (findFirst(gRNASeq,"TTT") > -1); # Triple Ts (triple Us) are bad because they're an RNApol stop codon https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1581-4

                        gRNATable[len(gRNATable)-1][1] = 'Possible backup (low aggregate off-target)'; # Edit this gRNA's status
                        gRNATable[len(gRNATable)-1][8] = str(offTargetScores[0]); # Edit this gRNA's aggregate off-target score
                        gRNATable[len(gRNATable)-1][9] = str(offTargetScores[1]); # Edit this gRNA's max pairwise off-target score
                        gRNATable[len(gRNATable)-1][12] = str(newGRNA.homopolymer); # Edit this gRNA's homopolymer
                        gRNATable[len(gRNATable)-1][13] = str(newGRNA.tripleT); # Edit this gRNA's tripleT
                        gRNATable[len(gRNATable)-1][15] = 'No recoded region necessary'; # Edit this gRNA's recoded status

                        containsCutSite = "" # store info on whether it contains a RE cut site
                        for site in filterCutSites: # for every cut site being filtered
                            if findFirst(site,gRNASeq) > -1 or findFirst(revComp(site),gRNASeq) > -1: # if cut site found,
                                containsCutSite = site; # notate site

                        if gc >= minGCContent and offTargetScores[0] >= minOffTargetScore and offTargetScores[1] <= maxOffTargetHitScore and not newGRNA.homopolymer and not newGRNA.tripleT and containsCutSite == "": # if total off-target score and max hit score are passable, and if no homopolymer or triple T, and if contains no cut sites,
                            gRNAs.append(newGRNA); # add this gRNA's information to list.
                            gRNATable[len(gRNATable)-1][1] = 'Valid'; # Edit this gRNA's status
                        else: # if failed off-target score culling.
                            newGRNA.label = newGRNA.label + "(backup)"; # notify backup on label
                            backupGRNAs.append(newGRNA); # add this gRNA's information to backup list.

                            gRNATable[len(gRNATable)-1][1] = "Possible backup " # contain messages
                            if gc < minGCContent: # if failed GC check,
                                gRNATable[len(gRNATable)-1][1] += '(low GC content); '; # Edit this gRNA's status
                            if offTargetScores[1] > maxOffTargetHitScore: # if failed max pairwise off-target check,
                                gRNATable[len(gRNATable)-1][1] += '(high max pairwise off-target); '; # Edit this gRNA's status
                            if newGRNA.homopolymer: # if failed max pairwise off-target check,
                                gRNATable[len(gRNATable)-1][1] += '(4-homopolymer); '; # Edit this gRNA's status
                            if newGRNA.tripleT: # if failed max pairwise off-target check,
                                gRNATable[len(gRNATable)-1][1] += '(triple T);'; # Edit this gRNA's status
                            if containsCutSite: # if failed cut site check
                                gRNATable[len(gRNATable)-1][1] += '(contains restriction cut site: ' + containsCutSite + ')'; # Edit this gRNA's status



                elif gc < minGCContent*0.8 : # if rejected due to low gc,
                    gRNATable.append(['Unlabeled','Rejected (low GC content)',enzyme,str(gRNAIndexes).replace(',',' to'),strandString,str(gc),'Not evaluated',onTargetMethod,'Not evaluated','Not evaluated',offTargetMethod,str(findFirst(gRNASeq.replace('T','A'),"AAAAAAAAAA") > -1),'Not evaluated','Not evaluated',gRNASeq,'Not recoded','-']); # starts storing info
                else: # if rejected due to low gc,
                    gRNATable.append(['Unlabeled','Rejected (>9 consecutive A/Ts)',enzyme,str(gRNAIndexes).replace(',',' to'),strandString,str(gc),'Not evaluated',onTargetMethod,'Not evaluated','Not evaluated',offTargetMethod,str(findFirst(gRNASeq.replace('T','A'),"AAAAAAAAAA") > -1),'Not evaluated','Not evaluated',gRNASeq,'Not recoded','-']); # starts storing info






            if target3Prime: # if targeting 3' end
                i -= 1; # advance indexer
            else: # if going for 5' end
                i += 1; # advance indexer

        inUTR = False; # will be used to explain status in table

        if len(gRNAs) == 0: # if no gRNAs found,
            log = log + "Warning: no acceptable gRNA with at least " + str(minGCContent*100) + "% GC content, " + str(minOnTargetScore) + " on-target score, " + str(minOffTargetScore) + " off-target total score, no 4-homopolymer sequences, and no TTT sequences found on gene " + gene.label + ".\n" + "Will use backup gRNA with highest GC content, if there are any backups.\n\n"; # add warning to log
        else: # if gRNAs found,
            newList = []; # new list will only contain gRNAs within acceptable range
            if ( gRNAs[0].index[0] > gene.index[1]-(3*codingGene) and target3Prime ) or ( gRNAs[0].index[1] < gene.index[0] and not target3Prime ): # if most downstream gRNA starts after start of stop codon, or most upstream is before gene start
                inUTR = True; # note gRNA is in UTR
                for g in gRNAs: # loop through gRNAs
                    if ( g.index[0] > gene.index[1]-(3*codingGene) and target3Prime ) or ( g.index[1] < gene.index[0] and not target3Prime ): # if this gRNA is still in the UTR or stop codon,
                        newList.append(g); # add it to the new list



            else: # if most downstream gRNA is in gene,
                for g in gRNAs: # loop through gRNAs
                    if targetCenter or (gRNAs[0].index[0] - g.index[0] <= maxDistanceBetweenGRNAS) or (gBlockDefault and gene.index[1]-(3*codingGene) - g.index[0] < maxTier1GBlockSize-2*gBlockOverlapSize): # if this gRNA is within the max distance from the most downstream gRNA or gBlocks are the default and the gRNA is within the cheapest gBlock range from the end of the gene, or if doing a knock-out,
                        newList.append(g); # add it to the new list




            gRNAs = newList; # keep the new list of gRNAs as the main list
            gRNAs.sort(reverse=True,key=lambda g: (g.gc,g.onTarget)); # sorts gRNA list according to custom function (GC content first, on-target score second)
            count = 1; # counter for numbering candidates according to their quality
            for g in gRNAs: # loop through gRNAs ordered by GC content
                g.label = g.label + str(count); # add number to gRNA label
                count += 1; # advance counter

            bestGRNA = gRNAs[0]; # will store best gRNA
            gRNAExtreme = gRNAs[0]; # will find most upstream gRNA
            for g in gRNAs: # loop through gRNAs
                if targetCenter: # if targeting center for KO
                    center = int((gene.index[1]-gene.index[0])/2) # find center of gene
                    if g.index[0]-center < gRNAExtreme.index[0]-center: # if more close to center than previous most close
                        gRNAExtreme = g; # set this gRNA as most close
                elif target3Prime: # if targeting 3',
                    if g.index[0] < gRNAExtreme.index[0]: # if more upstream than previous most upstream
                        gRNAExtreme = g; # set this gRNA as most upstream


                else: # if targeting 5',
                    if g.index[0] > gRNAExtreme.index[0]: # if more downstream than previous most downstream
                        gRNAExtreme = g; # set this gRNA as most downstream


            log = log + str(len(gRNAs)) + " acceptable gRNAs were selected automatically on gene " + gene.label + ". \ngRNA 1 has GC content of " + str(bestGRNA.gc*100) + "%, on-target score of " + str(bestGRNA.onTarget)  + ", \nand aggregated off-target score of " + str(bestGRNA.offTarget[0]) + " (Method: " + offTargetMethod + ", Max. Hit Score: " + str(bestGRNA.offTarget[1]) + ", Num. hits: " + str(bestGRNA.offTarget[2]) + ").\n\n"; # add warning to log

        geneGB.features = geneGB.features + gRNAs; # add gRNAs to gene GenBank object features list
        countBackups = 0; # counts how many backup gRNAs are included
        allGRNAS = gRNAs; # will store both valid and backup gRNAs
        if len(gRNAs) > 0: # if there is at least one gRNA
            for g in backupGRNAs: # for every possible backupGRNAs
                if ( target3Prime and g.index[0] > gRNAExtreme.index[0] ) or ( not target3Prime and g.index[0] < gRNAExtreme.index[0] ): # if this backup is downstream of most extreme gRNA,
                    geneGB.features.append(g); # add to features list
                    allGRNAS.append(g); # add to full gRNA list
                    countBackups +=1; # advances counter


        else: # if there are no gRNAs,
            maxGC = 0; # will track max gc content of gRNA
            maxOnTarget = 0 # will track on-target score of best gRNA
            for g in backupGRNAs: # for every possible backupGRNAs
                geneGB.features.append(g); # add to features list
                allGRNAS.append(g); # add to full gRNA list
                countBackups +=1; # advances counter
                if g.gc > maxGC or (g.gc == maxGC and g.onTarget > maxOnTarget):
                        # if this gRNA has a greater GC content,
                    gRNAExtreme = g; # set as most extreme gRNA
                    maxGC = g.gc # set max gc for evaluation
                    maxOnTarget = g.onTarget # track on-target score of best gRNA
                    bestGRNA = g # set as best gRNA



        for gArr in gRNATable: # loop over all gRNAs in table
            found = False; # determine whether the gRNA is in labeled list
            for g in allGRNAS: # loop over list of labeled gRNAs
                if gArr[14] == g.seq: # if not included in annotations
                    found = True; # set found to True
                    break; # break inner loop


            if not found: # if not included in annotations
                gArr[15] = 'Not recoded'; # Edit this gRNA's recoded status
                if gArr[1][0:8] == "Rejected": # if status is rejected,
                    gArr[1] = gArr[1] + " unlabeled (due to rejection)"; # Add unlabeled explanation to status
                else: # if not rejected,
                    if inUTR: # if downstream gRNA in UTR
                        gArr[1] = gArr[1] + " unlabeled (lead gRNA in UTR and this one in gene)"; # Add unlabeled explanation to status
                    else: # if downstream gRNA in gene,
                        gArr[1] = gArr[1] + " unlabeled (more than " + str(maxDistanceBetweenGRNAS) + " bp upstream of lead gRNA)"; # Add unlabeled explanation to status




        for g in allGRNAS: # for every gRNA annotated
            for gArr in gRNATable: # find this gRNA in table
                if g.seq == gArr[14]: # if the sequence is the same,
                    gArr[0] = g.label.replace(',',' '); # add this label to this row




        gRNATableNew = []; # will store new ordered list
        for g in gRNATable: # loop over list,
            if g[1][0:5] == "Valid": # if status is valid
                gRNATableNew.append(g); # add to new list


        gRNATableNew = sorted(gRNATableNew,key=operator.itemgetter(0)); # sort according to name (gRNA priority)
        for g in gRNATable: # loop over list,
            if g[1][0:14] == "Kept as backup": # if status is Kept as backup
                gRNATableNew.append(g); # add to new list


        for g in gRNATable: # loop over list,
            if g[1][0:15] == "Possible backup": # if status is Possible backup
                gRNATableNew.append(g); # add to new list


        for g in gRNATable: # loop over list,
            if g[1][0:23] == "Rejected (low on-target": # if status is Rejected due to on-target
                gRNATableNew.append(g); # add to new list


        for g in gRNATable: # loop over list,
            if g[1][0:12] == "Rejected (>9": # if status is Rejected due to >9 consecutive A/Ts
                gRNATableNew.append(g); # add to new list


        for g in gRNATable: # loop over list,
            if g[1][0:16] == "Rejected (low GC": # if status is Rejected due to gc content
                gRNATableNew.append(g); # add to new list


        gRNATableNew = [ [gene.label] + g for g in gRNATableNew ] # loop over final list adding gene to all gRNAs

        refCodon = "upstream of stop"; # string used in sgRNA table depending on end being targeted
        if not target3Prime: # if going for 5' end,
            refCodon = "downstream of start"; # string used in sgRNA table depending on end being targeted

        gRNATableString = "\n".join([",".join(g) for g in gRNATableNew]); # join array into csv string
        gRNATableString = "-,Values for rejected gRNAs, Rejected, "+enzyme+", "+str(searchStart)+" to "+str(searchEnd)+", +/-, <" + str(minGCContent) + ", <"+str(minOnTargetScore)+", "+onTargetMethod+", Not evaluated, Not evaluated, -, True, Not evaluated, Not evaluated, -, Not recoded, -\n" + gRNATableString; # Add rejected threshold
        gRNATableString = "-,Values for backup gRNAs, Backup, "+enzyme+", "+str(searchStart)+" to "+str(searchEnd)+", +/-, >=" + str(minGCContent*0.8) + ", >="+str(minOnTargetScore)+", "+onTargetMethod+", <"+str(minOffTargetScore)+", <"+str(maxOffTargetHitScore)+", "+offTargetMethod+", False, True, True, -, Recoded if " + refCodon + " codon, >=threshold\n" + gRNATableString; # Add backup threshold
        gRNATableString = "-,Values for valid gRNAs, Valid, "+enzyme+", "+str(searchStart)+" to "+str(searchEnd)+", +/-, >=" + str(minGCContent) + ", >="+str(minOnTargetScore)+", "+onTargetMethod+", >="+str(minOffTargetScore)+", >="+str(maxOffTargetHitScore)+", "+offTargetMethod+", False, False, False, -, Recoded if " + refCodon + " codon, >=threshold\n" + gRNATableString; # Add valid threshold
        gRNATableString = "Gene, Label, Status, Enzyme, Position, Strand, GC_content, On-target_score, On-target_method, Aggregated_off-target_score, Max_pairwise_off-target_score, Off-target_method, >9_consecutive_A/T, 4-Homopolymer, Triple_T, Sequence, Recoded_sequence, Recoded_sequence_pairwise_off-target_score\n" + gRNATableString; # Add column heads

        log = log + str(countBackups) + " backup gRNAs with possible off-target effects annotated.\n\n"
        if len(backupGRNAs) + len(gRNAs) == 0: # If there were absolutely no gRNAs under these settings,
            log = log + "ERROR: no gRNAs found. Please modify your criteria or select and annotate one manually.\n\n"; # say so
        else: # if there are any gRNAs,
            bestGRNA.label = bestGRNA.label + " (chosen)" # there should be a best one

    return {"out":gRNAExtreme, "log":log, "gRNATable":gRNATableString}; # returns gRNA and log


"""
Finds gRNA already annotated on gene. Filters for given restriction enzyme cut
sites. If no gRNA found on gene, nothing will happen (logs this).
"""
def findGRNA(geneGB, gene, filterCutSites=[cut_FseI,cut_AsiSI,cut_IPpoI,cut_ISceI,cut_AflII,cut_AhdI,cut_BsiWI,cut_NheI]):
    log = ""; # init log
    gRNAs = geneGB.findAnnsLabel("gRNA 1", True); # List of all gRNAs
    gRNAExtreme = GenBankAnn(); # init var to hold gRNA
    if len(gRNAs) == 0: # if no gRNAs found
        gRNAs = geneGB.findAnnsLabel("gRNA1", True); # List of all gRNAs
        if len(gRNAs) == 0: # if no gRNAs found
            gRNAs = geneGB.findAnnsLabel("gRNA"); # List of all gRNAs
            if len(gRNAs) == 0: # if no gRNAs found
                log = log + "Warning: no gRNA found on gene " + gene.label + ", selecting one automatically.\n\n"; # add warning to log
            else: # if gRNAs found,
                gRNAExtreme = gRNAs[0]; # will store gRNA most upstream
                for g in gRNAs: # loops across gRNAs
                    if ( target3Prime and g.index[0] < gRNAExtreme.index[0] ) or ( not target3Prime and g.index[0] > gRNAExtreme.index[0] ): # if more upstream
                        gRNAExtreme = g; # replace as most upstream





    gRNATableString = None
    if len(gRNAs) > 0: # if gRNAs found
        gRNAExtreme = gRNAs[0]; # will store gRNA most upstream
        gRNAExtreme = copy.deepcopy(gRNAExtreme); # fixes referencing issue. We want this to be a genuinenly new annotation
        for site in filterCutSites: # for every cut site being filtered
            if findFirst(site,gRNAExtreme.seq) > -1 or findFirst(revComp(site),gRNAExtreme.seq) > -1: # if cut site found,
                log = log + "Warning: gRNA sequence for gene " + gene.label + ": \n" + gRNAExtreme + "\ncontains restriction site " + site + "\n\n"; # add warning to log


        gRNAExtreme.label = gene.label + " gRNA (chosen)"; # renames gRNA according to this program's convention
        geneGB.features.append(gRNAExtreme);

        log = log + "gRNA for gene " + gene.label + " found on gene.\n\n"; # logs this process finished
        gRNATableString = "gRNAs not evaluated if they are user-defined.\nIf you want to check their scores, run the gene in automatic mode!\n"; # add disclaimer

    return {"out":gRNAExtreme, "log":log, "gRNATable":gRNATableString}; # returns gRNA and log


def adjustExonicRange(geneGB, gene, target3Prime, searchRange):
    exons = []
    for ann in geneGB.features: # loop through all annotations
        if ann.type == "exon" and gene.index[0] <= ann.index[0] and ann.index[1] <= gene.index[1]: # if an exon inside of this gene,
            exons.append(ann);

    exons = sorted( exons, key=lambda exon: exon.index[target3Prime], reverse=target3Prime )
        # sort exons by descending end for 3' editing, ascending start for 5'

    rangeCounter = -searchRange[0] # how much synthesis left to cover
    newRange = 0 # where we really are on gene
    curExon = exons[0] # current exon
    sectionStart = curExon.index[target3Prime]
    while rangeCounter > 0 and len(exons) > 0:
            # while still within synthesis cap and not done traversing gene
        curExon = exons[0] # current exon
        sectionStart = curExon.index[target3Prime]
        for exon in exons:
            if (1-target3Prime*2) * curExon.index[not target3Prime] > (1-target3Prime*2) * exon.index[target3Prime] and (1-target3Prime*2) * exon.index[not target3Prime] > (1-target3Prime*2) * curExon.index[not target3Prime]:
                    # if current exon overlaps with other one (for 3'/5', respectively) and the start/end of the new exon is more upstream/downstream than the current one (for 3'/5', respectively)
                curExon = exon

        change = min( abs(sectionStart - curExon.index[not target3Prime]), rangeCounter ) # take as much as synthesis allows from this section of exons
        rangeCounter -= change
        newRange = change + abs(gene.index[target3Prime] - sectionStart)

        newExons = []
        for exon in exons:
            if (1-target3Prime*2) * curExon.index[not target3Prime] < (1-target3Prime*2) * exon.index[target3Prime]:
                newExons.append(exon)

        exons = newExons

    return -newRange
