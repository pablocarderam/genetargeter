#TEMP
#gRNA = "gttagaaatggacgactaat"; #TODO: choose gRNA
#LHR = "tcacagatgatgaaattgttcatgctctaaaattgaatggtataaatttaggtaaaaacgatttatataaatatatgcataaacaagatatgaaatcaaattataaaaaaataatacaaacatcaaaagtaataaaccaatcaaatgataatactattctattaacaaatgattgtataagatatctatctttagttagactttatttaaatcgacataaatataaaatcatattaatagatgaaattcctatttttaatttaaacaattctgttcatgacgaattaaatagttttttaattggtaaagcaaagtcatttaattatataataagaaatcattttccaaataatacagtcctaattatttcacatcatgcaaatactttgtcttgttgtgactatatttatgtattaagaaagggagaaataacttatcgttgtagttacgaagatgtaaaaacgcaatctgaattatcacatttgttagaaatggacgact"; #TODO: choose LHR
#RHR = "cggtatattattacatgaataaacttacacacacccatgcatacataatacatacatacatacatacatacatacatacatatatatatatatatatatatatataattttaaccgtttataaagaaaaaattgctttaaaagaaagattacaattttatatatttacctcagaattaatctatagatcatatatattttattttgtttttctgtataaatcacaattaaaagaagtgccttataataaaacattgtgaactataaggtagtttttttttttttttttctttattttaattcccattgatatgctccctaagcatgttaaaaaaagatgcaatataatatgtaatatttttattcttttaagaataatgaaaaatttcttcacctttttttaagtgtaatgaatataatatatatatatatatatatatatatatgttgatcatattaataaattcctaattgaaaatgtgtaaacaaaaaaaaaaaaaaaaattatttgtgtttgtaccgaaagtaataagatataacttgagttaaattaaaaaataaaataaagatataaccatcataatatgcataataatatatatatatatatatatatcaattttatattaaaaaatataattaaagattatataaaatatataaaaaacacttatgcttacctgcttaagtacaaataaaacacgtatatg"; #TODO: choose RHR

from py.genetargeter.GeneTargeterMethods import *;
mrp1 = GenBank();
mrp1.load("../Input/pf3d7_0112200-mrp1.gb");
rhr = chooseRHR(mrp1, "MRP1");
lhr = chooseLHR(mrp1, "MRP1");
mrp1.features.append(lhr);
mrp1.features.append(rhr);
mrp1.save("../output/MRP1_test.gb"); # saves to file

pSN054 = GenBank();
pSN054.load("../Input/psn054-updated_v5-tagged-jn_final.gb"); # load plasmid sequence from GenBank format
pSN054_ARMED = insertTargetingElementsPSN054(pSN054, gRNA, LHR, RHR, "MRP1"); # inserts
pSN054_ARMED.save("../output/MRP1_test.gb"); # saves to file

from py.genetargeter.GeneTargeterMethods import *;
p = pSN054TargetGene("PF3D7_0302600 (ABC Transporter)", "pf3d7_0302600-abcb4.gb");
#p = pSN054TargetGene("PF3D7_0319700 (ABCI3)", "pf3d7_0319700-abci3.gb");
p = pSN054TargetGene("PF3D7_0523000 (MDR1)", "pf3d7_0523000-mdr1.gb");
p = pSN054TargetGene("PF3D7_0810200 (ABCK1)", "pf3d7_0810200-abck1.gb");
p = pSN054TargetGene("PF3D7_0813700 (ABCF1)", "pf3d7_0813700-abcf1.gb");
p = pSN054TargetGene("PF3D7_1145500 (ABCB3)", "pf3d7_1145500-abcb3.gb");
p = pSN054TargetGene("PF3D7_1145500 (ABCB3)", "pf3d7_1145500_abcb3.gb");
p = pSN054TargetGene("PF3D7_0106900 (IspD)", "pf3d7_0106900-ispd.gb", HRannotated=True);

from py.genetargeter.GeneTargeterMethods import *;
txt = open("input/genes/pf3d7_0106900-ispd.gb"); # Access given file
d = txt.read(); # Read file
p = pSN054TargetGene("PF3D7_0106900 (IspD)", d, HRannotated=True, useFileStrs=True);

pSN054BuildfFromFile()
