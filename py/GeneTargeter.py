"""
Calls main method to build targeted vector with arguments from the command line.
TODO: GUI
"""

from sys import argv; # Needed for receiving user input
from GeneTargeterMethods import *;

script, geneName, geneFile = argv[0:3]; # From console, argv returns script name, arguments
HRann = False; # default HRannotated value
if len(argv) == 4: # if only the script and three arguments passed,
    HRann = True; # signal HRs are annotated

pSN054TargetGene(geneName, geneFile, HRannotated=HRann); # call result

print geneName + " targeted. Results in output folder."

# Usage:
# python py/GeneTargeter.py "PF3D7_0106900 (IspD)" "pf3d7_0106900-ispd.gb"
# python py/GeneTargeter.py "PF3D7_0106900 (IspD)" "pf3d7_0106900-ispd.gb" True
