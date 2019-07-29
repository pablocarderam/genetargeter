# GeneTargeter

<img src="static/assets/roundLogoNilesLab.png" width="400px" title="Niles Lab Logo">

GeneTargeter creates custom gene-editing constructs developed by the [Niles Lab at MIT](http://web.mit.edu/nileslab/), designed for knock-out or conditional knock-down in [_Plasmodium falciparum_](http://www.who.int/mediacentre/factsheets/fs094/en/),  by delivering a 3' or 5' UTR post-transcriptional regulatory element payload to a specific given gene.

Most applications can be served well through the web application running at [genetargeter.mit.edu](genetargeter.mit.edu)

NOTE: the app sometimes takes a good 30 seconds to load if nobody's used it in a while.

However, if you really want to run it locally, you can do so as follows:

## Installation

Clone the repo and run setup from within the main directory:

```bash
git clone git://github.com/schacon/grit.git
cd GeneTargeter
python setup.py install
```

## Usage

To run GeneTargeter as a command line application, use the following command:

```bash
genetargeter PATH/TO/FILE.gb ParameterFile.txt PATH/TO/OUTPUTFOLDER NAME_OF_GENE
```

The last parameter is optional; if empty, the name will be taken from the name of the file. Additionally, instead of the name of a single file, a directory may be listed. In this case, all files in the folder will be processed and each gene's name
will be taken from the name of its file.

## Contributing
Pull requests are welcome, but we'd really appreciate your feedback about it! Please contact [pablocarderam at gmail dot com](pablocarderam@gmail.com) and open an issue first to discuss what you would like to change.

## Acknowledgements
This project has been possible thanks to the invaluable help of many others, particularly Drs. Jacquin Niles, Sumanta Dey, and Lisl Escherick.

On-target Cas9 gRNA scoring code and CFD scoring matrix obtained December 27, 2016 from [Doench et al. (2016)](https://dx.doi.org/10.1038/nbt.3437) supplementary material.

On-target Cas12 gRNA scoring coefficients obtained from [Kim et al. (2017)](https://dx.doi.org/10.1038/nmeth.4104) supplementary material, equivalent to the [CINDEL](http://big.hanyang.ac.kr/cindel/) online tool. Exact on-target scores may vary between GeneTargeter and [CINDEL](http://big.hanyang.ac.kr/cindel/) due to minor differenes in gRNA free energy calculations.

Self-folding Gibbs free energy of gRNAs (used to calculate [Kim et al. (2017)](https://dx.doi.org/10.1038/nmeth.4104) [CINDEL](http://big.hanyang.ac.kr/cindel/) scores) calculated using RNAfold obtained from [ViennaRNA v2.3.3](http://www.tbi.univie.ac.at/RNA/index.html) [(Lorenz et al., 2011)](https://dx.doi.org/10.1186/1748-7188-6-26).

Off-target Zhang Lab scoring matrix obtained December 27, 2016 from the online [MIT CRISPR tool](http://crispr.mit.edu/about) [(Hsu et al., 2016)](https://dx.doi.org/10.1038/nbt.2647).

Off-target gRNA database built from the _P. falciparum_ genome first published by [ Gardner et al. (2002)](https://dx.doi.org/10.1038/nature01097) and subsequently edited by the [PlasmoDB](http://plasmodb.org/plasmo/) community.

Codon frequency tables obtained from the [High-performance Integrated Virtual Environment-Codon Usage Tables (HIVE-CUT) database](https://hive.biochemistry.gwu.edu/cuts/about) [(Athey et al., 2017)](https://dx.doi.org/10.1186/s12859-017-1793-7).


## Authors

I'm Pablo Cárdenas, a member of the [Niles Lab](http://web.mit.edu/nileslab/) at [MIT Biological Engineering](be.mit.edu). Follow my science antics at [@pcr_guy on Twitter!](https://twitter.com/pcr_guy)

Cheers and thanks for using!

## License
Well [MIT](https://choosealicense.com/licenses/mit/), of course!
