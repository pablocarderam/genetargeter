# GeneTargeter

<img src="static/assets/roundLogoNilesLab.png" width="100px" title="Niles Lab Logo">

GeneTargeter creates custom gene-editing constructs developed by the [Niles Lab at MIT](http://web.mit.edu/nileslab/), designed for knock-out or conditional knock-down in [_Plasmodium falciparum_](http://www.who.int/mediacentre/factsheets/fs094/en/),  by delivering a 3' or 5' UTR post-transcriptional regulatory element payload to a specific given gene.

Most applications can be served well through the web application running at:

[genetargeter.mit.edu](genetargeter.mit.edu)

NOTE: the app sometimes takes a good 30 seconds to load if nobody's used it in a while.

However, if you really want to run it locally, you can do so as follows:

## Installation

Clone the directory and run

```bash
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

## Authors
I'm Pablo Cárdenas, a member of the [Niles Lab](http://web.mit.edu/nileslab/) at [MIT Biological Engineering](be.mit.edu). Follow my science antics at [@pcr_guy on Twitter!](https://twitter.com/pcr_guy)

This project has been possible thanks to the invaluable help of many others, particularly Jacquin Niles, Sumanta Dey, and Lisl Escherick.

Cheers and thanks for using!

## License
Well [MIT](https://choosealicense.com/licenses/mit/), of course!
