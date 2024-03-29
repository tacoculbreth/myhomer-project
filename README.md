# myhomer
myhomer is simple motif enrichment program to identify statistically significant enrichments of known motifs.

# Install Instructions
## 1. Install miniconda 
Installation requires multiple python packages. For using the packages more fluently, please install [miniconda](https://developers.google.com/earth-engine/guides/python_install-conda#mac) for your environment.
Note, `logomaker` and `scipy` required seperate installation using `pip` command. 

## 2. Install other packages
After installing miniconda, please install `requirements.txt` in the repository. 
You can follow these commands for all packages:
```
conda install --file requirements.txt
pip install logomaker
pip install scipy
```
## 3. Download data
By default, the program check motifs aginst MEME-formatted HOCOMOCOv11 database from [The MEME Suite](https://meme-suite.org/meme/doc/download.html) to GRCm38 mouse reference genome. 
- For motifs database, please download 'Motif Databases' under 'Databases' to get the latest version
- For reference genome, please download fasta and associated files from CSE 185 Public genome directory. Alternatively, you can download the mouse genome from [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc/mouse)

***Motifs database and reference genome must be formatted like this inside myhomer-project directory***
- `myhomer-project/motif_databases/`
- `myhomer-project/reference/`

# Basic Usage
For testing myhomer on your example,`cd myhomer` from the `myhomer-project` directory where necessary scripts are located. 

The basic usage is:
```
python3 myhomer.py [peaks_file] [reference_genome] [output_directory]
```

# Input
Required input to myhomer is
1. Homer Peaks file (Default: `peaks.txt`)
2. Reference fasta file (Default: GRCm38 or mm10)

# Contributors
This repository was generated by Hwa Yeon Lee, with inspiration from HOMER, SEA and CSE 185 course. 


# Testing
Run using HOMER peaks file with reference genome:
```
python3 myhomer.py ../tests/tagdirs/Klf4/peaks.txt ../reference/GRCm38.fa ../output/Klf4/
```
