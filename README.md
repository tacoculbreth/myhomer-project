# myhomer 
This is part of project for CSE 185, Spring '23. It implements motif finding against set of known motifs. 

(Work in Progress)

# Install Instructions
Installation requires the pyfaidx and pysam libraries to be installed. You can install these with pip:
```
pip install biopython pandas
```

Once required libraries are installed, you can install `myhomer` with the following command:
```
python3 setup.py install --prefix=$HOME
```

Note: if you do not have root access, you can run the commands above with additional options to install locally:
```
pip install --user biopython pandas
python3 setup.py install --user
```

If the install was successful, typing `myhomer --help` should show a useful message.

# Basic Usage
The basic usage of `myhomer` is:
```
myhomer peaks.txt ref.fasta out_dir
```

To run `myhomer` on a small test example:
```
prefix=example
myhomer example-files/peaks.txt \
    /example-files/test.fa
    /motifs/{prefix}
```

# Input
Required input to myhomer is
1. Homer Peaks file (peaks.txt)
2. Reference fasta file

Default values for size is 100. And there is option to mask. 

# File format

# Contributors
This repository was generated by Hwa Yeon Lee, with inspiration from Homer and many other projects. 
Please submit a pull request with any corrections or suggestions. 

# Testing
To run tests:
```
# do something
```