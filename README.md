# MAHDS_addition
This repository contains some additional software for testing and assessing MAHDS.

## Contents
(It is assumed that all paths begin with the root of the repository)

### fmax
Location: `fmax/`.  
Installation: `cd fmax && make`.  
Input (command line parameters): path to alignment in FASTA format (required), paths to configuration files (optional). *Note: by default configuration files are located in `fmax/install/cfgdir`*.  
Output (stdout): Fmax value.  

### artif_seq_gen
Location: `testing_utils/artificial_seq_tests/artif_seq_gen.py`.  
Installation: not required.  
Input (command line parameters): path to configuration file (optional). *Note: by default configuration file is `testing_utils/artificial_seq_tests/asg_config.json`*.  
Output (files): files in FASTA format containing sequences with given parameters.  

### EMBL-EBI_test
Location: `testing_utils/EMBL-EBI_methods/EMBL-EBI_test.py`.  
Installation: not required.  
Input (command line parameters): path to configuration file (optional). *Note: by default configuration file is `testing_utils/EMBL-EBI_methods/eet_config.json`*.  
Output (files): files in FASTA format containing multiple alignments built by services provided by EMBL-EBI.  

### rnd_seq_gen
Location: `testing_utils/rnd_sequences/rnd_seq_gen.py`.  
Installation: not required.  
Input (command line parameters): path to configuration file (optional). *Note: by default configuration file is `testing_utils/rnd_sequences/rsg_config.json`*.  
Output (files): files in FASTA format containing sequences with given parameters.  

### get_di_stat
Location: `testing_utils/stat_del-ins/get_di_stat.py`.  
Installation: not required.  
Input (command line parameters): path to alignment in FASTA format (required).  
Output (stdout): number of gaps and gap openings in given alignment.  
