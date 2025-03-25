# Mixed Linear Integer optimizers of primer selection

Adam R. Rivers 
March 25, 2025
## Problem

Designing multiplexed pcr reactions is challenging because as the number of primers grows the chance of primer dimerization grows exponentially.  Primer pairs with a Gibbs Free energy (delta G) of less than -9000 cal/mol are considered to be at risk for dimerizing.  Initially, we tried to maximize the global delta G of the selected primers. As we thought about it more we realized that this mattered less than minimizing the number of pairs with a deltaG below the threshold deltaG.

This repository contains some example optimizations for doing this with synthetic data or user supplied data

## cp-sat_minimize_num_low_deltaG.py:  

This script uses a CP-SAT model (from OR-Tools) to select one candidate primer pair per block
and minimizes the number of blocks that have a primer with a DeltaG < -9000.0 cal/mol with any
selected primer from a different block.

It can run in two modes:
  1. Generation mode: Provide --kmers, --num_blocks, and --candidates_per_block.
     A  deltaG matrix and candidate blocks will be generated using the first 
     (num_blocks * candidates_per_block) kmer fasta records. Kmers should be primer length ~20nt. 
     I used bbtools kmercountexact.sh to generate my kmer fasta. 
  2. File input mode: Provide --deltaG_csv and --block_json to load the deltaG matrix and candidate
     block definitions from files.

After solving, the script prints the selected candidate pairs for each block, and computes:
    - The mean DeltaG for primers within each selected pair (intra-block DeltaG).
    - The mean DeltaG for all interactions between selected primers from different blocks (inter-block DeltaG).

Usage examples:
  Generation mode:
    `cp-sat_minimize_num_low_deltaG_cli.py --kmers data/foot-and-mouth-dir/GCF_002816555.k21.fa --num_blocks 20 --candidates_per_block 5`

  File input mode:
    `cp-sat_minimize_num_low_deltaG_cli.py --deltaG_csv deltaG_matrix.csv --block_json block_list.json`

# Data 
    * `fmv_blocks_from_Ao`: 
        1. `deltaG_matrix.csv`: a distance matrix of primers  selected by primer3 in spaced blocks.
        2. `block_list.json`: a list of pairs selected.
    * `foot-and-mouth-dir`:
        1. `foot-and-mouth.sh`: a script to  create kemrs from a genome
        2. `GCF_002816555.1_ASM281655v1_genomic.fna`: genome
        3. `GCF_002816555.k21.fa`: fasta of k=21 mers from `bbtools kmercountexact.sh`

# Other files
    * `util.py` some utility functions for testing
    * `prompts` text used to generate the initial ortools code from Qodo gen