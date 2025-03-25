import sys
import os
from typing import List, Tuple
from itertools import combinations, islice

import numpy as np

from Bio import SeqIO
import primer3


#!/usr/bin/env python3
"""
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
    cp-sat_minimize_num_low_deltaG_cli.py --kmers data/foot-and-mouth-dir/GCF_002816555.k21.fa --num_blocks 20 --candidates_per_block 5

  File input mode:
    cp-sat_minimize_num_low_deltaG_cli.py --deltaG_csv data/fmv_blocks_from_Ao/deltaG_matrix.csv --block_json data/fmv_blocks_from_Ao/block_list.json
"""

import argparse
import json
import numpy as np
import random
import sys
from ortools.sat.python import cp_model
from Bio import SeqIO
from itertools import islice

def generate_random_candidates(num_blocks: int, candidates_per_block: int) -> list:
    """
    Generate random candidate primer pairs for a specified number of blocks.

    Each block contains a specified number of candidate pairs, where each pair
    is a tuple of two distinct primer indices randomly selected from the total
    number of candidates.

    Parameters:
        num_blocks (int): The number of blocks to generate candidates for.
        candidates_per_block (int): The number of candidate pairs per block.

    Returns:
        list: A list of lists, where each inner list contains tuples representing
        candidate primer pairs for a block.

    Raises:
        ValueError: If num_blocks or candidates_per_block are not positive integers.
    """
    if not isinstance(num_blocks, int) or not isinstance(candidates_per_block, int):
        raise ValueError("Both num_blocks and candidates_per_block must be integers.")
    if num_blocks <= 0 or candidates_per_block <= 0:
        raise ValueError("Both num_blocks and candidates_per_block must be positive integers.")
    
    return [
        [
            tuple(np.random.choice(num_blocks * candidates_per_block, size=2, replace=False))
            for _ in range(candidates_per_block)
        ]
        for _ in range(num_blocks)
    ]

DEFAULT_THRESHOLD = -9000.0

def compute_inter_bad_indicator(candidate1: List[int], candidate2: List[int], deltaG_matrix: np.ndarray, threshold: float = DEFAULT_THRESHOLD) -> int:
    """
    Compute an indicator for a bad interaction between two candidate pairs.
    
    Parameters:
    - `candidate1`: List of primer indices representing the first candidate pair.
    - `candidate2`: List of primer indices representing the second candidate pair.
    - `deltaG_matrix`: 2D numpy array containing deltaG values for primer interactions.
    - `threshold`: Float value representing the deltaG threshold for bad interactions.
    
    There are 4 interactions (each primer from candidate1 with each primer from candidate2).
    Return 1 if any of these deltaG values is below the threshold; otherwise 0.
    """
    if not isinstance(deltaG_matrix, np.ndarray):
        raise TypeError('deltaG_matrix must be a NumPy array')
    for p in candidate1:
        for q in candidate2:
            if deltaG_matrix[p, q] < threshold:
                return 1
    return 0

def build_and_solve_model(blocks: list, deltaG_matrix: np.array, deltag_cutoff: float) ->Tuple[list[Tuple[np.int64]], int]:
    """
        Build and solve a CP-SAT model to select one candidate primer pair per block,
        minimizing the number of blocks with a primer having a deltaG below the specified cutoff
        when interacting with primers from all other blocks.

        Parameters:
        - blocks: List of blocks, each containing candidate primer pairs.
        - deltaG_matrix: 2D numpy array of deltaG values for primer interactions.
        - deltag_cutoff: Float value representing the deltaG threshold for bad interactions.

        Returns:
        - Tuple containing:
        - selected_candidates: List of tuples, each with a block index and the selected primer pair.
        - num_faults: Integer count of blocks with flagged bad inter-block interactions.
        """
    num_blocks = len(blocks)
    candidates_per_block = len(blocks[0])
    
    model = cp_model.CpModel()
    
    # Decision variables: x[b,j] == 1 if candidate j in block b is selected.
    x = {}
    for b in range(num_blocks):
        for j in range(candidates_per_block):
            x[(b, j)] = model.NewBoolVar(f"x_{b}_{j}")
    
    # Each block must have exactly one candidate selected.
    for b in range(num_blocks):
        model.Add(sum(x[(b, j)] for j in range(candidates_per_block)) == 1)
    
    # Fault variable for each block: fault[b] = 1 if the selected candidate in block b
    # interacts badly (deltaG < -9.0) with a selected candidate in any other block.
    fault = {}
    for b in range(num_blocks):
        fault[b] = model.NewBoolVar(f"fault_{b}")
    
    # Auxiliary variables for bad interactions among blocks.
    y = {}
    for b1 in range(num_blocks):
        for b2 in range(b1 + 1, num_blocks):
            for j in range(candidates_per_block):
                for k in range(candidates_per_block):
                    if compute_inter_bad_indicator(blocks[b1][j], blocks[b2][k], deltaG_matrix, threshold=deltag_cutoff):
                        y[(b1, b2, j, k)] = model.NewBoolVar(f"y_{b1}_{b2}_{j}_{k}")
                        # Link y with the selected candidate variables.
                        model.Add(x[(b1, j)] + x[(b2, k)] - 1 <= y[(b1, b2, j, k)])
                        model.Add(y[(b1, b2, j, k)] <= x[(b1, j)])
                        model.Add(y[(b1, b2, j, k)] <= x[(b2, k)])
                        # Flag both blocks if a bad interaction is found.
                        model.Add(fault[b1] >= y[(b1, b2, j, k)])
                        model.Add(fault[b2] >= y[(b1, b2, j, k)])
    
    # Objective: minimize the total number of blocks flagged.
    model.Minimize(sum(fault[b] for b in range(num_blocks)))
    
    # Solve model.
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = os.cpu_count()  # Dynamically set based on available CPU cores
    solver.parameters.log_search_progress = True
    solver.parameters.max_time_in_seconds = 300.0
    status = solver.Solve(model)
    
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        selected_candidates = []
        for b in range(num_blocks):
            for j in range(candidates_per_block):
                if solver.Value(x[(b, j)]) == 1:
                    selected_candidates.append((b, blocks[b][j]))
        blocks_with_faults = sum(solver.Value(fault[b]) for b in range(num_blocks))
        print("Solution found!")
        print(f"Total blocks with a bad inter-block interaction: {blocks_with_faults}")
        print("Selected candidate pairs per block (block, (primer1, primer2)):")
        for sol in selected_candidates:
            print(sol)
        return selected_candidates, blocks_with_faults
    else:
        print("No solution found.")
        sys.exit(1)

def compute_mean_deltaG(selected_candidates: List[Tuple[int, Tuple[int, int]]], deltaG_matrix: np.ndarray) -> Tuple[float, float]:
    """
    Computes the mean DeltaG values for selected primer candidates.

    This function calculates two types of mean DeltaG:
    1. Intra-block DeltaG: The average DeltaG between the two primers in each selected candidate pair.
    2. Inter-block DeltaG: The average DeltaG for all interactions between primers from different blocks,
       where each pair of selected candidates contributes four interactions.

    Args:
        selected_candidates (List[Tuple[int, Tuple[int, int]]]): A list of tuples where each tuple contains
            an integer block index and a tuple of two integers representing selected primer indices.
        deltaG_matrix (np.ndarray): A 2D numpy array representing the DeltaG values between primers.

    Returns:
        Tuple[float, float]: A tuple containing the mean intra-block DeltaG and the mean inter-block DeltaG.

    Raises:
        ValueError: If deltaG_matrix is not a 2D numpy array or if selected_candidates is not a list of tuples.
    """
    if not isinstance(deltaG_matrix, np.ndarray) or deltaG_matrix.ndim != 2:
        raise ValueError("deltaG_matrix must be a 2D numpy array.")
    if not all(isinstance(candidate, tuple) and len(candidate) == 2 for candidate in selected_candidates):
        raise ValueError("selected_candidates must be a list of tuples.")
    if not selected_candidates:
        return float('nan'), float('nan')

    intra_values = [deltaG_matrix[p, q] for _, (p, q) in selected_candidates]
    intra_mean = np.mean(intra_values)

    inter_values = [deltaG_matrix[p1, q1] for (_, (p1, q1)), (_, (p2, q2)) in combinations(selected_candidates, 2) for p1 in (p1, p2) for q1 in (q1, q2)]
    inter_mean = np.mean(inter_values)

    return intra_mean, inter_mean

def load_deltaG_matrix(csv_path: str) -> np.ndarray:
    """
    Loads a symmetric deltaG matrix from a CSV file.

    Parameters:
        csv_path (str): The path to the CSV file containing the deltaG matrix.

    Returns:
        np.ndarray: The loaded symmetric deltaG matrix.

    Raises:
        ValueError: If the loaded matrix is not symmetric.
        OSError: If there is an error reading the file.
    """
    try:
        # Assuming the csv file is formatted with comma separation.
        matrix = np.loadtxt(csv_path, delimiter=',')
        if not np.allclose(matrix, matrix.T):
            raise ValueError("The loaded matrix is not symmetric. Shape is: {}".format(matrix.shape))
        return matrix
    except OSError as e:
        print(f"File error reading {csv_path}: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Data error in {csv_path}: {e}")
        sys.exit(1)

def load_blocks(json_path: str) -> dict:
    """
    Loads block candidate definitions from a JSON file specified by the given path.

    Args:
        json_path (str): The file path to the JSON file containing block definitions.

    Returns:
        dict: A dictionary containing the block candidate definitions.

    Raises:
        SystemExit: If the file is not found or if there is an error decoding the JSON.
    """
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            blocks = json.load(f)
        return blocks
    except FileNotFoundError:
        print(f"File not found: {json_path}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from {json_path}: {e}")
        sys.exit(1)

def calculate_deltaG_matrix(primers: List[str]) -> np.ndarray:
    """
    Calculate a symmetric 2D DeltaG matrix for a list of primers.

    This function computes the DeltaG values for all possible heterodimer
    combinations of the given primers using the primer3 library. The resulting
    matrix is symmetric, with each element representing the DeltaG value for
    the interaction between two primers.

    Parameters:
        primers (List[str]): A list of primer sequences as strings.

    Returns:
        np.ndarray: A symmetric 2D array of DeltaG values with shape (n x n),
        where n is the number of primers.

    Raises:
        ValueError: If any primer is not a non-empty string containing only
        the nucleotides A, T, C, and G.
    """
    # Validate primers list
    valid_nucleotides = {'A', 'T', 'C', 'G'}
    for primer in primers:
        if not primer or not set(primer).issubset(valid_nucleotides):
            raise ValueError("Primers must be non-empty strings containing only A, T, C, G.")

    num_primers = len(primers)
    deltaG_matrix = np.zeros((num_primers, num_primers))
    indices = np.triu_indices(num_primers)

    for i, j in zip(*indices):
        try:
            result = primer3.calc_heterodimer(primers[i], primers[j])
            deltaG_value = result.dg
        except Exception as e:
            print(f"Error calculating DeltaG for primers {i} and {j}: {e}")
            deltaG_value = float('inf')  # or another appropriate default value
        deltaG_matrix[i, j] = deltaG_value

    deltaG_matrix += np.triu(deltaG_matrix, 1).T
    return deltaG_matrix

def data_from_kmers(kmers: str, num_blocks: int, candidates_per_block: int) -> Tuple[np.ndarray, List[List[Tuple[int, int]]]]:
    """
    Generate a DeltaG matrix and random candidate primer pairs from a k-mers file.

    Parameters:
        kmers (str): Path to a file containing k-mer sequences in uncompressed FASTA format.
        num_blocks (int): Number of blocks to divide the sequences into.
        candidates_per_block (int): Number of candidate primer pairs per block.

    Returns:
        Tuple[np.ndarray, List[List[Tuple[int, int]]]]:
            - A symmetric 2D DeltaG matrix of the sequences.
            - A list of lists containing tuples of random candidate primer pairs for each block.

    Raises:
        ValueError: If `kmers` is not a valid file path or file-like object.
    """
    if not isinstance(kmers, (str, bytes)) or not os.path.exists(kmers):
        raise ValueError("`kmers` must be a valid file path or file-like object.")
    
    try:
        seqs = SeqIO.parse(kmers, 'fasta')
        if not any(seqs):
            raise ValueError("The kmers file does not contain valid FASTA sequences.")
    except Exception as e:
        print(f"Error parsing kmers file: {e}")
        return None, None
    
    totseqs = num_blocks * candidates_per_block * 2
    seqlist = [str(rec.seq) for rec in islice(seqs, totseqs)]
    deltag_matrix = calculate_deltaG_matrix(seqlist)

    block_list = generate_random_candidates(num_blocks, candidates_per_block)
    return deltag_matrix, block_list

def parse_arguments():
    """
    Parse command-line arguments for the CP-SAT model script.

    This function sets up an argument parser to handle two modes of operation:
    1. Generation mode: Requires --num_blocks, --kmers, and --candidates_per_block
    to generate candidate primer pairs.
    2. File input mode: Requires --deltaG_csv and --block_json to load existing
    deltaG matrix and block definitions.

    Additional parameters include:
    - --deltag_cutoff: The cutoff DeltaG value in cal/mol (default: -9000.0).

    The function also validates the existence and format of input files and ensures
    that required arguments are provided for each mode.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="CP-SAT model to minimize bad deltaG interactions between selected primer pairs."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    # Mode 1: Generation using parameters.
    group.add_argument(
        "--num_blocks", type=int, help="Number of blocks for candidate generation."
    )
    # Mode 2: File input (deltaG matrix and block definitions).
    group.add_argument(
        "--deltaG_csv", type=str,
        help="Path to a CSV file with the deltaG matrix (file input mode)."
    )

    # Additional parameters: if generation mode is used.
    parser.add_argument("--kmers", type=str, help="a fasta file of kmers around k=21 for simulation")
    parser.add_argument("--candidates_per_block", type=int, help="Candidates per block for generation.")
    # For file mode, block JSON must be provided.
    parser.add_argument("--block_json", type=str, help="Path to a JSON file with block candidate definitions.")
    parser.add_argument("--deltag_cutoff", type=float, default=-9000.0, help="The cutoff G value in cal/mol")

    args = parser.parse_args()

    if args.kmers:
        try:
            with open(args.kmers, 'r') as file:
                # Attempt to parse the file to ensure it's in the correct format
                SeqIO.parse(file, "fasta")
        except Exception as e:
            print(f"Error: The kmers file is not in the correct format or is invalid. Details: {e}")
            sys.exit(1)

    if args.deltaG_csv and not os.path.isfile(args.deltaG_csv):
        print(f"Error: The file {args.deltaG_csv} does not exist or is not accessible.")
        sys.exit(1)

    if args.block_json and not os.path.isfile(args.block_json):
        print(f"Error: The file {args.block_json} does not exist or is not accessible.")
        sys.exit(1)

    if args.num_blocks is not None or args.kmers is not None or args.candidates_per_block is not None:
        if not (args.num_blocks and args.kmers and args.candidates_per_block):
            parser.error("--num_blocks, --kmers, and --candidates_per_block must be provided together in generation mode.")
    
    if args.deltaG_csv is not None:
        if not args.block_json:
            parser.error("--block_json must be provided when using file input mode.")

    return args

if __name__ == "__main__":
    args = parse_arguments()

    # Mode selection based on provided arguments.
    if args.deltaG_csv:
        # File mode - both deltaG matrix and block JSON paths are required.
        if not args.block_json:
            print("Error: --block_json must be provided when using file input mode.")
            sys.exit(1)
        deltaG_matrix = load_deltaG_matrix(args.deltaG_csv)
        blocks = load_blocks(args.block_json)
        print("Using file input mode.")
    else:
        # Generation mode.
        # Check required parameters.
        if args.num_blocks is None or args.candidates_per_block is None or args.kmers is None:
            print("Error: --num_blocks, --kmers and --candidates_per_block must be provided in generation mode.")
            sys.exit(1)
        deltaG_matrix, blocks = data_from_kmers(kmers = args.kmers, num_blocks=args.num_blocks, candidates_per_block=args.candidates_per_block)
        print("Using generation mode with random data.")
    
    # For reproducibility.
    np.random.seed(42)
    random.seed(42)
    
    # Build and solve the model.
    selected_candidates, faults = build_and_solve_model(blocks, deltaG_matrix, deltag_cutoff=args.deltag_cutoff)
    
    # Compute mean deltaG statistics.
    intra_mean, inter_mean = compute_mean_deltaG(selected_candidates, deltaG_matrix)
    print("\nMean intra-block DeltaG (within selected candidate pairs):", intra_mean)
    print("Mean inter-block DeltaG (across selected primers from different blocks):", inter_mean)
    
    # The script returns the selected candidate pairs and the computed means.
