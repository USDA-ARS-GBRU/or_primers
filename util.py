
import numpy as np
import primer3

def calculate_deltaG_matrix(primers):
    """
    Calculate a symmetric 2D DeltaG matrix using the calc_heterodimer
    function of primer3-py for all primers against each other.

    Parameters:
      primers (list): List of primer sequences as strings.

    Returns:
      np.ndarray: A symmetric 2D array of DeltaG values (shape n x n)
    """
    num_primers = len(primers)
    deltaG_matrix = np.zeros((num_primers, num_primers))

    for i in range(num_primers):
        for j in range(i, num_primers):
            # Run the DeltaG calculation using primer3-py calc_heterodimer.
            result = primer3.calcHeterodimer(primers[i], primers[j])
            deltaG_value = result.dg  # Assuming the function returns an object with 'dg' for DeltaG.
            
            # Update both (i, j) and (j, i) to ensure the matrix is symmetric.
            deltaG_matrix[i, j] = deltaG_value
            deltaG_matrix[j, i] = deltaG_value

    return deltaG_matrix

import random

def generate_random_dna_strings(n, l):
    """
    Generate a list containing n strings, where each string has length l.
    Each string is composed of the characters 'G', 'A', 'T', and 'C' chosen at random with replacement.

    Parameters:
      n (int): Number of strings to generate.
      l (int): Length of each string.

    Returns:
      List[str]: A list of n randomly generated DNA strings, each of length l.
    """
    letters = "GATC"
    dna_strings = [''.join(random.choices(letters, k=l)) for _ in range(n)]
    return dna_strings

# Example usage:
    n = 5  # number of strings
    l = 20  # length of each string
    result = generate_random_dna_strings(n, l)
    print("Generated DNA strings:")
    for dna in result:
        print(dna)


# Example usage:
    primers = [
        "TAGCATCACCAGACGCACAC",
        "CACCACTAGATAAGCACGGA",
        "AATTACTCGGGAAAGTTAGC",
        "TGTACCAGGCCTAGCCAGGC",
        "CAACTTGCGAGTGACCCCGC",
        "TGCTCCACTCTAATTCTTTG",
        "AGACGAACAGGTGAATTAGA",
        "TTCTTCGTTACTGAGACGTC"
    ]
    deltaG_matrix = calculate_deltaG_matrix(primers)
    print("DeltaG Matrix:")
    print(deltaG_matrix)
