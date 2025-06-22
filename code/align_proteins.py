#!/usr/bin/env python

# To install Biopython: pip install biopython
from Bio import Align
import Bio.Align

def calculate_protein_similarity(seq1: str, seq2: str) -> float:
    """
    Calculates the percentage identity between two protein sequences by direct positional comparison.

    This function assumes the sequences are of the same length.
    For sequences of different lengths or to account for insertions/deletions/substitutions,
    more advanced algorithms (e.g., Needleman-Wunsch or Smith-Waterman) are required.

    Args:
        seq1: The first protein sequence (string of amino acid single-letter codes).
        seq2: The second protein sequence (string of amino acid single-letter codes).

    Returns:
        A float representing the percentage similarity (identity) between the two sequences.
        Returns -1.0 if sequences are of different lengths to indicate an error or
        invalid comparison.
    """
    # Convert sequences to uppercase to ensure case-insensitivity
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # Check if sequences have the same length
    if len(seq1) != len(seq2):
        print("Error (Simple Identity): Sequences must be of the same length for direct identity comparison.")
        return -1.0  # Indicate an error

    matches = 0
    total_length = len(seq1)

    # Iterate through the sequences and count identical amino acids
    for i in range(total_length):
        if seq1[i] == seq2[i]:
            matches += 1

    # Calculate percentage similarity
    if total_length > 0:
        similarity_percentage = (matches / total_length) * 100
    else:
        similarity_percentage = 0.0 # Handle empty sequences

    return similarity_percentage

def calculate_protein_similarity_biopython(seq1: str, seq2: str) -> float:
    """
    Calculates the percentage identity between two protein sequences using Biopython's
    pairwise aligner (Needleman-Wunsch algorithm for global alignment).

    This method is more robust for biological sequences as it accounts for gaps
    (insertions/deletions) and calculates identity based on the best alignment.

    Args:
        seq1: The first protein sequence (string of amino acid single-letter codes).
        seq2: The second protein sequence (string of amino acid single-letter codes).

    Returns:
        A float representing the percentage identity of the best global alignment.
        Returns 0.0 if an alignment cannot be found or sequences are empty.
    """
    aligner = Bio.Align.PairwiseAligner()
    # Set the substitution matrix for protein alignment (e.g., BLOSUM62)
    # If not specified, default is usually identity matrix or a simple PAM/BLOSUM.
    # For more accurate biological similarity, you'd load a specific matrix.
    # For *identity*, we usually set match score to 1 and mismatch to 0.
    # However, Biopython's aligner calculates identity based on the alignment it finds.
    # For simple identity, we can use an identity matrix or customize scores.
    # For percentage identity, a simple approach is to count matches in the alignment.

    # A common way to get identity is to count identical residues in the alignment.
    # We can use biopython's default scoring or set it up explicitly for identity.
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0 # Mismatches get 0 score
    aligner.open_gap_score = -0.5 # Penalty for opening a gap
    aligner.extend_gap_score = -0.1 # Penalty for extending a gap

    # Perform global alignment
    alignments = aligner.align(seq1.upper(), seq2.upper())

    if not alignments:
        return 0.0 # No alignment found

    # Get the best alignment (Biopython returns an iterator of alignments)
    best_alignment = alignments[0]

    # The alignment object has attributes like `identity`, `score` etc.
    # For `identity`, we typically look at the number of identical characters
    # divided by the length of the shorter sequence or the aligned length.
    # Biopython's `format_alignment` gives a visual. We need to calculate identity.

    # Calculate identity based on the aligned sequences.
    # Iterate through the aligned sequences (represented as strings with '-')
    aligned_seq1 = str(best_alignment).split('\n')[0]
    aligned_seq2 = str(best_alignment).split('\n')[2] # The third line is the second sequence

    # Ensure sequences are valid and not empty
    if not aligned_seq1 or not aligned_seq2:
        return 0.0

    matches = 0
    # The length of the alignment can be different from original sequence lengths due to gaps.
    # We calculate identity based on the length of the alignment.
    alignment_length = len(aligned_seq1) # Both aligned sequences will have the same length

    for i in range(alignment_length):
        if aligned_seq1[i] == aligned_seq2[i] and aligned_seq1[i] != '-':
            matches += 1

    if alignment_length > 0:
        identity_percentage = (matches / alignment_length) * 100
    else:
        identity_percentage = 0.0

    return identity_percentage


# --- Example Usage ---
if __name__ == "__main__":
    print("--- Simple Positional Identity Calculation ---")
    # Example 1: Identical sequences
    protein_seq_a = "MVLSPADKTNVKAAWGKVGAHAG"
    protein_seq_b = "MVLSPADKTNVKAAWGKVGAHAG"
    similarity = calculate_protein_similarity(protein_seq_a, protein_seq_b)
    if similarity != -1.0:
        print(f"Similarity between '{protein_seq_a}' and '{protein_seq_b}': {similarity:.2f}%")
    print("-" * 30)

    # Example 2: Partially similar sequences
    protein_seq_c = "MVLSPADKTNVKAAWGKVGAHAG"
    protein_seq_d = "MVLSPEAKTNVKAAYGKVGAHAS" # Some differences
    similarity = calculate_protein_similarity(protein_seq_c, protein_seq_d)
    if similarity != -1.0:
        print(f"Similarity between '{protein_seq_c}' and '{protein_seq_d}': {similarity:.2f}%")
    print("-" * 30)

    # Example 3: Completely different sequences
    protein_seq_e = "ABCDEFGHIJKLM"
    protein_seq_f = "ZYXWVUTSRQPON"
    similarity = calculate_protein_similarity(protein_seq_e, protein_seq_f)
    if similarity != -1.0:
        print(f"Similarity between '{protein_seq_e}' and '{protein_seq_f}': {similarity:.2f}%")
    print("-" * 30)

    # Example 4: Sequences of different lengths (will result in an error message)
    protein_seq_g = "MVLSPAD"
    protein_seq_h = "MVLSPADKTNVKAAWGKVGAHAG"
    similarity = calculate_protein_similarity(protein_seq_g, protein_seq_h)
    if similarity == -1.0:
        print("Comparison failed due to unequal sequence lengths.")
    print("-" * 30)

    # Example 5: Empty sequences
    protein_seq_i = ""
    protein_seq_j = ""
    similarity = calculate_protein_similarity(protein_seq_i, protein_seq_j)
    if similarity != -1.0:
        print(f"Similarity between '{protein_seq_i}' and '{protein_seq_j}': {similarity:.2f}%")
    print("=" * 50)


    print("\n--- Biopython Pairwise Alignment Identity Calculation ---")
    # Example 1: Identical sequences (Biopython)
    protein_seq_a_bio = "MVLSPADKTNVKAAWGKVGAHAG"
    protein_seq_b_bio = "MVLSPADKTNVKAAWGKVGAHAG"
    similarity_bio = calculate_protein_similarity_biopython(protein_seq_a_bio, protein_seq_b_bio)
    print(f"Biopython Similarity between '{protein_seq_a_bio}' and '{protein_seq_b_bio}': {similarity_bio:.2f}%")
    print("-" * 30)

    # Example 2: Partially similar sequences with gaps (Biopython)
    protein_seq_c_bio = "MVLSPADKTNVKAAWGKVGAHAG"
    protein_seq_d_bio = "MVLSPEAKTNVKAAYGKVGAHAS-G" # Differences and a gap
    similarity_bio = calculate_protein_similarity_biopython(protein_seq_c_bio, protein_seq_d_bio)
    print(f"Biopython Similarity between '{protein_seq_c_bio}' and '{protein_seq_d_bio}': {similarity_bio:.2f}%")
    print("-" * 30)

    # Example 3: Different lengths (Biopython handles this via alignment)
    protein_seq_g_bio = "MVLSPAD"
    protein_seq_h_bio = "MVLSPADKTNVKAAWGKVGAHAG"
    similarity_bio = calculate_protein_similarity_biopython(protein_seq_g_bio, protein_seq_h_bio)
    print(f"Biopython Similarity between '{protein_seq_g_bio}' and '{protein_seq_h_bio}': {similarity_bio:.2f}%")
    print("-" * 30)

    # Example 4: Empty sequences (Biopython)
    protein_seq_i_bio = ""
    protein_seq_j_bio = ""
    similarity_bio = calculate_protein_similarity_biopython(protein_seq_i_bio, protein_seq_j_bio)
    print(f"Biopython Similarity between '{protein_seq_i_bio}' and '{protein_seq_j_bio}': {similarity_bio:.2f}%")
    print("-" * 30)
