"""
Modified and Improved from the example from 
https://medium.freecodecamp.org/programming-the-genome-with-crispr-bd567a214e2a
"""

"""
Example function to find Cas9 targets in a given sequence. Recomend
starting with the main cas9_target_finder function and then reading
the helper functions.
"""


def dna_compliment(sequence):
    """
    Helper function that returns the compliment of a DNA sequence.
    """
    compliment = ''
    for nucleotide in sequence:
        if nucleotide == 'A':
            compliment += 'T'
        elif nucleotide == 'T':
            compliment += 'A'
        elif nucleotide == 'C':
            compliment += 'G'
        elif nucleotide == 'G':
            compliment += 'C'
    print "DNA Compliement: " + compliment
    return compliment

def reverse_dna_compliment(sequence):
    """
    Helper function that returns the reverse compliment of a DNA sequence.
    Useful for returning results from the non reference (bottom) strand
    while still keeping the 5' to 3' convention.
    """
    print "Reverse DNA Compliment: " + dna_compliment(sequence)[::-1]
    return dna_compliment(sequence)[::-1]

def enough_seq_context(pam_start, seq_len, guide_len, pam_len, strand):
    """
    Helper function that makes sure the program does not go out of bounds
    trying to find target sequences at the start or end of the sequence.
    """
    if strand == 'FWD' and pam_start > guide_len and pam_start + pam_len <= seq_len:
        return True
    if strand == 'BOT' and pam_start + pam_len + guide_len < seq_len:
        return True
    return False

def pam_matches(sequence, pam_seq):
    """
    Helper function that checks if sequence matches PAM.
    """
    if len(sequence) != len(pam_seq):
        return False
    for i, pam_char in enumerate(pam_seq):
        # all nucleatotides match the N constraint
        if pam_char == 'N':
            pass
        elif sequence[i] != pam_char:
            return False
    return True

# define constant variables for standard S. Pyogenes Cas9
# for a more flexible program these could be set as function parameters
PAM = 'NGG'
GUIDE_LEN = 20
PAM_LEN = len(PAM)
CUT_DIFF = 3 # difference between PAM start and cut pos

def cas9_target_finder(sequence):
    """
    Main function that returns all Cas9 targets in a given DNA sequence.
    Assumes that sequence is an all uppcase DNA sequence.
    (only 'A', 'T', 'C', or 'G').
    """
    seq_len = len(sequence)
    targets = []
    for i in range(seq_len):
        # find targets on reference (top) strand
        # check for out of bounds and correct PAM
        if (enough_seq_context(i, seq_len, GUIDE_LEN, PAM_LEN, 'FWD') and
                pam_matches(sequence[i:i+PAM_LEN], PAM)):
            targets.append({
                'cut_pos': i - CUT_DIFF,
                'pam_seq': sequence[i:i+PAM_LEN],
                'target_seq': sequence[i-GUIDE_LEN:i],
                'strand': 'forward'})
        # find targets on non-reference (bottom) strand
        # check for out bounds and correct PAM
        elif (enough_seq_context(i, seq_len, GUIDE_LEN, PAM_LEN, 'REV') and
              pam_matches(reverse_dna_compliment(sequence[i:i+PAM_LEN]), PAM)):
            # We take the reverse compliment of the pam_seq and target_seq to
            # keep the 5` to 3` convention
            targets.append({
                'cut_pos': i + PAM_LEN + CUT_DIFF,
                'pam_seq': reverse_dna_compliment(sequence[i:i+PAM_LEN]),
                'target_seq':
                    reverse_dna_compliment(sequence[i+PAM_LEN:i+PAM_LEN+GUIDE_LEN]),
                'strand': 'reverse'})
    print "cas9 target finder: " + targets
    return targets


if __name__ == '__main__':
    input = "input: " + "CCACGGTTTCTGTAGCCCCATACTTTGGATG"
    print input
    # cas9_target_finder("CCACGGTTTCTGTAGCCCCATACTTTGGATG")
    # dna_compliment("CCACGGTTTCTGTAGCCCCATACTTTGGATG")
    reverse_dna_compliment("CCACGGTTTCTGTAGCCCCATACTTTGGATG")