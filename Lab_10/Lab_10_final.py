#!/usr/bin/env python
# coding: utf-8

# # Lab 10
# **Creator:** Lauren Enriquez <br>
# **Date:** November 13, 2019
# 
# The script will be run like this:<br>
# ```
# script.py scoringmatrix inputfasta
# ```
#    
# The format of this lab allows for the user to input the scoringmatrix and the inputfasta.
# 
# To test the script, the following files were used: <br>
# scoringmatrix: `BLOSUM62`<br>
# inputfasta: `example.fa`<br>

# In[7]:


# %load Todo_nw_aligner.py
#!/usr/bin/python

"""
Needleman-Wunsch Aligner
Bioengineering 131/231, Fall 2018

Command-line script to read in the contents of a multi-FASTA containing two
sequences and a score matrix, perform a global alignment, and print the
resulting alignment matrix and optimal alignment to STDOUT.

Your script should accept two positional command-line arguments:
(1) the path to the scoring matrix file
(2) the path to the multi-FASTA to be aligned.

It will be run like this:
   script.py scoringmatrix inputfasta
"""

import os
import sys

class NWAligner:
    def __init__(self, score_matrix_fname):
        self.score_matrix, self.gap_penalty = self.load_score_matrix(score_matrix_fname)

    @staticmethod
    def load_score_matrix(fname):
        """
        Input: (String) A path to a scoring matrix file.
        Output: (Dictionary) A nested dictionary where keys are strings
                and elements are scores as integers.
    
        Example:
    
        >>> matrix, gap_penalty = NWAligner.load_score_matrix('/home/bioe131/BLOSUM62')
        >>> matrix['A']['A']
        4
        >>> matrix['W']['W']
        11
        >>> gap_penalty
        -4

        """

        score_matrix = {}
        gap_penalty = None

        with open(fname) as fp:
            for line_numb, line in enumerate(fp):
                # ignore comments in matrix file
                if line.startswith("#"):
                    continue
                vals = line.split(" ") 
                if line_numb == 0:
                    for v in vals:
                        v = v.split("\n")[0]
                        temp_dict = {}
                        for a in vals:
                            a = a.split("\n")[0]
                            temp_dict[a] = None
                            score_matrix[v] = temp_dict
                            keys = list(score_matrix.keys())
                elif len(vals) == 1:
                    gap_penalty = int(vals[0])
                else:
                    index = 0
                    cKey = keys[line_numb-1]
                    for key in keys:
                        score_matrix[cKey][key] = int(vals[index])
                        index += 1

                ### TODO ###
                # Parse matrix file line-by-line and load into nested dictionaries.
                #
                # Last line of matrix contains the gap penalty which must be pulled
                # out and returned.

        return score_matrix, gap_penalty

    @staticmethod
    def load_FASTA(fname):
        """
        Input: (String) A path to a FASTA file containing exactly two sequences.
        Output: (List) A list containing two strings: one for each sequence.

        Example:

        >>> seqs = NWAligner.load_FASTA('example.fa')
        >>> seqs[0]
        'YAADSKATPGNPAFHQDEIFLARIAFIYQMWDGGQLKLIDYAPHHVMCEE'
        >>> seqs[1]
        'WVGQPNMKVQHWSNMKACCVKFITWTFIAPEKHACKWTETAYQADCDIIW'
        >>> len(seqs)
        2

        """

        seqs = []

        ### TODO ###
        # Load FASTA file and return list of sequences.
        # Throw an error if there are more than two sequences in the file.
        from Bio import SeqIO
        input_hand = open(fname, 'r')
        for record in SeqIO.parse(input_hand, 'fasta'):
            seqs.append(str(record.seq))
        if len(seqs) > 2:
            raise Exception("Greater than two sequences in the file")
        else:
            return seqs

    def align(self, seq_x, seq_y, print_matrix = False):
        """
        Input: (Strings) Two sequences to be aligned (seq_x and seq_y).
               (Boolean) If print_matrix is True, print the dynamic programming
                         matrix before traceback.
        Output: (Tuple of strings) Two sequences, aligned.

        Example:

        >>> aligner = NWAligner('BLOSUM62')
        >>> seqs = aligner.load_FASTA('example.fa')
        >>> aligner.align(seqs[0], seqs[1])
        ('YAAD-SKATPGNPAF---HQDEIF--L-AR--IA-FIYQM-WDGGQLK-LIDYAPH-HVM-C---E-------E---',
         'W---VGQ--P-N--MKVQH----WSNMKA-CCV-KFI---TW------TFI--APEKH--ACKWTETAYQADCDIIW')

        """

        ###
        ### INITIALIZATION
        ###

        # create two empty matrices with sizes based on the input sequences.
        # one contains the dynamic programming matrix, the other contains
        # pointers we'll use during traceback
        matrix = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]
        pointers = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]

        ### TODO ###
        # Fill the top row of the matrix with scores for gaps.
        # Fill the first column of the matrix with scores for gaps.
        for i in range(len(matrix)):
            matrix[i][0] = matrix[i-1][0] + self.gap_penalty
            pointers[i][0] = -1
        for i in range(len(matrix) -1):
            matrix[0][i] = matrix[0][i-1] + self.gap_penalty
            pointers[0][i] = 1
        ###
        ### RECURSION
        ###

        # fill the dynamic programming and pointer matrices
        for x in range(1, len(seq_x) + 1):
            for y in range(1, len(seq_y) + 1):
                match_score = self.score_matrix[seq_x[x - 1]][seq_y[y - 1]]
                diagMatch = matrix[x-1][y-1] + match_score
                leftGap = matrix[x-1][y] + self.gap_penalty
                upGap = matrix[x][y-1] + self.gap_penalty
                matrix[x][y] = max(diagMatch, leftGap, upGap)
                if diagMatch == matrix[x][y]:
                    pointers[x][y] = 0
                elif leftGap == matrix[x][y]: 
                    pointers[x][y] = -1
                else:
                    pointers[x][y] = 1


                ### TODO ###
                # Take the maximum score of three possibilities:
                #   1) The element in the matrix diagonal from this one
                #      plus the score of an exact match
                #   2) The element to the left plus a gap penalty
                #   3) The element above plus a gap penalty
                # ... and set the current element (matrix[x][y]) equal to that
                #
                # Keep track of which of these choices you made by setting
                #   the same element (i.e., pointers[x][y]) to some value that
                #   has meaning to you.

        # print the dynamic programming matrix
        if print_matrix:
            for x in range(len(seq_x) + 1):
                print (" ".join(map(lambda i: str(int(i)), matrix[x])))

        ###
        ### TRACEBACK
        ###

        # starting from the bottom right corner, follow the pointers back
        x, y = len(seq_x), len(seq_y)

        # fill these lists with the aligned sequences
        align_x = []
        align_y = []

        while x > 0 or y > 0:
            move = pointers[x][y]

            ### TODO ###
            # Follow pointers back through the matrix to the origin.
            # Depending on which "move" you made at each element in the
            #   matrix, you'll either align seq_x to seq_y, seq_x to a gap, or
            #   seq_y to a gap.
            
            #align seq_y to a gap
            if move == -1:
                align_x.append("-")
                align_y.append(seq_y[y-1])
                x -= 1
                
            #align seq_x to seq_y   
            elif move == 0:
                align_x.append(seq_x[x-1])
                align_y.append(seq_y[y-1])
                x -= 1
                y -= 1
                
            #align seq_x to a gap
            else:
                align_x.append(seq_x[x-1])
                align_y.append("-")
                y-= 1
                
        # flip the alignments, as they're reversed
        return ("".join(align_x[::-1]), "".join(align_y[::-1]))

###                                      ###
### NO NEED TO EDIT CODE BELOW THIS LINE ###
###                                      ###

#USER INPUT
scoringmatrix = input("Enter scoringmatrix: ")
inputfasta = input("Enter inputfasta: ")
sys.argv = ["-f", '%s' % scoringmatrix, '%s' % inputfasta]


if __name__ == '__main__':
    def usage():
        print ('usage: %s matrixfilename stringfilename')
        sys.exit(1)

    if len(sys.argv) != 3:
        usage()

    for fname in sys.argv[1:]:
        if not os.path.isfile(fname):
            print ('Can not open %s' % (fname,))
            usage()

    aligner = NWAligner(sys.argv[1])
    seqs = aligner.load_FASTA(sys.argv[2])
    result = aligner.align(seqs[0], seqs[1])

    print('>seq1\n%s\n>seq2\n%s' % (result[0], result[1]))


# In[8]:


get_ipython().system('ipython nbconvert --to=python Lab_10_final.ipynb')

