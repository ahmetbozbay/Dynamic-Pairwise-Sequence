# Ahmet Bozbay 150119861
# Anıl Batuhan Aslan 150119656
# Çağan Kurt 150119677

for i in range(1,6):
    with open('test'+str(i)+'.seq', 'r') as f:
        data = f.readlines()
        firstdnapiece=data[0]
        seconddnapiece = data[-1]
    def dynamic_sequence_alignment(seq1, seq2, match_score=3.0, mismatch_score=-1.0, open_gap_score=-1, extend_gap_score=0.5):
        m = len(seq1)
        n = len(seq2)
        # Initialize the matrix with all zeros
        dp = [[0 for _ in range(n+1)]
              for _ in range(m+1)]
        # Fill in the matrix using the dynamic programming approach
        for i in range(1, m+1):
            for j in range(1, n+1):
                # Calculate the scores for aligning the two characters
                match = dp[i-1][j-1] + match_score \
                    if seq1[i-1] == seq2[j-1] \
                    else dp[i-1][j-1] + mismatch_score
                delete = dp[i-1][j] + open_gap_score + extend_gap_score
                insert = dp[i][j-1] + open_gap_score + extend_gap_score
                # Take the maximum of the three scores
                dp[i][j] = max(match, delete, insert)
        # Initialize variables to store the alignment
        align1 = ""
        align2 = ""
        # Traceback to retrieve the optimal alignment
        i = m
        j = n
        while i > 0 or j > 0:
            if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
            elif i > 0 and dp[i][j] == dp[i-1][j] + open_gap_score + extend_gap_score:
                align1 = seq1[i-1] + align1
                align2 = "-" + align2
                i -= 1
            else:
                align1 = "-" + align1
                align2 = seq2[j-1] + align2
                j -= 1
        # Return the alignments and the score
        return align1, align2, dp[m][n]
    # Test the function
    align1, align2, score = dynamic_sequence_alignment(firstdnapiece, seconddnapiece, match_score=3.0, mismatch_score=-1.0, open_gap_score=-1, extend_gap_score=0.5)
    print(" *-*-*-*-*-*-*-*"'         Test'+str(i)+'.seq               '"*-*-*-*-*-*-*-*\n")
    print(' First given Dna align is  : ',align1,'\n Second given Dna align is : ',align2,'\n Calculated dynamic sequence alignment socre is : ',score,"\n\n")