import sys, fasta, profileHMM, viterbi


#### Main

if __name__ == "__main__":


    if len(sys.argv) !=3:
        sys.stderr.write("""
        Usage: python seedAlign.fa protSeq.fa

""")
        sys.exit(-1)

    seedAlignFN = sys.argv[1]
    protSeqFN = sys.argv[2]

    # create profile HMM based on seed alignment
    seedAlignL = [seq for header,seq in fasta.load(seedAlignFN)]
    profile = profileHMM.Hmm(seedAlignL)
    emissionProbs = profile.emissions
    transitionProbs = profile.transitions
    states = profile.states
    score, bestpath = viterbi.logOdds(protSeqFN, 5, transitionProbs, emissionProbs, states)
    print("The score is", score)
    print("The best path is ", bestpath)

    # load db
    sseq = fasta.load(protSeqFN)[0][1]

