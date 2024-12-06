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
    seq = fasta.load(protSeqFN)[0][1]
    print(seq)

    # create profile HMM based on seed alignment
    seedAlignL = [seq for header,seq in fasta.load(seedAlignFN)]
    profile = profileHMM.Hmm(seedAlignL)
    emissionProbs = profile.emissions
    transitionProbs = profile.transitions
    print((transitionProbs['I57>D58']))
    posTrans = profile.posTransitions
    states = profile.states
    score, bestpath = viterbi.logOdds(seq, 5, transitionProbs, emissionProbs, states, posTrans)
    print("The score is", score)
    print("The best path is ", bestpath)

    # load db
    sseq = fasta.load(protSeqFN)[0][1]

