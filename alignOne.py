import sys, fasta, profileHMM


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

    # load db
    sseq = fasta.load(protSeqFN)[0][1]

