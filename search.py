import sys, fasta, profileHMM, viterbi


#### Main

if __name__ == "__main__":


    if len(sys.argv) !=4:
        sys.stderr.write("""
        Usage: python seedAlign.fa protDB.fa outfile.txt

""")
        sys.exit(-1)

    seedAlignFN = sys.argv[1]
    dbSeqFN = sys.argv[2]
    outFN = sys.argv[3]

    # create profile HMM based on seed alignment
    seedAlignL = [seq for header,seq in fasta.load(seedAlignFN)]
    profile = profileHMM.Hmm(seedAlignL)
    emissionProbs = profile.emissions
    transitionProbs = profile.transitions
    posTrans = profile.posTransitions
    states = profile.states
    


    # load db
    dbL = fasta.load(dbSeqFN)

    outL=[]
    for hd,sseq in dbL:

        # some stuff to get score here!
        score, bestpath = viterbi.logOdds(sseq, 5, transitionProbs, emissionProbs, states, posTrans)
        print("The score is", score)
        print("The best path is ", bestpath)
        outL.append((score,hd))

    outL.sort(reverse=True) # sort by score, high to low

    # write sorted output to file
    f=open(outFN,'w')
    for sc,hd in outL:
        print(sc,hd,file=f)
    f.close()
    
