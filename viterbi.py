import math, profileHMM, random
#alignment = ["IV..EN", "IV...D", "LSKYEN", "IS..PD", "I....D"]  
#profile = profileHMM.Hmm(alignment)
#protSeq = "ASELNK"
#transD = profile.transitions
#emisD = profile.emissions
#states = profile.states

def viterbi(protSeq, transD, emisD, states):
    viterbi = {}
    backtrack = {}
    aa = protSeq[0]
    for state in states:
        viterbi[state] = [_ for _ in range(len(protSeq))]
        backtrack[state] = [_ for _ in range(len(protSeq))]
        if state + aa in emisD.keys():
            emisProb = emisD[state+aa]
        else:
            emisProb = 1
        viterbi[state][0] = math.log(0.1) * emisProb
        backtrack[state][0] = '-'
    for i in range(1,len(protSeq)):
        aa = protSeq[i]
        for curstate in states:
            maxProb = -float('inf')
            bestTrans = ''
            #probs = []
            for prevstate in states:
                if prevstate+curstate in transD.keys():
                    if curstate+aa in emisD.keys():
                        emisProb = emisD[curstate+aa]
                    else:
                        emisProb = 1
                    
                    prob = viterbi[prevstate][i-1] + math.log(transD[prevstate+curstate]*emisProb)
                    if prob>maxProb:
                        maxProb = prob
                        bestTrans = prevstate
            viterbi[curstate][i] = maxProb
            backtrack[curstate][i] = bestTrans
    finalProbs = []
    for state, probs in viterbi.items():
        finalProbs.append(probs[len(protSeq)-1])
    return max(finalProbs), backtrack, viterbi

def getBacktrack(protSeq, backtrack, viterbi, states):
    stateSeq = ''
    length = len(protSeq)
    finalMax = -float('inf')
    for state in states:
        score = viterbi[state][length - 1]
        if score > finalMax:
            finalMax = score
            curState = state
            prevState = backtrack[state][length-1]
    for i in range(length-1):
        stateSeq = curState + stateSeq
        curState = prevState
        prevState = backtrack[prevState][length - i - 2]
    return stateSeq
        
def logOdds(protSeq, reps, transD, emisD, states):
    realScore, backtrack, viter = viterbi(protSeq, transD, emisD, states)
    bestPath = getBacktrack(protSeq, backtrack, viter, states)
    scores = []
    for i in range(reps):
        shuffled = ''
        shuffled = shuffled.join(random.sample(protSeq, len(protSeq)))
        score, _, _ = viterbi(shuffled, transD, emisD, states)
        scores.append(score)
    final = realScore - (sum(scores)/reps)
    return final, bestPath

#score = viterbi(protSeq, transD, emisD, states)
#print(score)