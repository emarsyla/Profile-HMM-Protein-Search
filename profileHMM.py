class Hmm:
    def __init__(self, alignment):
        self.transitions, self.emissions = self.makeProfile(alignment)
        numMatches = self.getNumMatches(alignment)
        self.states = self.getStates(numMatches)

    def getNumMatches(self, alignment):
        # Match state if > 50% of column is aligned
        numMatches = 0
        numAligned = len(alignment)
        seqlen = len(alignment[0])
        for i in range(seqlen):
            colscore = 0
            for j in range(numAligned):
                if alignment[j][i] != '.':
                    colscore += 1
            if colscore > (0.5 * numAligned):
                numMatches += 1
        return numMatches
            
    def getTransitions(self, numMatches):
        possibleTransitions = []
        for i in range(numMatches + 1):
            M = 'M' + str(i)
            I = 'I' + str(i)
            D = 'D' + str(i)
            nextM = 'M' + str(i+1)
            nextI = 'I' + str(i)
            nextD = 'D' + str(i+1)
            if i == 0:
                M = 'B0'
                D = 'DNE'
            if i == numMatches:
                nextM = 'E'
                nextD = 'DNE'
            for first in [M, I, D]:
                for second in [nextM, nextI, nextD]:
                    trans = first + second
                    if 'DNE' not in trans:
                        possibleTransitions.append(trans)
        return possibleTransitions
    
    def getEmissions(self, numMatches):
        aalist = ['G', 'P', 'A', 'V', 'L', 'I', 'M', 'C', 'F', 'Y', 'W', 'H', 'K', 'R', 'Q', "N", 'E', 'D', 'S', 'T']
        posEms = []
        posEmissions = []
        if numMatches >= 1:
            posEms.append('I0')
            for i in range(numMatches):
                M = 'M' + str(i + 1)
                I = 'I' + str(i + 1)
                for emis in [M, I]:
                    posEms.append(emis)
        for aa in aalist:
            for emis in posEms:
                posEmissions.append(emis + aa)
        return posEmissions
    
    def getStates(self, numMatches):
        states = ['B0', 'I0']
        for i in range(numMatches):
            states.append('M' + str(i+1))
            states.append('D' + str(i+1))
            states.append('I' + str(i+1))
        states.append('E')
        return states
                    
    def makeProfile(self, alignment):
        # List of amino acids
        aalist = ['G', 'P', 'A', 'V', 'L', 'I', 'M', 'C', 'F', 'Y', 'W', 'H', 'K', 'R', 'Q', "N", 'E', 'D', 'S', 'T']
        insertEmisList = []
        # Init emission and transition dictionaries
        emissionsProbs = {}
        transitionsProbs = {}
        # Initialize variables
        numAligned = len(alignment)
        seqlen = len(alignment[0])
        numMatches = 0
        seqPos = 0
        modelPos = 0
        curStates = ['B0' for _ in range(numAligned)]
        transitions = []
        insertColScore = 0
        # Iterate until we reach the end of the sequence
        while seqPos <= seqlen:
            prevStates = curStates
            curStates = []
            col = []
            colScore = 0
            if seqPos == seqlen:
                # We are at end of string
                modelPos += 1

                for aa in aalist:
                    # Emis prob for previous insert state
                    numEmisofAA = insertEmisList.count(aa) + 1
                    numEmisofAllAA = insertColScore + 20
                    emis = 'I' + str(modelPos - 1) + aa
                    prob = numEmisofAA/numEmisofAllAA
                    emissionsProbs[emis] = prob

                curStates = ['E' for _ in range(numAligned)]
                for i in range(numAligned):
                        if prevStates[i] == curStates[i]:
                            if 'I' in prevStates[i] and 'I' in curStates[i]:
                                transitions.append(prevStates[i] + curStates[i])
                        else:
                            transitions.append(prevStates[i] + curStates[i])
                
                # Get pos transitions
                possibleTransitions = []
                prevM = 'M' + str(modelPos - 1)
                prevI = 'I' + str(modelPos - 1)
                prevD = 'D' + str(modelPos - 1)
                curM = 'E'
                curI = 'I' + str(modelPos - 1)
                
                for first in [prevM, prevI, prevD]:
                    for second in [curM, curI]:
                        trans = first + second
                        if 'DNE' not in trans:
                            possibleTransitions.append(trans)
                
                # Find transition probabilities
                numStatesToTransTo = 2
                for transition in possibleTransitions:
                    init = transition[:2]
                    newList = [trans for trans in transitions if (trans[:2] == init)]
                    countTransSpecState = newList.count(transition) + 1
                    countTransAnyState = len(newList) + numStatesToTransTo
                    transitionsProbs[transition] = countTransSpecState/countTransAnyState
                transitions = []

            else:
                for j in range(numAligned):
                    col.append(alignment[j][seqPos])
                    if alignment[j][seqPos] != '.':
                        colScore += 1
                if colScore > (0.5 * numAligned):
                    # Match state, this means we are moving through the model
                    numMatches += 1
                    modelPos += 1
                    # Find states within match state
                    for i in range(numAligned):
                        if col[i] == '.':
                            curStates.append('D' + str(modelPos))
                        else:
                            curStates.append('M' + str(modelPos))
                    # Find emission probabilities
                    for aa in aalist:
                        # Emis prob for current M state
                        numEmisofAA = col.count(aa) + 1
                        numEmisofAllAA = colScore + 20
                        emis = 'M' + str(modelPos) + aa
                        prob = numEmisofAA/numEmisofAllAA
                        emissionsProbs[emis] = prob

                        # Emis prob for previous insert state
                        numEmisofAA = insertEmisList.count(aa) + 1
                        numEmisofAllAA = insertColScore + 20
                        emis = 'I' + str(modelPos - 1) + aa
                        prob = numEmisofAA/numEmisofAllAA
                        emissionsProbs[emis] = prob

                    # Get transitions in this state    
                    for i in range(numAligned):
                        if prevStates[i] == curStates[i]:
                            if 'I' in prevStates[i] and 'I' in curStates[i]:
                                transitions.append(prevStates[i] + curStates[i])
                        else:
                            transitions.append(prevStates[i] + curStates[i])

                    # Get all possible transitions
                    possibleTransitions = []
                    prevM = 'M' + str(modelPos - 1)
                    prevI = 'I' + str(modelPos - 1)
                    prevD = 'D' + str(modelPos - 1)
                    curM = 'M' + str(modelPos)
                    curI = 'I' + str(modelPos - 1)
                    curD = 'D' + str(modelPos)
                    if modelPos == 1:
                        prevM = 'B0'
                        prevD = 'DNE'
                    for first in [prevM, prevI, prevD]:
                        for second in [curM, curI, curD]:
                            trans = first + second
                            if 'DNE' not in trans:
                                possibleTransitions.append(trans)
                    
                    # Find transition probabilities
                    numStatesToTransTo = 3
                    for transition in possibleTransitions:
                        init = transition[:2]
                        newList = [trans for trans in transitions if (trans[:2] == init)]
                        countTransSpecState = newList.count(transition) + 1
                        countTransAnyState = len(newList) + numStatesToTransTo
                        transitionsProbs[transition] = countTransSpecState/countTransAnyState
                    insertColScore = 0
                    insertEmisList = []
                    transitions = []

                else:
                    # Non match state
                    insertEmisList.append(col)
                    insertColScore = insertColScore +colScore
                    for i in range(numAligned):
                        if col[i] == '.':
                            curStates.append(prevStates[i])
                        else:
                            curStates.append('I' + str(modelPos))
                    for i in range(numAligned):
                        if prevStates[i] == curStates[i]:
                            if 'I' in prevStates[i] and 'I' in curStates[i]:
                                transitions.append(prevStates[i] + curStates[i])
                        else:
                            transitions.append(prevStates[i] + curStates[i])
            seqPos += 1
        return transitionsProbs, emissionsProbs
            
                


#alignment = ["IV..EN", "IV...D", "LSKYEN", "IS..PD", "I....D"]  
#profileHMM = Hmm(alignment)
#print(profileHMM.states)