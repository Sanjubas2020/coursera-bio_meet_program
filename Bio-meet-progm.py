# Bio-meet-propm.-pyhton-codes
Biology Meets Programming: Bioinformatics for Beginners Completed

Week 1
vibrio_cholera = ("ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC")
print (len(vibrio_cholera))


def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 
    
  
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i + len(Pattern)] == Pattern:
            count = count + 1
    return count

Text ="ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"
Pattern = "TGATCA"
# Finally, print the result of calling PatternCount on Text and Pattern.
# Don't forget to use the notation print() with parentheses included!
print (PatternCount(Text, Pattern))


def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        for i in range(n-k+1):
            if Text[i:i+k] == Pattern:
                freq[Pattern] = freq[Pattern] + 1
    return freq
    
    
    def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern]= PatternCount(Text,Pattern)
        
    return freq

def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
            words.sort()
    return words
    
    
def Reverse(Pattern):
  reversed_Pattern=Pattern[::-1]
  return reversed_Pattern

def Complement(Pattern):
    my_complements = {"A" : "T", "T":"A", "C":"G", "G":"C"}
    new_string = "".join([my_complements[Pattern[i]] for i in range(len(Pattern))])
    return new_string

def ReverseComplement(Pattern): 
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern
def Reverse(Pattern):
  reversed_Pattern=Pattern[::-1]
  return reversed_Pattern
def Complement(Pattern):
    my_complements = {"A" : "T", "T":"A", "C":"G", "G":"C"}
    new_string = "".join([my_complements[Pattern[i]] for i in range(len(Pattern))])
    return new_string
    
    
    def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    for i in range(len(Genome)-len(Pattern)+1):
        if Pattern == Genome[i:i+len(Pattern)]:
            positions.append(i)
   
    return positions
    
    import sys                             # needed to read the genome
input = sys.stdin.read().splitlines() #
v_cholerae = input[1] 

def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    for i in range(len(Genome)-len(Pattern)+1):
        if Pattern == Genome[i:i+len(Pattern)]:
            positions.append(i)
    return positions
print (PatternMatching('CTTGATCAT', v_cholerae))

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count


Text = "AACTCTATACCTCCTTTTTGTCGAATTTGTGTGATTTATAGAGAAAATCTTATTAACTGAAACTAAAATGGTAGGTTTGGTGGTAGGTTTTGTGTACATTTTGTAGTATCTGATTTTTAATTACATACCGTATATTGTATTAAATTGACGAACAATTGCATGGAATTGAATATATGCAAAACAAACCTACCACCAAACTCTGTATTGACCATTTTAGGACAACTTCAGGGTGGTAGGTTTCTGAAGCTCTCATCAATAGACTATTTTAGTCTTTACAAACAATATTACCGTTCAGATTCAAGATTCTACAACGCTGTTTTAATGGGCGTTGCAGAAAACTTACCACCTAAAATCCAGTATCCAAGCCGATTTCAGAGAAACCTACCACTTACCTACCACTTACCTACCACCCGGGTGGTAAGTTGCAGACATTATTAAAAACCTCATCAGAAGCTTGTTCAAAAATTTCAATACTCGAAACCTACCACCTGCGTCCCCTATTATTTACTACTACTAATAATAGCAGTATAATTGATCTGA"
count_2 = PatternCount(Text, "CTTGATCAT")
print (count_1 + count_2)


weeek 2

def SymbolArray(Genome, symbol):
    # type your code here
    array = {}
    n = len(Genome)
    ExtendGenome = Genome +Genome[0:n//2]
    for i in range(n):
        array [i] = PatternCount (ExtendGenome[i:i+n//2], symbol)
    return array


# Reproduce the PatternCount function here.
def PatternCount(Text, Pattern):
    # type your code here
    count = 0 
    n = len(Pattern)
    for i in range(len(Text)):
        if Text[i:i+n] == Pattern:
            count += 1
    return count
    
    # Input:  Strings Genome and symbol
Genome = "AAAAGGGG"
symbol = "A"
# Output: FasterSymbolArray(Genome, symbol)
def FasterSymbolArray(Genome, symbol):
    array = {}
    n= len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    
    #look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])
    
    for i in range(1, n):
        #start by setting the current array value equal to the previous array value
        array[i] = array[i - 1]
        
        #the current array value can differ from the previous a value by at most 1
        if ExtendedGenome[i - 1] == symbol:
            array[i] = array[i] - 1
        if ExtendedGenome[i + (n//2) - 1] == symbol:
            array[i] = array[i] + 1
    return array
   
   
   # Input:  A String Genome
# Output: The skew array of Genome as a list.
def SkewArray(Genome):
    # your code here
    n = len(Genome)
    Skew = [0]
    Symbol1 = 'G'
    Symbol2 = 'C'
    for i in range(1,n+1):
        Skew.append(Skew[i-1])
        if Genome[i-1] == Symbol1:
            Skew[i] = Skew[i-1] + 1
        elif Genome[i-1] == Symbol2:
            Skew[i] = Skew[i-1] - 1
    return Skew
    
    def SkewArray(Genome):
    array = [0]
    Skew = 0
    for i in Genome:
        if i == 'A' or i == 'T':
            Skew += 0
            array.append(Skew)
        if i == 'C':
            Skew -= 1
            array.append(Skew)
        if i == 'G':
            Skew += 1
            array.append(Skew)
    return array

def MinimumSkew(Genome):
    array = SkewArray(Genome)
    positions = []
    count = 0
    minarray = min(array)
    for i in array:
        if i == minarray:
            positions.append(count)
        count +=1
    return positionsre
    return skew
    
    # Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    count = 0
    for i, j in zip(p, q):
       if i != j:
           count += 1
    return count


def HammingDistance(p, q):
    count=0
    t1=p
    t2=q
    for i in range(len(p)):
        if t1[i]!=t2[i]:
            count=count+1
    return count
def ApproximatePatternMatching(Text, Pattern,d):
    positions = []
    count = 0
    t=Text
    p=Pattern
    k = len(p)
    l = len(t)
    for i in range(0,l-k+1):
        if HammingDistance((t[i:k+i]), p) <= d:
            positions.append(count)
        count += 1
    return positions
    
    
    week 3


def Count(Motifs):
    count = {} # initializing the count dictionary
    # your code here
    return count


def Count(Motifs):

    count = {} # initializing the count dictionary
    # your code here
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count
# Input: A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    # insert your code here
    profile1 = Count(Motifs)
    for key in "ACGT":
        for key in profile1:  
            profile[key] = [x / t for x in profile1[key]]
    return profile
    # Insert your Count(Motifs) function here.
def Count(Motifs):
    k = len(Motifs[0])
    count = {'A': [0]*k, 'T' : [0]*k, 'G' : [0]*k, 'C' : [0]*k}
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] +=1
    return count

# Input: A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    consensus = ''
    freq_letter = ''
    count = Count(Motifs)
    for j in range(k):
        list = [count['A'][j], count['T'][j], count['C'][j], count['G'][j]]
        list.sort()
        m = list[-1]
        for symbol in 'ACTG':
            if count[symbol][j] == m:
                freq_letter = ''.join(symbol)
        consensus+= freq_letter[0]
    return consensus
    
    def Consensus(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    consensus = ''
    freq_letter = ''
    count = Count(Motifs)
    for j in range(k):
        list = [count['A'][j], count['T'][j], count['C'][j], count['G'][j]]
        list.sort()
        m = list[-1]
        for symbol in 'ACTG':
            if count[symbol][j] == m:
                freq_letter = ''.join(symbol)
        consensus+= freq_letter[0]
    return consensus
# Copy your Count(Motifs) function here.
def Count(Motifs):
    k = len(Motifs[0])
    count = {'A': [0]*k, 'T' : [0]*k, 'G' : [0]*k, 'C' : [0]*k}
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] +=1
    return count
# Input: A set of k-mers Motifs

def Score(Motifs):
    consensus = Consensus(Motifs)
    count = 0
    for motif in Motifs:
        for index in range(len(motif)):
            if motif[index] != consensus[index]:
                count += 1
    return count
# Output: The score of these k-mers.

def Pr(Text, Profile):
    # insert your code here
    pr = 1
    t = len(Text)
    for i in range(t):
        pr = pr*Profile[Text[i]][i]
    return pr
    
    
    
def Pr(Text, profile):

    n = len(Text)

    p = 1

    for i in range (n):

        p *= profile[Text[i]][i]

    return p 
def ProfileMostProbableKmer(text, k, profile):

    n = len(text)

    PrMap = {}

    most_pr_kmer= []

    for i in range(n-k+1):

        k_mer = text[i:i+k]    

        PrMap[k_mer] = Pr(k_mer,profile)

        m = max(PrMap.values())     

    for k_mer in PrMap:

        if PrMap[k_mer] == m:

            most_pr_kmer.append(k_mer)

    return most_pr_kmer[0]     
    
    
    
    def Count(Motifs):
    count = {} # initializing the count dictionary
    # your code here
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

def Consensus(Motifs):
    # insert your code here
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = Count(Motifs)
    for i in 'ACTG':
        for j in range(k):
            profile[i][j] = profile[i][j]/t  
    return profile

def Score(Motifs):
    # Insert code here
    score = 0
    k = len(Motifs[0])
    count = Count(Motifs)
    max_symbol = Consensus(Motifs)
    sum1 = 0
    for j in range(k):
        m = 0
        for symbol in "ATCG":
            if count[symbol][j] > m:
                sum1 += count[symbol][j]
    for j in range(k):
        m = 0
        for symbol in "AGTC":
            if count[symbol][j] > m:
                m = count[symbol][j]
        score += m  
    return sum1-score

def Pr(Text, Profile):
    p=1
    k = len(Profile["A"])
    for i in range(len(Text)):
        p=p*Profile[Text[i]][i]
    return p


def ProfileMostProbablePattern(text,k,profile):
    p=-1
    result=text[0:k]
    for i in range(len(text)-k+1):
        seq=text[i:i+k]
        pr=Pr(seq,profile)
        if pr>p:
            p=pr
            result=seq
    return result

def GreedyMotifSearch(Dna,k,t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
    
    
    def Count(Motifs):

    count = {}

    k = len(Motifs[0])

    for symbol in "ACGT":

        count[symbol] = []

        for j in range(k):

            count[symbol].append(0)

    t = len(Motifs)

    for i in range(t):

        for j in range(k):

            symbol = Motifs[i][j]

            count[symbol][j] += 1

    return count

 

def Consensus(Motifs):

  

    k = len(Motifs[0])

    count = Count(Motifs)

    consensus = ""

    for j in range(k):

        m = 0

        frequentSymbol = ""

        for symbol in "ACGT":

            if count[symbol][j] > m:

                m = count[symbol][j]

                frequentSymbol = symbol

        consensus += frequentSymbol

    return consensus

 

def Profile(Motifs):

    t = len(Motifs)

    k = len(Motifs[0])

    profile = Count(Motifs)

    for i in 'ACTG':

        for j in range(k):

            profile[i][j] = profile[i][j]/t  

    return profile

 

def Score(Motifs):

    # Insert code here

    score = 0

    k = len(Motifs[0])

    count = Count(Motifs)

    max_symbol = Consensus(Motifs)

    sum1 = 0

    for j in range(k):

        m = 0

        for symbol in "ATCG":

            if count[symbol][j] > m:

                sum1 += count[symbol][j]

    for j in range(k):

        m = 0

        for symbol in "AGTC":

            if count[symbol][j] > m:

                m = count[symbol][j]

        score += m  

    return sum1-score

 

def Pr(Text, Profile):

    p=1

    k = len(Profile["A"])

    for i in range(len(Text)):

        p=p*Profile[Text[i]][i]

    return p

 

#Finally solved it

def ProfileMostProbablePattern(text,k,profile):

    p=-1

    result=text[0:k]

    for i in range(len(text)-k+1):

        seq=text[i:i+k]

        pr=Pr(seq,profile)

        if pr>p:

            p=pr

            result=seq

    return result

 

def GreedyMotifSearch(Dna,k,t):

    BestMotifs = []

    for i in range(0, t):

        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])

    for m in range(n-k+1):

        Motifs = []

        Motifs.append(Dna[0][m:m+k])

        for j in range(1, t):

            P = Profile(Motifs[0:j])

            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):

            BestMotifs = Motifs

    return BestMotifs

Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

 

# set t equal to the number of strings in Dna and k equal to 15

k = 15

t = 10

# Call GreedyMotifSearch(Dna, k, t) and store the output in a variable called Motifs

Motifs = GreedyMotifSearch(Dna, k, t)

# Print the Motifs variable

print(Motifs)

# Print Score(Motifs)

print(Score(Motifs))



def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {} # initializing the count dictionary
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count
    
    
    # Input: A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    cont=CountWithPseudocounts(Motifs)
    for i in range(k):
        su=0
        for symbol in "ACGT":
            su=su+cont[symbol][i]
        for symbol in "ACGT":
            cont[symbol][i] = cont[symbol][i]/su
    profile=cont
    return profile

# Input: A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
# HINT: You need to use CountWithPseudocounts as a subroutine of ProfileWithPseudocounts
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {} # output variable
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count
    
    
    
    def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
# Copy all needed subroutines here.  These subroutines are the same used by GreedyMotifSearch(),
# except that you should replace Count(Motifs) and Profile(Motifs) with the new functions
# CountWithPseudocounts(Motifs) and ProfileWithPseudocounts(Motifs).
def Score(Motifs):
    score = 0
    consensus = Consensus(Motifs)
    k = len(consensus)
    for j in range(len(Motifs)):
        for symbol in range(k):
            if Motifs[j][symbol] != consensus[symbol]:
                score += 1
    return score
def CountWithPseudocounts(Motifs):
    count = {} 
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = CountWithPseudocounts(Motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j]/(t+4))
    return profile
def Consensus(Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "AGCT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus
# Then copy your ProfileMostProbablePattern(Text, k, Profile) and Pr(Text, Profile) functions here.
def Pr(Text, profile):
    p = 1
    for i in range(len(Text)):
        p *= profile[Text[i]][i]
    return p
def ProfileMostProbablePattern(Text, k, Profile):
    most_prob = Text[0:k] 
    p_max = Pr(Text[0:k], Profile)
    for i in range(1, len(Text) - k + 1):
         if Pr(Text[i:i+k], Profile) > p_max:
                p_max = Pr(Text[i:i+k], Profile)
                most_prob = Text[i:i+k]        
    return most_prob
    
    
    def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
# Copy all needed subroutines here.  These subroutines are the same used by GreedyMotifSearch(),
# except that you should replace Count(Motifs) and Profile(Motifs) with the new functions
# CountWithPseudocounts(Motifs) and ProfileWithPseudocounts(Motifs).
def Score(Motifs):
    score = 0
    consensus = Consensus(Motifs)
    k = len(consensus)
    for j in range(len(Motifs)):
        for symbol in range(k):
            if Motifs[j][symbol] != consensus[symbol]:
                score += 1
    return score
def CountWithPseudocounts(Motifs):
    count = {} 
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = CountWithPseudocounts(Motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j]/(t+4))
    return profile
def Consensus(Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "AGCT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus
# Then copy your ProfileMostProbablePattern(Text, k, Profile) and Pr(Text, Profile) functions here.
def Pr(Text, profile):
    p = 1
    for i in range(len(Text)):
        p *= profile[Text[i]][i]
    return p
def ProfileMostProbablePattern(Text, k, Profile):
    most_prob = Text[0:k] 
    p_max = Pr(Text[0:k], Profile)
    for i in range(1, len(Text) - k + 1):
         if Pr(Text[i:i+k], Profile) > p_max:
                p_max = Pr(Text[i:i+k], Profile)
                most_prob = Text[i:i+k]        
    return most_prob
Dna= ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC','CACGTCAATCAC','CAATAATATTCG']
k = 3
t= 5
print (GreedyMotifSearchWithPseudocounts(Dna, k, t))


def Motifs(Profile, Dna):
    motifs = []
    t = len(Dna)
    k = 4
    for i in range(t):
        motif = ProfileMostProbablePattern(Dna[i], k, Profile)
        motifs.append(motif)
    return motifs
# Insert your ProfileMostProbablePattern(Text, k, Profile) and Pr(Pattern, Profile) functions here.
def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    return p

def ProfileMostProbablePattern(Text, k, Profile):
    p_dict = {}
    for i in range(len(Text)- k +1):
        p = Pr(Text[i: i+k], Profile)
        p_dict[i] = p
    m = max(p_dict.values())
    keys = [k for k,v in p_dict.items() if v == m]
    ind = keys[0]
    return Text[ind: ind +k]
    
    
    import random

def RandomMotifs(Dna, k, t):

    t = len(Dna)
    l = len(Dna[0])
    RandomMotif =[]
    for i in range(t):
        r = random.randint(1,l-k) # 1 is not added as it is inclusive of last element also
        RandomMotif.append(Dna[i][r:r+k])
    return RandomMotif
    
    import random

def randomMotifs(dna,k,t):

    kmm = []
    sc = []
    k = 3
    D = {}
    for i in range(0,len(dna)):
        km = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        D[i] = km
    for m in range(0,t):
        ran = random.randint(0,len(D[0])-1)
        kmm += [D[m][ran]]

    return kmm


def ProfileWithPseudocounts(Motifs):


    t = len(Motifs)
    k = len(Motifs[0])
    profile = CountWithPseudocounts(Motifs) # output variable
    for symbol in profile:
        for kk in range(0,len(profile[symbol])):
            profile[symbol][kk] = profile[symbol][kk]/(len(Motifs) + 4)

    return profile


def CountWithPseudocounts(Motifs):

    count = {}
    for i in 'ACGT':
        count[i] = []
        for ii in range(len(Motifs[0])):
            count[i].append(1)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count



def Score(Motifs):


    count = 0
    L = Consensus(Motifs)
    for i in Motifs:
        for chr1, chr2 in zip(i,L):
            if chr1 != chr2:
                count += 1
    return count


def Consensus(Motifs):



    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus


def Count(Motifs):

    count = {}
    for i in 'ACGT':
        count[i] = []
        for ii in range(len(Motifs[0])):
            count[i].append(0)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count


def RandomMotifs(dna,k,t):

    kmm = []
    sc = []
    D = {}
    for i in range(0,len(dna)):
        km = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        D[i] = km
    for m in range(0,t):
        ran = random.randint(0,len(D[0])-1)
        kmm += [D[m][ran]]

    return kmm

def Motifs(pf,dna):

    k = len(pf['A'])
    D = []
    for i in range(0,len(dna)):
        km = []
        sc = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        for i in km:
            sc += [Pr(i,pf)]
        D += [km[sc.index(max(sc))]]

    return D


def Pr(Text, Profile):

    p = 1
    for i in range(0,len(Text)):
        p *= Profile[Text[i]][i]

    return p


def RandomizedMotifSearch(Dna, k, t):

    M = RandomMotifs(Dna, k, t)
    BestMotifs = M

    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs
            
            
            
 def Normalize(P):
    d = {}
    for k,v in P.items():
        d[k] = P[k]/sum(P.values())
    return d

d = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}


import random

def WeightedDie(Probabilities):
    n = random.uniform(0, 1)
    for p in Probabilities:
        n -= Probabilities[p]
        if n <= 0:
            return p
            
            # first, import the random package
import random
from operator import itemgetter
# then, copy Pr, Normalize, and WeightedDie below this line
def WeightedDie(d):
    ran = random.uniform(0, 1)
    #print(ran,d)
    tot = 0
    for k, v in sorted(d.items(),key=itemgetter(1)):
        if tot <= ran < v + tot:
            return k
        tot += v
def Normalize(P):
    D = {}
    for k,v in P.items():
        D[k] = P[k]/sum(P.values())
    return D
def Pr(Text, Profile):
    p = 1
    for i in range(0,len(Text)):
        p *= Profile[Text[i]][i]
    return p
# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    # your code here
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)

    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)
    
    
    import random

 

# Input:  Integers k, t, and N, followed by a collection of strings Dna

# Output: GibbsSampler(Dna, k, t, N)

def GibbsSampler(Dna, k, t, N):

    BestMotifs = [] # output variable

    # your code here

    Motifs = RandomMotifs(Dna, k, t)

    BestMotifs = Motifs

    for j in range(1,N):

        i = random.randint(0,t-1)

        ReducedMotifs = []

        for j in range(0,t):

            if j != i:

                ReducedMotifs.append(Motifs[j])

            

        Profile = ProfileWithPseudocounts(ReducedMotifs)

        Motif_i = ProfileGeneratedString(Dna[i], Profile, k)

        

        Motifs[i] = Motif_i

    

        if Score(Motifs) < Score(BestMotifs):

                BestMotifs=Motifs

          

    

    return BestMotifs

 

# place all subroutines needed for GibbsSampler below this line

# Input:  A list of strings Dna, and integers k and t

# Output: RandomMotifs(Dna, k, t)

# HINT:   You might not actually need to use t since t = len(Dna), but you may find it convenient

def RandomMotifs(Dna, k, t):

    # place your code here.

    s = len(Dna[0])

    rm = []

    

    for i in range(0,t):

        init_index = random.randint(1,s-k)

        rm.append(Dna[i][init_index:init_index+k])    

    

    return rm

 

# Input:  A set of kmers Motifs

# Output: ProfileWithPseudocounts(Motifs)

def ProfileWithPseudocounts(Motifs):

    t = len(Motifs)

    k = len(Motifs[0])

    profile = {} # output variable

    # your code here

    c = CountWithPseudocounts(Motifs)

    for n in 'ACGT':

        p = []

        for i in range(0,k):

            p.append(c[n][i]/(t+4))

        profile[n] = p

    return profile

 

# Input:  A set of kmers Motifs

# Output: CountWithPseudocounts(Motifs)

def CountWithPseudocounts(Motifs):

    t = len(Motifs)

    k = len(Motifs[0])

    # insert your code here

    count = {} # initializing the count dictionary

    for symbol in "ACGT":

        count[symbol] = []

        for j in range(k):

             count[symbol].append(1)

    for i in range(t):

        for j in range(k):

             symbol = Motifs[i][j]

             count[symbol][j] += 1

    return count 

 

#tests in which of the intervals defined by list ar the number r lies

def testinterval(ar,r):

    ar.sort()

    

    if r<= ar[0]:

      return ar[0]

      

    for i in range(1,len(ar)-1):

      if ar[i-1]<r<=ar[i]:

        return ar[i]

        

    if ar[len(ar)-2]< r:

      return ar[len(ar)-1]

      

# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers

# Output: A randomly chosen k-mer with respect to the values in Probabilities

def WeightedDie(Probabilities):

    # your code here

    

    sumprob = {}

    

    s = 0

    

    for p in Probabilities:

        s += Probabilities[p]

        sumprob[p] = s

    

    

    revprob = {}

    

    for q in sumprob:

      revprob[sumprob[q]] = q

    

    w = list(sumprob.values())

    

    r = random.uniform(0,1)

   

    

    

    kmer = revprob[testinterval(w,r)]

    

    return kmer

    

# Input:  A string Text, a profile matrix Profile, and an integer k

# Output: ProfileGeneratedString(Text, profile, k)

def ProfileGeneratedString(Text, profile, k):

    # your code here

    n = len(Text)

    probabilities = {} 

 

    for i in range(0,n-k+1):

        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)

    

    probabilities = Normalize(probabilities)

    return WeightedDie(probabilities)

 

# Input:  String Text and profile matrix Profile

# Output: Pr(Text, Profile)

def Pr(Text, Profile):

    # insert your code here

    p = 1

    for i in range(0,len(Text)):

        p *= Profile[Text[i]][i]

 

    return p

    

    

# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)

# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities

def Normalize(Probabilities):

    # your code here

    result = {}

    sum = 0

    for m in Probabilities:

        sum += Probabilities[m]

        

    for n in Probabilities:

        result[n]= Probabilities[n]/sum

    

    return result  

    

# Input:  A set of k-mers Motifs

# Output: The score of these k-mers.

def Score(Motifs):

    # Insert code here

    k = len(Motifs[0])

    t = len(Motifs)

    cs = ConsensusWithPseudocounts(Motifs)

    

    score = 0

    for j in range(0,k):

        for i in range(0,t):

            if Motifs[i][j] != cs[j]:

                score += 1

        

    return score

 

# Input:  A set of kmers Motifs

# Output: A consensus string of Motifs.

def ConsensusWithPseudocounts(Motifs):

    # insert your code here

    k = len(Motifs[0])

    count = CountWithPseudocounts(Motifs)

    

    consensus = ""

    for j in range(k):

        m = 0

        frequentSymbol = ""

        for symbol in "ACGT":

            if count[symbol][j] > m:

                m = count[symbol][j]

                frequentSymbol = symbol

        consensus += frequentSymbol

        

    return consensus    

    

k = 8

t = 5

N = 100


Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA","GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG","TAGTACCGAGACCGAAAGAAGTATACAGGCGT","TAGATCAAGTTTCAGGTGCACGTCGGTGAACC","AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

    

print(GibbsSampler(Dna, k, t, N))  

            
            
            




    
