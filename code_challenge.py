import re
import operator
import random
import itertools
import bisect

bases = ['A', 'C', 'G', 'T']


def RunFunWithFile(fun, file):
    out = []
    with open(file) as f:
        args = []
        for line in f:
            line = line.rstrip("\n")
            if line != '':
                if re.search('^[+-]?\d+$', line):
                    line = int(line)
                args.append(line)
                # print(args)
            else:
                result = fun(*args[:-1])
                estimate = args[-1]
                if isinstance(result, list):
                    result = ' '.join(map(str, result))
                elif isinstance(result, set):
                    estimate = set(estimate.split())
                else:
                    result = str(result)
                # print(result)
                if isinstance(result, str) and re.search('^[+-]?\d+$', result):
                    result = int(result)
                out.append((result, estimate))
                args = []
    return out


def PatternCount(Text, Pattern, cmp=operator.eq):
    """Counting Pattern in the Text. Substrings are comparing by cmp
    function, it is the equvalence by default."""
    count = 0
    for i in range (0, len (Text) - len (Pattern) + 1):
        if cmp(Text[i:i + len (Pattern)], Pattern):
            count = count + 1
    return count

def ComputingFrequencies(Text, k):
    """Compute a frequencies of k-mers in the string Text"""

    Frequencies = {}
    
    for i in range (len (Text) - k + 1):
        Pattern = Text[i:(i + k)]
        if not Pattern in Frequencies:
            Frequencies[Pattern] = 1
        else:
            Frequencies[Pattern] = Frequencies[Pattern] + 1

    return Frequencies
    
def FrequentWords(Text, k):
    """Finding a most frequent k-mers in the string Text"""

    Count = ComputingFrequencies(Text, k)
    FrequentPatterns = set()
    

    maxCount = max(Count.values())
    
    for Pattern in Count:
        if Count[Pattern] == maxCount:
            FrequentPatterns.add(Pattern)

    return FrequentPatterns


def ReverseComplement(Text):
    """Return reverse complement of the Text"""

    return Text[::-1].translate(''.maketrans('ATGC', 'TACG'))
        
def PatternPositions(Pattern, Genome):
    """All starting positions where Pattern appears as a substring of Genome"""

    return [format(m.start()) for m in re.finditer('(?=' + Pattern + ')', Genome)]
        
    
def ClumpFind(Genome, k, L, t):
    """All distinct k-mers forming (L, t)-clumps in Genome"""

    Clumps = set()
    for i in range(len(Genome) - L):
        Count = ComputingFrequencies(Genome[i:i + L], k)
        for Pattern in Count:
            if Count[Pattern] >= t:
                Clumps.add(Pattern)

    return Clumps

def Skew(Genome):
    """Returns Skew sequence of Genome"""

    s = [0]

    for Letter in Genome:
        if Letter == 'G':
            s.append(s[-1] + 1)
        elif Letter == 'C':
            s.append(s[-1] - 1)
        else:
            s.append(s[-1])

    return s

def MinSkewPositions(Genome):
    """Find a position in a genome minimizing the skew"""

    s=Skew(Genome)
    m=min(s)
    return [i for i, v in enumerate(s) if v == m]

def HammingDistance(p, q):
    """Compute the Hamming distance between two strings"""

    if len(p) != len(q):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(p, q))
               
    
def ApproximatePatternMatch(Pattern, Text, d):
    """Find all approximate occurrences of a pattern in a string"""

    out = []
    for i in range(0, len(Text) - len(Pattern) + 1):
        if HammingDistance(Text[i:i + len(Pattern)], Pattern) <= d:
            out.append(i)

    return out

def ApproximatePatternCount(Text, Pattern, d):
    """Counting the total number of occurences of Patter in Text with
    d mismatches at most."""

    return PatternCount(Text, Pattern,
                        cmp=lambda a, b:
                        operator.le(HammingDistance(a, b), d))

def Neighbors(Pattern, d):
    """All neighborhoods of Pattern with at most d mismatches"""

    if d == 0:
        return [Pattern]

    if len(Pattern) == 1:
        return bases

    Neighborhood = set()

    SuffixNeighbors = Neighbors(Pattern[1:], d)

    for Text in SuffixNeighbors:
        if HammingDistance(Pattern[1:], Text) < d:
            for b in bases:
                Neighborhood.add(b + Text)
        else:
            Neighborhood.add(Pattern[0] + Text)

    return Neighborhood

def FrequentWordsWithMismatches(Text, k, d, Reverse=False):
    """All most frequent k-mers with up to d mismatches in Text"""

    FrequentPatterns = set()
    Neighborhoods = []
    Frequencies = {}

    for i in range(0, len(Text) - k + 1):
        Neighborhoods.append(list(Neighbors(Text[i:i + k], d)))
        if Reverse:
            Neighborhoods.append(
                list(Neighbors(ReverseComplement(Text[i:i + k]), d)))

    Neighborhoods = sum(Neighborhoods, [])
    
    for Pattern in Neighborhoods:
        if not Pattern in Frequencies:
            Frequencies[Pattern] = 1
        else:
            Frequencies[Pattern] = Frequencies[Pattern] + 1

    maxCount = max(Frequencies.values())
    
    for Pattern in Frequencies:
        if Frequencies[Pattern] == maxCount:
            FrequentPatterns.add(Pattern)

    return FrequentPatterns
    
def FrequentWordsWithMismatchesReverseComplement(Text, k, d):
    """All k-mers Pattern maximizing the sum
    Count_d(Text, Pattern)+ Count_d(Text, ReverseComplement(Pattern)) """

    return FrequentWordsWithMismatches(Text, k, d, Reverse=True)

def MotifEnumeration(k, d, *Dna):
    """A brute force algorithm for motif finding"""

    Patterns = set()

    for Text in Dna:
        for i in range(len(Text) - k + 1):
            for kmer in Neighbors(Text[i:i + k], d):
                Count = 0
                for Text in Dna:
                    if PatternCount(Text, kmer,
                                    lambda a, b:
                                    HammingDistance(a, b) <= d) > 0:
                        Count = Count + 1
                if Count >= len(Dna):
                    Patterns.add(kmer)


    return Patterns

def MakeAllKmers(k):
    """Returns all k-mers. Not works for k > 10"""

    if k == 1:
        return bases
    else:
        out = set()
        for kmer in MakeAllKmers(k - 1):
            for n in bases:
                out.add(n + kmer)
    return out

def DistanceBetweenPatternAndStrings(Pattern, *Dna):
    """The sum for each Text in Dna of the minimum Hamming distances
    between Pattern and all k-mers in Text"""

    k = len(Pattern)
    distance = 0

    for Text in Dna:
        hd = float('inf')

        for i in range(len(Text) - k + 1):
            hdNew = HammingDistance(Text[i:i + k], Pattern)
            if hd > hdNew:
                hd = hdNew

        distance = distance + hd

    return distance

def MedianString(k, *Dna):
    """The brute force solution for the Median String Problem"""

    Median = ""
    distance = float('inf')

    for Pattern in MakeAllKmers(k):

        d = DistanceBetweenPatternAndStrings(Pattern, *Dna)
        if distance > d:
            distance = d
            Median = Pattern

    return Median

def Probability(Text, Profile):
    """ """

    p = 1
    
    for n, i in zip(Text, range(len(Text))):
        p = p * Profile[n][i]

    return p

def MostProbableKmer(Text, k, Profile):
    """ """

    pmax = 0
    kmer = Text[0:k]
    for i in range(len(Text) - k + 1):
        p = Probability(Text[i:i + k], Profile)

        if p > pmax:
            pmax = p
            kmer = Text[i:i + k]

    return kmer

def MostProbableKmerStr(Text, k, *ProfileStr):
    """Wrapper around MostProbableKmer for testing"""
    Profile = {}
    for n, i in zip(bases, range(4)):
        Profile[n] = list(map(float, str.split(ProfileStr[i])))

    return MostProbableKmer(Text, k, Profile)

def GreedyMotifSearchLaplace(k, t, *Dna):
    return GreedyMotifSearchGeneral(k, t, True, *Dna)

def GreedyMotifSearch(k, t, *Dna):
    return GreedyMotifSearchGeneral(k, t, False, *Dna)

def GreedyMotifSearchGeneral(k, t, Laplace, *Dna):
    """ """

    BestMotifs = []
    bmScore = float('inf')
    
    for Text in Dna:
        BestMotifs.append(Text[0:k])

#    for m in range(t):

    m = 0

    for i in range(len(Dna[m]) - k + 1):
        Motifs = []
        Motifs.append(Dna[m][i:i + k])


        for j in set(range(t)) - set([m]):
            Profile = GetProfile(Motifs, Laplace)
            Motifs.append(MostProbableKmer(Dna[j], k, Profile))


        Score = DistanceBetweenPatternAndStrings(Consensus(Profile), *Motifs)
        if Score < bmScore:
            BestMotifs = Motifs
            bmScore = Score

    return BestMotifs

def GetProfile(Motifs, Laplace=False):
    """ """

    
    Profile = {}
    Freq = {}
    
    for n in bases:
        Profile[n] = []

        
    for i in range(len(Motifs[0])):
        for n in bases:
            Freq[n] = 0
        t = len(Motifs)
        for j in range(t):
            n = Motifs[j][i]
            Freq[n] = Freq[n] + 1

        if Laplace:
            for n in Freq:
                Freq[n] = Freq[n] + 1
                t = 2 * t
        for n in Freq:
            Profile[n].append(Freq[n]/t)

    return Profile

def Consensus(Profile):
    """ """

    c = ''
    for i in range(len(Profile['A'])):
        p = 0
        for n in Profile:
            if p < Profile[n][i]:
                p = Profile[n][i]
                m = n
        c = c + m
    return c

def Score(Motifs):
    """ """

    return DistanceBetweenPatternAndStrings(Consensus(GetProfile(Motifs, True)),
                                            *Motifs)

def RandomMotifs(k, *Dna):
    """ """
    Motifs = []
    #random.seed()
    
    for Text in Dna:
        
        i = random.choice(range(len(Text) - k + 1))
        Motifs.append(Text[i:i + k])

    return Motifs



def RandomizedMotifSearch(k, t, *Dna):
    """ """

    BestMotifs = RandomMotifs(k, *Dna)

    Motifs = BestMotifs[:]
    
    while True:

        Profile = GetProfile(Motifs, True)

        Motifs = []

        for Text in Dna:
            Motifs.append(MostProbableKmer(Text, k, Profile))
        #print(Score(Motifs))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs[:]
        else:
            break
                          
    return BestMotifs

def RandomizedMotifSearchTimes(times, k, t, *Dna):
    return MotifSearchTimes(RandomizedMotifSearch, times, k, t, *Dna)

def GibbsSamplerTimes(times, k, t, N, *Dna):
    return MotifSearchTimes(GibbsSampler, times, k, t, N, *Dna)
                          
def MotifSearchTimes(fun, times, *args):
    """ """

    result = []
    
    for i in range(times):
        result.append(fun(*args))

    #print(list(map(Score, result)))
    return min(result, key=Score)
    #return result

def GibbsSampler(k, t, N, *Dna):
    """ """

    BestMotifs = RandomMotifs(k, *Dna)

    Motifs = BestMotifs[:]
#    print(Motifs)

    for j in range(1, N):

        
        #random.seed()
        i = random.choice(range(0, len(Motifs)))
        before = Motifs[0:i]
        after = Motifs[i + 1:len(Motifs)]
        Profile = GetProfile(before + after, True)

        Motifs = before + [ProfileRandomlyKmer(Dna[i], k, Profile)] + after
#        Motifs = before + [MostProbableKmer(Dna[i], k, Profile)] + after

        #print(i)
        #print(Motifs)

        score = Score(Motifs)
        bestScore = Score(BestMotifs)
#        print(score)
#        print(bestScore)
        if score < bestScore:
            BestMotifs = Motifs[:]

    return BestMotifs


def Random(P):
    """ """

    cumdist = list(itertools.accumulate(P))
    x = random.random() * cumdist[-1]

    return bisect.bisect(cumdist, x)

def ProfileRandomlyKmer(Text, k, Profile):
    """ """

    P = []

    for i in range(len(Text) - k + 1):

        P = P + [Probability(Text[i:i + k], Profile)]

    i = Random(P)
    #print(P)
    #print(i)
    return Text[i:i + k]
