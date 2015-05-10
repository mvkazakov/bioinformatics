import unittest
from code_challenge import *

class TestCodeChallenge(unittest.TestCase):

    def test_Skew(self):
        
        for p in RunFunWithFile(Skew, 'Skew.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)

    def test_MinSkewPositions(self):
        
        for p in RunFunWithFile(MinSkewPositions, 'MinSkewPositions.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)

    def test_HammingDistance(self):
        
        for p in RunFunWithFile(HammingDistance, 'HammingDistance.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
    def test_ApproximatePatternMatch(self):
        
        for p in RunFunWithFile(ApproximatePatternMatch, 'ApproximatePatternMatch.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
    def test_PatternCount(self):
        
        for p in RunFunWithFile(PatternCount, 'PatternCount.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
    def test_ApproximatePatternCount(self):
        
        for p in RunFunWithFile(ApproximatePatternCount, 'ApproximatePatternCount.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
    def test_Neighbors(self):
        
        for p in RunFunWithFile(Neighbors, 'Neighbors.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
    def test_FrequentWordsWithMismatches(self):
        
        for p in RunFunWithFile(FrequentWordsWithMismatches, 'FrequentWordsWithMismatches.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
    def test_FrequentWordsWithMismatchesReverseComplement(self):
        
        for p in RunFunWithFile(FrequentWordsWithMismatchesReverseComplement, 'FrequentWordsWithMismatchesReverseComplement.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
    def test_MotifEnumeration(self):
        
        for p in RunFunWithFile(MotifEnumeration, 'MotifEnumeration.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
    def test_MedianString(self):
        
        for p in RunFunWithFile(MedianString, 'MedianString.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
    def test_DistanceBetweenPatternAndStrings(self):
        
        for p in RunFunWithFile(DistanceBetweenPatternAndStrings, 'DistanceBetweenPatternAndStrings.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
    def test_MostProbableKmer(self):
        
        for p in RunFunWithFile(MostProbableKmerStr, 'MostProbableKmer.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
    def test_GreedyMotifSearch(self):
        
        for p in RunFunWithFile(GreedyMotifSearch, 'GreedyMotifSearch.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
    def test_GreedyMotifSearchLaplace(self):
        
        for p in RunFunWithFile(GreedyMotifSearchLaplace, 'GreedyMotifSearchLaplace.txt'):
            with self.subTest(p=p):
                self.assertEqual(*p)
        
if __name__ == '__main__':
    unittest.main()
