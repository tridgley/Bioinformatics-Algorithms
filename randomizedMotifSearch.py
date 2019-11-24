#!/usr/bin/env python3
# Trevor Ridgley (tridgley)
# BME-206 Assignment 3: CRISPR Promoter Sequences
"""
Created on Tue Oct  8 00:42:12 2019

@author: tridgley

Usage: cat p1860Crisprs.txt paeCrisprs.txt | python3 randomizedMotifSearch.py > crisp1.txt -i 1000
       -k 13 -p 1 [-g] [-m] [-r]
       
Note to grader: Running with -g option results in considerably longer runtime because 100 internal
                iterations are performed for each user-specified iteration calling search method.
                Recommend -i 100 or -i 1000
                
Input: A fasta file or concatenation of fasta files containing fasta records.
Output: The Consensus promoter sequence comprising the fasta records, and its encoding cost.
                
randomizedMotifSearch.py reads genetic sequences upstream of bacterial CRISPR genes
    and searches for conserved promoter motifs of length -k over -i iterations using one of two 
    different approaches:
    1) Randomized search establishes a baseline profile using a random selection of motifs
       with length k, then calculates the entropy score for each collection of n motifs. The 
       randomized search will continue until the entropy score no longer improves (precision
       configurable using the PRECISION constant below). At this point, another iteration is
       attempted according to the -i option set by user.
    2) Gibbs search (-g option) establishes a baseline profile by randomly changing just one of
       the motifs with length k from the colletion of n motifs. The new candidate motif is selected
       from the deleted row using a special Gibbs die prior to re-profiling. While the entropy scores 
       are calculated the same way, we allow scores to worsen since potentially a better entropy
       can be achieved after some waiting. To avoid incremental improvements that would result
       in extraordinary runtimes, there is a hard cutoff of internal iterations that can be
       configured using the INTERNAL_GIBBS_ITERATIONS constant below.

    A baseline entropy score and motif can be established using randomly shuffled sequences via
       -r option
    More detailed view of the predicted motifs for each fasta record can be viewed using the
       -m option
"""
import sys, random, math, numpy as np

PRECISION = 0 # Determines the rate of descent for random search, higher = more iterations
INTERNAL_GIBBS_ITERATIONS = 100

class UpstreamSequences:
    """Collection of upstream sequences that contain hidden promoter candidates"""
    def __init__(self, k):
        self.headerList = list()
        self.sequenceList = list()
        self.kmer = k
        # self.gProfile = None
        
    def printLengths(self):
        for i in self.sequenceList:
            print(len(i))
        
    def addSequence(self, h, s):
        """Add DNA sequences from fasta file to the list of candidates to search"""
        self.headerList.append(h)
        self.sequenceList.append(s.upper())
    
    def search(self, pseudocount, gibbsFlag, candidates):
        """Search for conserved promoter sequences in a collection of fasta sequences
        
        Params:
            pseudocount = A real value that reflects confidence in the data
                          Higher values shift the distribution toward uniform
            gibbsFlag = Use the Gibbs method (True) instead of random search (False)
            candidates = Provides the starting point when using Gibbs method from which
                         just one sequence is selected for iteration.
        
        The profile data structure is a dictionary of base keys with Z as the sum:
        A: list(a1, a2, ..., ak)
        C: list(c1, c2, ..., ck)
        G: list(g1, g2, ..., gk)
        T: list(t1, t2, ..., tk)
        Z: list(sum(col1), sum(col2), ..., sum(colk)) 
        """
        n = len(self.sequenceList) # This is the number of fasta records to avoid recalculation
        m = len(self.sequenceList[0]) # This is the number of bases in each record to avod recalculation
        currentEntropy, minEntropy = (2*self.kmer, 2*self.kmer) # Max entropy = 2 * kmer len
        profile = dict()
        candidateKmers = list() # These are just integers of the starting motif coords
        gibbsIterations = int()
        randSeq = int()
        
        # Randomly seed the motif positions
        for i in range(n):
            # Iterate the k fasta sequences
            rand = random.randint(0, m-self.kmer)
            candidateKmers.append(rand)
           
        # Iterate the model with more data until the score shows no improvement
        while True:
            if gibbsFlag:
                if gibbsIterations < INTERNAL_GIBBS_ITERATIONS:
                    # Randomly remove a motif to calculate the new dice
                    randSeq = random.randint(0, n-1)
                    profile = self.updateGibbsProfile(candidateKmers, pseudocount, n, randSeq)
                    candidateKmers[randSeq] = self.gibbsRandom(profile, randSeq, m)
                    # print(candidateKmers)
                    # Then add the new motif back in to make another profile
                    profile = self.updateGibbsProfile(candidateKmers, pseudocount, n, n+1)
                    currentEntropy = self.entropyScore3(profile)
                    # print(currentEntropy)
                    if currentEntropy < minEntropy:
                        minEntropy = currentEntropy
                else:
                    break
                gibbsIterations += 1
            else:
                profile = self.updateProfile(candidateKmers, pseudocount, n)
                candidateKmers = self.motifGivenProfile(profile, n, m)
                currentEntropy = self.entropyScore3(profile)
                # print(currentEntropy)
                if round(currentEntropy, PRECISION) < round(minEntropy, PRECISION):
                    minEntropy = currentEntropy
                    profile = self.updateProfile(candidateKmers, pseudocount, n)
                else:
                    minEntropy = min(currentEntropy, minEntropy)
                    break
        # print('>>>' + str(minEntropy))
        return(minEntropy, candidateKmers)
    
    def gibbsRandom(self, prof, delSeq, m):
        """Roll the dice but with kmer score probabilities for choosing a replacement motif"""
        sideCounts = list()
        for i in range(m-self.kmer):
            # Iterate the kmer start in this sequence
            currentKmerScore = float(1)
            motif = self.sequenceList[delSeq][i:i+self.kmer]
            for j in range(0, self.kmer):
                # Now iterate the length of particular kmer and score it
                currentKmerScore *= prof[motif[j]][j]/prof['Z'][j]
            sideCounts.append(currentKmerScore)
        dieSum = sum(sideCounts)
        kmerDie = list()
        for i in range(m-self.kmer):
            # This should generate a die of length dieSum
            kmerDie += [i for j in range(int(10000*sideCounts[i] / dieSum))]     
        rand = random.randint(0, len(kmerDie)-1)
        motifToAdd = kmerDie[rand]
        return motifToAdd
        
    def updateGibbsProfile(self, motifs, pseudocount, n, skip):
        """Update the profile with the best motifs, but ignore the random choice for Gibbs.
        
        Params:
            motifs = the starting indices of kmers to use for scoring
            pseudocount = populate the profile with these values
            n = the number of fasta records
            skip = the fasta record left-out for this Gibbs iteration
                   (a value greater than n can be used to NOT skip)
        """
        profile = dict()
        for base in ['A','C', 'G', 'T']: # Z stores the total counts
            profile[base] = [pseudocount for i in range(self.kmer)]
            profile['Z'] = [pseudocount * 4 for i in range(self.kmer)]
        
        for i in range(n):
            if i != skip:
                # Iterate the n fasta sequences
                start = motifs[i]
                motif = self.sequenceList[i][start:start+self.kmer]
                for j in range(self.kmer):
                    # Then iterate the counts at each index
                    profile[motif[j]][j] += 1
                    profile['Z'][j] += 1
        return profile
            
    def updateProfile(self, motifs, pseudocount, n):
        """Update the profile with the best motifs after they were scored."""
        profile = dict()
        for base in ['A','C', 'G', 'T']: # Z stores the total counts
            profile[base] = [pseudocount for i in range(self.kmer)]
            profile['Z'] = [pseudocount * 4 for i in range(self.kmer)]
        for i in range(n): 
            # Iterate the n fasta sequences
            start = motifs[i]
            motif = self.sequenceList[i][start:start+self.kmer]
            for j in range(self.kmer):
                # Then iterate the counts at each index
                profile[motif[j]][j] += 1
                profile['Z'][j] += 1
        return profile
        
    def motifGivenProfile(self, prof, n, m):
        """Find best scoring (most probable) motifs given the profile
        
        Ex: Given profile (matrix) of A=1,C=1,G=2,T=1, score all 4mers in a row as 
            TTAC = (1/5)(1/5)(1/5)(1/5) = 0.0016
        """
        bestKmersInProfile = list()
        for seq in range(n): # range(0, 5):
            # First iterate each sequence in the set of fasta records
            bestKmerInRow = int()
            bestKmerScoreInRow = float()
            for i in range(m-self.kmer): # range(0, 191 - self.kmer):
                # Then iterate the kmer start in this sequence
                currentKmerScore = float(1)
                motif = self.sequenceList[seq][i:i+self.kmer]
                for j in range(0, self.kmer):
                    # Now iterate the length of particular kmer and score it
                    currentKmerScore *= prof[motif[j]][j]/prof['Z'][j]
                if currentKmerScore > bestKmerScoreInRow:
                    # Store the top probability and starting index of best kmer
                    bestKmerScoreInRow = currentKmerScore
                    bestKmerInRow = i
            bestKmersInProfile.append(bestKmerInRow)
        return bestKmersInProfile
    
    def entropyScore3(self, prof):
        """Calculate the entropy for a profile of n kmers
        
        Entropy (S) = SUM over k indices of the kmer profile
                      SUM over base in ACGT
                      -P(base)*log_2(Pbase)
        P(base) = Count at profile index i, j / Sum at profile index i
        """
        entropyTerms = list()
        for i in range(self.kmer):
            # First iterate each column of the kmer profile
            innerSum = float() # Pseudocounts can be non-int type
            for base in prof.keys():
                # Then iterate bases in ACGT
                baseProb = prof[base][i] / prof['Z'][i]
                try:
                    innerSum -= baseProb*math.log(baseProb,2)
                except:
                    # Handle log(0) cases when there are no base occurrences or pseudocounts
                    innerSum -= 0
            entropyTerms.append(innerSum)
        # print(sum(entropyTerms))
        return(sum(entropyTerms))
    
    def getConsensus2(self, motifs):
        """Get the consensus sequence from a set of sequences"""
        consensus = list()
        for i in range(self.kmer):
            # First iterate each position for all of the sequences
            baseCountDict = {key: 0 for key in ['A', 'C', 'G', 'T']}
            for j in range(len(motifs)):
                # Then iterate each promoter sequence found
                base = self.sequenceList[j][motifs[j]+i]
                baseCountDict[base] += 1
            invBaseCountDict = {value: key for key, value in baseCountDict.items()}            
            consensus.append(invBaseCountDict[max(baseCountDict.values())])
        return ''.join(consensus)

class FastAreader :
    """Class for parsing fasta sequence files."""
	
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        """Handle fasta data via file object or std input stream"""
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
 
    def readFasta (self):
        """Parse a fasta file
        
        Returns: Tuple with record header at index 0, sequence at index 1
        """
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
			
            header = ''
            sequence = ''
 
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
						
        yield header,sequence

class CommandLine() :
    '''Handle the command line, usage and help requests.
    
    Author: David Bernick
    History: dlb 08/20/2011 Created

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    '''
    
    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        
        Implements a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                              epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                              add_help = True, #default is True 
                                              prefix_chars = '-', 
                                              usage = '%(prog)s [options] -option1[default] <input >output' 
                                             )
        self.parser.add_argument('-g', '--g', action = 'store_true', help='Use Gibbs search instead of random.')       
        self.parser.add_argument('-i', '--i', action = 'store', default=1000, type=int, help='Number of times to find trajectory.')
        self.parser.add_argument('-k', '--k', action = 'store', type=int, default=13, help='Motif size to search for consensus.')
        self.parser.add_argument('-m', '--m', action = 'store_true', help='Print additional info (record headers and motifs).')
        self.parser.add_argument('-p', '--p', action = 'store', default=1, type=float, help='Provide the number of pseudocounts.')
        self.parser.add_argument('-r', '--r', action = 'store_true', help='Use shuffled entropy as control / baseline.')

        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''
    def __init__(self, msg):
        self.msg = msg 

def main(myCommandLine=None):
    '''Main program for finding a Crispr promoter sequence from fasta records and user-defined options'''
    
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else :
        myCommandLine = CommandLine(myCommandLine) # interpret the list passed from the caller of main as the commandline.

    try:
        if myCommandLine.args.i < 1:
            raise Usage('Iterations outside valid range. Try again using -i > 1...')
    except Usage as err:
        print (err.msg)
    try:
        if myCommandLine.args.k < 3:
            raise Usage('Motif length outside valid range. Try again using -k > 2...')
    except Usage as err:
        print (err.msg)
        
    try:
        if myCommandLine.args.p < 0:
            raise Usage('Pseudocounts outside valid range. Try again using -p >= 0...')
    except Usage as err:
        print (err.msg)

    mySeqs = UpstreamSequences(myCommandLine.args.k)
    myFastaReader = FastAreader() # No args for std in
    for head,seq in myFastaReader.readFasta():
        if myCommandLine.args.r:
            shuffledBases = list()
            shuffle = np.random.choice(len(seq), len(seq), False)
            for i in shuffle:
                shuffledBases.append(seq[i])
            shuffledSeq = ''.join(shuffledBases)
            mySeqs.addSequence(head, shuffledSeq)
        else:
            mySeqs.addSequence(head, seq)
    bestResult = 2*mySeqs.kmer # Because 2 is the highest entropy per coloumn w k columns
    promoters = list()
    result, kmers = (float(), list())
    for i in range(myCommandLine.args.i):
        result, kmers = mySeqs.search(myCommandLine.args.p, myCommandLine.args.g, kmers)
        if result < bestResult:
            promoters = kmers
        bestResult = min(result, bestResult)
    if myCommandLine.args.r:
        print('NOTE: Result uses shuffled sequences to serve as a baseline.')
    if myCommandLine.args.m:
        for i in range(len(promoters)):
            print(mySeqs.headerList[i].split()[0])
            pos = promoters[i]
            motif = mySeqs.sequenceList[i][pos: pos+mySeqs.kmer]
            rawSeq = mySeqs.sequenceList[i].lower()
            highlightedMotif = rawSeq.replace(motif.lower(), motif)
            print(highlightedMotif)
    print('Consensus Promoter: ' + str(mySeqs.getConsensus2(promoters)))
    print('Min Entropy: ' + str(bestResult))
    
if __name__ == "__main__":
    main()