#!/usr/bin/env python3
# Name: Trevor Ridgley (tridgley)
# BME-205 Week 2 Assignment - Search for the missing motifs in bacterial genomes
# 
# Note to grader:
# I collaborated with Carmelle Catamura discussing the kScoring calculation which I
#   previously did wrong, instead using binomial mean = np and std dev = sqrt(npq)

"""missingMotif.py - Find motifs that occur less frequently than expected in a DNA sequence

usage: python3 missingMotif.py < a.fasta > out.txt --minMotif 3 --maxMotif 8 --cutoff -5 [--kScoring] 

Input: A fasta file containing 1 or more DNA sequence records.
Output: Report of underrepresented motifs in the DNA sequence, sorted by motif length (descending)
        then by z-score (increasing, hence rarest first)
"""
import sys, scipy.stats, itertools #, time
from collections import defaultdict

class Genome():
    """Create a genome object containing motifs and their statistics."""
    def __init__(self, min, max):
        self.min = min
        self.max = max
        self.motifDict = dict() # defaultdict(int)
#        for i in range(min-2, max+1):
#            self.motifDict[i] = defaultdict(int)
        for i in range(min-2,min):
            self.motifDict[i] = defaultdict(int)
        # Dr. Bernick emailed the following code on Oct 4, 2019
        for i in range(min, max+1):
            self.motifDict[i] = dict()
            for k in itertools.product(['A', 'C', 'G', 'T'], repeat=i):
                kmer = ''.join(k)
                self.motifDict[i][kmer] = 0
        self.seqDict = dict()
        self.genSz = int()
        self.motifProbs = dict() # Actual probabilities
        self.nullProbs = dict() # Null model probabilities
        self.motifMeans = dict() # Null model counts
        self.motifStdDevs = dict()
        self.zScores = dict()
        self.kScores = dict()
#        self.classProbSum = [float() for i in range(self.max+1)]
        for i in range(min, max+1):
            self.motifProbs[i] = defaultdict(int)
            self.motifMeans[i] = dict()
            self.motifStdDevs[i] = dict()
            self.zScores[i] = dict()
            self.kScores[i] = dict()
            self.nullProbs[i] = dict()
            
    def revComp(self, seq):
        """Returns the reverse complement sequence of seq arguments."""
        temp = seq.replace('A', 'Z')
        temp = temp.replace('T', 'A')
        temp = temp.replace('Z', 'T')
        temp = temp.replace('C', 'Z')
        temp = temp.replace('G', 'C')
        temp = temp.replace('Z', 'G')
        return ''.join(temp[::-1])
            
    def findMotifs(self, s):
        """Find motifs from size min to max within each DNA fasta record."""
        self.updateGenomeSize(s)
        # for i in range(len(s) - (self.min - 2) + 1):
        for i in range(len(s) - self.max):
            # Iterate the sequence start index
            for j in range(i + self.min-2, i + self.max + 1):
                # Then iterate the sequence stop index
                motifLen = j - i
                if j <= len(s):
                    # Ensures that motifs at the seq 3' end are not missed
                    m = s[i:j]
                    if set(m) <= {'A', 'C', 'G', 'T'}:
                        # Check that motif is subset of canoncial bases, otherwise invalid motif...
                        self.motifDict[motifLen][m] += 1
                    
    def compressDict(self):
        """Reduce the motif dictionary to forward sequences only.
        
        While reverse-complement keys are deleted, the counts are summed first.
        """
        for i in self.motifDict.keys():
            for motif in list(self.motifDict[i].keys()):
                rc = self.revComp(motif)
                if self.motifDict[i].get(motif) and motif != rc:
                    self.motifDict[i][motif] += self.motifDict[i][self.revComp(motif)]
                    del self.motifDict[i][rc]
                    
    def calcNullProbs(self):
        """Calculate the theoretical probabilities using k-2 Markov null model.
        
        Algo: P(X | k-1, k-2) = P(k1) * P(k2) / P(k3)
              k1, k2 are k-1ers while k3 is k-2er
        """
        
        for i in range(self.min, self.max+1):
            # Iterate each kmer dict
            kMerCount = sum(self.motifDict[i].values()) # same as self.genSz
            for motif in self.motifDict[i].keys():
                # Then calculate the k-1 and k-2 approximations             
                k1 = motif[:-1] if motif[:-1] in self.motifDict[i-1].keys() else self.revComp(motif[:-1])
                k2 = motif[1:] if motif[1:] in self.motifDict[i-1].keys() else self.revComp(motif[1:])
                k3 = motif[1:-1] if motif[1:-1] in self.motifDict[i-2].keys() else self.revComp(motif[1:-1])
                kMerEstimate = self.motifDict[i-1][k1] * self.motifDict[i-1][k2] / self.motifDict[i-2][k3]
                try:
                    # Handle any possible divive by 0 cases
                    self.nullProbs[i][motif] = kMerEstimate / kMerCount
                except:
                    print('calcNullProbs: Divide by 0, setting prob to 0.0')
                    self.nullProbs[i][motif] = 0.0
                
    def calcKScores(self, i):
        """Calculate normalized z-scores for a more reasonable distribution.
        
        normalized stochastic var = pK / pNull
        normalized mean = probability of expected counts within the class k of null model
        normalized std dev = sqrt(N * pNull * (1-pNull))
        ** Note here in the code that normalizedMean = sum of all probilities in class K (roughly n*p)
        """
        normalizedClassDict = dict()
        for motif in self.motifDict[i].keys():
            try:
                normalizedClassDict[motif] = self.motifDict[i][motif] / self.motifMeans[i][motif]
            except:
                normalizedClassDict[motif] = 0.0
        normalizedMean = sum(normalizedClassDict.values()) / len(normalizedClassDict.values())
        # Calculate standard deviation using computing formula
        normalizedStdDev = (sum([(k-normalizedMean)**2 for k in normalizedClassDict.values()]) / 
                            (len(normalizedClassDict.values())-1))**(1/2)
        for motif in normalizedClassDict.keys():
            try:
                self.zScores[i][motif] = (normalizedClassDict[motif] - normalizedMean) / normalizedStdDev
            except:
                self.zScores[i][motif] = 0.0
                
# THIS WAS MY OLD EXTRA CREDIT CODE BEFORE CHATTING WITH CARMEL; THE VALUES WERE MORE
# TIGHTLY DISTRIBUTED AROUND MEAN OF 1
#        for motif in self.motifDict[i].keys():
#            normalizedAct = self.motifDict[i][motif] / self.motifMeans[i][motif]
#            normalizedMean = self.classProbSum[i]
#            normalizedStdDev = (normalizedMean *(1-normalizedMean/self.genSz))**(1/2)
#            try:
#                self.zScores[i][motif] = (normalizedAct - normalizedMean) / normalizedStdDev
#            except:
#                self.zScores[i][motif] = 0.0
                    
    def buildNullModel(self, kScore):
        """Estimate expected number of kmers and their distribution.
        
        Binomial distr mean = np, N = kmer class size (same as genome size), p = prob of kmer
                       std deviation = sqrt(np(1-p))
        """
        self.calcNullProbs()
        for i in range(self.min, self.max+1):
            for motif in self.motifDict[i].keys():
                # Mean = np
                self.motifMeans[i][motif] = self.genSz * self.nullProbs[i][motif]
                # Std dev = sqrt(npq)
                self.motifStdDevs[i][motif] = (self.motifMeans[i][motif] * (1-self.nullProbs[i][motif]))**(1/2)
                # Now calc the z-scores
                if not kScore:
                    try:
                        self.zScores[i][motif] = (self.motifDict[i][motif] - self.motifMeans[i][motif]) / \
                            self.motifStdDevs[i][motif]
                    except:
                        self.zScores[i][motif] = 0.0
            if kScore:
#                self.sumKClassProbs(i)
                self.calcKScores(i)
    
#    def sumKClassProbs(self, i):
#        """Calculate expected probability for each kmer only within SAME CLASS"""
#        self.classProbSum[i] = sum(self.nullProbs[i].values())
                    
    def updateGenomeSize(self, s):
        """Increment the genome size with each new fasta record added to genome.
        
        Note: For simplifying the calculations, N is the genome length - maxMotif.
              Doing this makes the number of motifs all the same.
        """
        self.genSz += len(s) - self.max

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
        
        self.parser.add_argument('--minMotif', action = 'store', default=3, type=int, help='Min length of seq motif to search for')
        self.parser.add_argument('--maxMotif', action = 'store', default=8, type=int, help='Max length of seq motif to search for')
        self.parser.add_argument('--cutoff', action = 'store', type=float, default=0, help='Statistical cutoff for the output report')
        self.parser.add_argument('--kScoring', action = 'store_true', help='Use scipy to calculate z score')
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
    '''Main program for finding and reporting rare kmers with user-defined arguments'''
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else :
        myCommandLine = CommandLine(myCommandLine) # interpret the list passed from the caller of main as the commandline.

    minMotif, maxMotif = (myCommandLine.args.minMotif, myCommandLine.args.maxMotif)

    try:
        # print (myCommandLine.args) # print the parsed argument string .. as there is nothing better to do
        if myCommandLine.args.minMotif < 3 or myCommandLine.args.minMotif > 8:
            raise Usage('minMotif outside valid range, changing to 3...')
    except Usage as err:
        minMotif = 3
        print (err.msg)
       
    try:
        if myCommandLine.args.maxMotif < myCommandLine.args.minMotif or myCommandLine.args.maxMotif > 8:
            raise Usage('maxMotif outside valid range, changing to 8...')
    except Usage as err:
        maxMotif = 8
        print (err.msg)
       
#    start = time.time()
    myGenome = Genome(minMotif, maxMotif)
    myFastaReader = FastAreader() # No args for std in
    for head,seq in myFastaReader.readFasta():
        myGenome.findMotifs(seq)
#    print("Time after finding motifs: {0}".format(time.time() - start))
    myGenome.compressDict()
#    print("N = {0}".format(myGenome.genSz))
    if myCommandLine.args.kScoring:
        print('sequence: reverse\tcount\tExpect\tZscore\tKscore')
    else:
        print('sequence: reverse\tcount\tExpect\tZscore')
    myGenome.buildNullModel(myCommandLine.args.kScoring)
    outputDict = defaultdict(list)
    for i in range(minMotif, maxMotif+1):
        for j in myGenome.motifDict[i].keys():
            # 0 = key, 1 = revcomp, 2 = count, 3 = means, 4 = zScore
            s1 = j if j < myGenome.revComp(j) else myGenome.revComp(j)
            if myCommandLine.args.kScoring and myGenome.zScores[i][j] < myCommandLine.args.cutoff:
                ks = scipy.stats.norm.cdf(myGenome.zScores[i][j])
                outputDict[i].append((s1, myGenome.revComp(s1), myGenome.motifDict[i][j], myGenome.motifMeans[i][j],
                                      myGenome.zScores[i][j], ks))
            elif myGenome.zScores[i][j] < myCommandLine.args.cutoff:
                outputDict[i].append((s1, myGenome.revComp(s1), myGenome.motifDict[i][j], myGenome.motifMeans[i][j],
                                      myGenome.zScores[i][j]))
    for i in range(maxMotif, minMotif - 1, -1):
        if myCommandLine.args.kScoring:
            outputDict[i].sort(key=lambda i : i[4])
            for j in outputDict[i]:
                print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}\t{5:0.2f}'.format(j[0], j[1], j[2], j[3], j[4], j[5]))
        else:
            outputDict[i].sort(key=lambda i : i[4])
            for j in outputDict[i]:
                print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}'.format(j[0], j[1], j[2], j[3], j[4]))
     
#    print("Runtime: {0}".format(time.time() - start))
    
if __name__ == "__main__":
    main()