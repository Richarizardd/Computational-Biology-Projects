from itertools import islice, tee
import numpy
import re

class HMM(object):
    ''' Simple Hidden Markov Model implementation.  User provides
        transition, emission and initial probabilities in dictionaries
        mapping 2-character codes onto floating-point probabilities
        for those table entries.  States and emissions are represented
        with single characters.  Emission symbols comes from a finite.  '''
    
    def __init__(self, A, E, I):
        ''' Initialize the HMM given transition, emission and initial
            probability tables. '''
        
        # put state labels to the set self.Q
        self.Q, self.S = set(), set() # states and symbols
        for a, prob in A.items():
            asrc, adst = a[0], a[1]
            self.Q.add(asrc)
            self.Q.add(adst)
        # add all the symbols to the set self.S
        for e, prob in E.items():
            eq, es = e[0], e[1]
            self.Q.add(eq)
            self.S.add(es)
        
        self.Q = sorted(list(self.Q))
        self.S = sorted(list(self.S))
        
        # create maps from state labels / emission symbols to integers
        # that function as unique IDs
        qmap, smap = {}, {}
        for i in range(len(self.Q)): qmap[self.Q[i]] = i
        for i in range(len(self.S)): smap[self.S[i]] = i
        lenq = len(self.Q)
        
        # create and populate transition probability matrix
        self.A = numpy.zeros(shape=(lenq, lenq), dtype=float)
        for a, prob in A.items():
            asrc, adst = a[0], a[1]
            self.A[qmap[asrc], qmap[adst]] = prob
        # make A stochastic (i.e. make rows add to 1)
        self.A /= self.A.sum(axis=1)[:, numpy.newaxis]
        
        # create and populate emission probability matrix
        self.E = numpy.zeros(shape=(lenq, len(self.S)), dtype=float)
        for e, prob in E.items():
            eq, es = e[0], e[1]
            self.E[qmap[eq], smap[es]] = prob
        # make E stochastic (i.e. make rows add to 1)
        self.E /= self.E.sum(axis=1)[:, numpy.newaxis]
        
        # initial probabilities
        self.I = [ 0.0 ] * len(self.Q)
        for a, prob in I.items():
            self.I[qmap[a]] = prob
        # make I stochastic (i.e. adds to 1)
        self.I = numpy.divide(self.I, sum(self.I))
        
        self.qmap, self.smap = qmap, smap
        
        # Make log-base-2 versions for log-space functions
        self.Alog = numpy.log2(self.A)
        self.Elog = numpy.log2(self.E)
        self.Ilog = numpy.log2(self.I)
    
    def jointProb(self, p, x):
        ''' Return joint probability of path p and emission string x '''
        p = list(map(self.qmap.get, p)) # turn state characters into ids
        x = list(map(self.smap.get, x)) # turn emission characters into ids
        tot = self.I[p[0]] # start with initial probability
        for i in range(1, len(p)):
            tot *= self.A[p[i-1], p[i]] # transition probability
        for i in range(0, len(p)):
            tot *= self.E[p[i], x[i]] # emission probability
        return tot
    
    def jointProbL(self, p, x):
        ''' Return log2 of joint probability of path p and emission
            string x.  Just like self.jointProb(...) but log2 domain. '''
        p = list(map(self.qmap.get, p)) # turn state characters into ids
        x = list(map(self.smap.get, x)) # turn emission characters into ids
        tot = self.Ilog[p[0]] # start with initial probability
        for i in range(1, len(p)):
            tot += self.Alog[p[i-1], p[i]] # transition probability
        for i in range(0, len(p)):
            tot += self.Elog[p[i], x[i]] # emission probability
        return tot
    
    def viterbi(self, x):
        ''' Given sequence of emissions, return the most probable path
            along with its probability. '''
        x = list(map(self.smap.get, x)) # turn emission characters into ids
        nrow, ncol = len(self.Q), len(x)
        mat   = numpy.zeros(shape=(nrow, ncol), dtype=float) # prob
        matTb = numpy.zeros(shape=(nrow, ncol), dtype=int)   # backtrace
        # Fill in first column
        for i in range(0, nrow):
            mat[i, 0] = self.E[i, x[0]] * self.I[i]
        # Fill in rest of prob and Tb tables
        for j in range(1, ncol):
            for i in range(0, nrow):
                ep = self.E[i, x[j]]
                mx, mxi = mat[0, j-1] * self.A[0, i] * ep, 0
                for i2 in range(1, nrow):
                    pr = mat[i2, j-1] * self.A[i2, i] * ep
                    if pr > mx:
                        mx, mxi = pr, i2
                mat[i, j], matTb[i, j] = mx, mxi
        # Find final state with maximal probability
        omx, omxi = mat[0, ncol-1], 0
        for i in range(1, nrow):
            if mat[i, ncol-1] > omx:
                omx, omxi = mat[i, ncol-1], i
        # Backtrace
        i, p = omxi, [omxi]
        for j in range(ncol-1, 0, -1):
            i = matTb[i, j]
            p.append(i)
        p = ''.join(map(lambda x: self.Q[x], p[::-1]))
        return omx, p # Return probability and path
    
    def viterbiL(self, x):
        ''' Given sequence of emissions, return the most probable path
            along with log2 of its probability.  Just like viterbi(...)
            but in log2 domain. '''
        x = list(map(self.smap.get, x)) # turn emission characters into ids
        nrow, ncol = len(self.Q), len(x)
        mat   = numpy.zeros(shape=(nrow, ncol), dtype=float) # prob
        matTb = numpy.zeros(shape=(nrow, ncol), dtype=int)   # backtrace
        # Fill in first column
        for i in range(0, nrow):
            mat[i, 0] = self.Elog[i, x[0]] + self.Ilog[i]
        # Fill in rest of log prob and Tb tables
        for j in range(1, ncol):
            for i in range(0, nrow):
                ep = self.Elog[i, x[j]]
                mx, mxi = mat[0, j-1] + self.Alog[0, i] + ep, 0
                for i2 in range(1, nrow):
                    pr = mat[i2, j-1] + self.Alog[i2, i] + ep
                    if pr > mx:
                        mx, mxi = pr, i2
                mat[i, j], matTb[i, j] = mx, mxi
        # Find final state with maximal log probability
        omx, omxi = mat[0, ncol-1], 0
        for i in range(1, nrow):
            if mat[i, ncol-1] > omx:
                omx, omxi = mat[i, ncol-1], i
        # Backtrace
        i, p = omxi, [omxi]
        for j in range(ncol-1, 0, -1):
            i = matTb[i, j]
            p.append(i)
        p = ''.join(map(lambda x: self.Q[x], p[::-1]))
        return omx, p # Return log probability and path
        
def pairwise(iterable):
    """ Create iterator over adjacent pairs of elements in given
        iterator. """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def read_fasta_char_by_char(fn):
    """ Read FASTA file, yielding nucleotides one by one """
    with open(fn) as fh:
        for ln in fh:
            if ln[0] == '>':
                continue
            for c in ln.strip():
                yield c

def read_labels_char_by_char(fn):
    """ Read label file, yielding labels one by one """
    with open(fn) as fh:
        for ln in fh:
            for c in ln.strip():
                yield c

def read_genome_and_training(franken_fn, training_fn):
    """ Given filename for Frankengenome and filename for training data,
        return iterator over labeled positions of the Frankengenome,
        yielding pairs like ('A', 0) """
    return zip(read_fasta_char_by_char(franken_fn),
               read_labels_char_by_char(training_fn))

def read_training_pairs(training_fn):
    """ Return iterator over adjacent label pairs in given label file """
    return pairwise(read_labels_char_by_char(training_fn))


def evaluate_on_test_data(prediction, test_data_fn):
    """ Given an iterator over all the predictions (one per
        genome position), and given the name of the test data
        file, return a pair giving the number of correct and
        the total number of predictions. """
    ncorrect, ntot = 0, 0
    for pred, act in zip(islice(prediction, 50000, None),  # skip first 50K predictions
                         read_labels_char_by_char(test_data_fn)):
        if pred == act:
            ncorrect += 1
        ntot += 1
    return ncorrect, ntot
    
def main():
    # Reads in training set and testing set
    train_fn = "trainingData1.txt"; test_fn = "testData1.txt"; franken_fn = "frankengene1.fasta"
    genome_ziptrain = read_genome_and_training(franken_fn, train_fn)
    state_train = open(train_fn).read().rstrip()
    _0 = float(state_train.count('0')); _1 = float(state_train.count('1'))
    
    # Building transition states
    _00 = float(len([m.start() for m in re.finditer('(?=00)', state_train)]))/_0
    _01 = float(len([m.start() for m in re.finditer('(?=01)', state_train)]))/_0
    _10 = float(len([m.start() for m in re.finditer('(?=10)', state_train)]))/_1
    _11 = float(len([m.start() for m in re.finditer('(?=11)', state_train)]))/_1
    transition_states = {"00":_00, "10":_10, "01":_01, "11":_11}
    
    # Building transition nucleotide
    genome0 = ""; genome1 = ""
    for base, state in genome_ziptrain:
        if state == '0':
            genome0 += base
        if state == '1':
            genome1 += base
            
    _0A = float(genome0.count('A'))/_0; _1A = float(genome1.count('A'))/_1
    _0C = float(genome0.count('C'))/_0; _1C = float(genome1.count('C'))/_1
    _0G = float(genome0.count('G'))/_0; _1G = float(genome1.count('G'))/_1
    _0T = float(genome0.count('T'))/_0; _1T = float(genome1.count('T'))/_1
    transition_nucleotides = {"0A":_0A, "0C":_0C, "0G":_0G, "0T":_0T,
                              "1A":_1A, "1C":_1C, "1G":_1G, "1T":_1T}
                          
    # Builds initial case
    transition_init = {"0":0.50, "1":0.50}
    
    # Builds model
    franken_hmm = HMM(transition_states, transition_nucleotides, transition_init)
    
    # Prediction
    genome = read_fasta_char_by_char(franken_fn)
    prediction = franken_hmm.viterbiL(genome)
    outfile = open("prediction.txt", "w")
    outfile.write(prediction[1])
    outfile.close()
    #ncorrect, ntot = evaluate_on_test_data(prediction[1], test_fn)
    #print "% Correct: ", float(ncorrect)/float(ntot)*100,"%"
main()