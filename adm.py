#!/usr/bin/env python3

import sys
import getopt
import base
import numpy
import io
import random
import numpy as np
import copy

to_int = dict(list(zip("ACGT", list(range(4)))))
nucs="ACGT"

myf="%.5f"

pseudo_count=1
verbose=False

class adm(object):

    def __init__(self, t_probs, i_probs=None):
        if i_probs is None:
            if type(t_probs) == np.ndarray:
                rows,cols=t_probs.shape                                 # When called as adm(transitions)
                self.k = cols
                t = normalize_transition_matrix(t_probs)
                self.transition_probabilities = t[:,1:]                 # Note: initial probabilities are not included here
                self.initial_probabilities = generate_all_initial_probabilities2(t)
            else:
                k = t_probs                                             # When called as adm(k)
                self.k = k
                self.transition_probabilities = numpy.empty((16, k-1))
                self.initial_probabilities = numpy.empty((4, k))
                self.transition_probabilities.fill(0.25)
                self.initial_probabilities.fill(0.25)
        else:
            assert numpy.isfinite(t_probs).all()
            assert numpy.isfinite(i_probs).all()
            self.k = len(i_probs[0])                 # When called as adm(t_probs, i_probs)
            assert self.k - 1 == len(t_probs[0])
            assert t_probs.shape == (16,self.k-1)
            assert i_probs.shape == (4, self.k)
            self.transition_probabilities = t_probs
            self.initial_probabilities = i_probs
        self.shape = (16, self.k)
        
    def invariantxx(self):
        assert self.transition_probabilities.shape == (16, self.k-1)
        assert self.initial_probabilities.shape == (4, self.k)

        for i in range(self.k):
            assert abs(sum(self.initial_probabilities[:, i]) - 1.0) < 0.00001
            
        for i in range(self.k-1):
            for a in range(4):
                assert abs(sum(self.transition_probabilities[4*a:4*(a+1), i]) - 1.0) < 0.00001

    def invariant(self):
        if self.transition_probabilities.shape != (16, self.k-1):
            return False
        if self.initial_probabilities.shape != (4, self.k):
            return False

        for i in range(self.k):
            if abs(sum(self.initial_probabilities[:, i]) - 1.0) >= 0.00001:
                return False
            
        for i in range(self.k-1):
            for a in range(4):
                s = sum(self.transition_probabilities[4*a:4*(a+1), i])
                if abs(s - 1.0) >= 0.00001 and s != 0.0:
                    return False

        # If we arrive at a state, we must also be able to leave it
        for i in range(self.k-1):
            for a in range(4):
                if sum(self.transition_probabilities[4*a:4*(a+1), i]) == 0.0 and self.initial_probabilities[a, i] > 0.0:
                    return False

        if self.initial_probabilities.max() > 1.0 or self.initial_probabilities.min() < 0.0:
            return False
            
        if self.transition_probabilities.max() > 1.0 or self.transition_probabilities.min() < 0.0:
            return False
            
        return True

    def representation(self):
        result = np.zeros((16,self.k))
        result[0:4,0] = self.initial_probabilities[:,0]
        result[:, 1:] = self.transition_probabilities
        return result
                
    # Get the initial probability of 'a'
    def get_initial(self, a):
        if type(a) is str:
            a = to_int[a]
        return self.initial_probabilities[a][0]

    def set_initial(self, a, p):
        if type(a) is str:
            a = to_int[a]
        self.initial_probabilities[a][0] = p

    def get_transition(self, pos, a, b):
        if type(a) is str:
            a = to_int[a]
        if type(b) is str:
            b = to_int[b]
        row = 4*a + b
        return self.transition_probabilities[row][pos]

    def set_transition(self, pos, a, b, p):
        if type(a) is str:
            a = to_int[a]
        if type(b) is str:
            b = to_int[b]
        row = 4*a + b
        self.transition_probabilities[row][pos] = p

    def __str__(self):
        result=[]
        temp1=io.StringIO()
        temp2=io.StringIO()
        for i in range(4):
            base.printmatrix(self.transition_probabilities[i*4:(i+1)*4,:], temp1, format=myf)
            temp1.write("\n")
        base.printmatrix(self.initial_probabilities, temp2, format=myf)
        result.append("Model width is %i" % self.k)
        result.append("Transition probabilities are")
        result.append(temp1.getvalue())
        result.append("Initial probabilities are:")
        result.append(temp2.getvalue())
        temp1.close()
        temp2.close()
        return "\n".join(result)

    def str2(self, fmt="%.6f"):
        temp=io.StringIO()
        for row in range(16):
            for col in range(self.k-1):
                temp.write("%.6f\t" % self.transition_probabilities[row,col])
            temp.write("ADM_DI\t%s\n" % (nucs[row//4] + nucs[row%4]))
        for row in range(4):
            for col in range(self.k):
                temp.write("%.6f\t" % self.initial_probabilities[row,col])
            temp.write("ADM_MONO_%s\n" % nucs[row])
        s=temp.getvalue()
        temp.close()
        return s

transform = [15, 11, 7, 3, 14, 10, 6, 2,
             13, 9, 5, 1, 12, 8, 4, 0]

def reverse_complement_adm(adm1):
    k = adm1.k
    i=base.reverse_complement_pwm(adm1.initial_probabilities)
    t = numpy.zeros((16, k-1))
    for j in range(k-1):
        for ab in range(16):
            a = ab // 4
            b = ab % 4
            divisor = adm1.initial_probabilities[b, j+1]
            if divisor != 0.0:
                t[transform[ab], k-j-2] = adm1.transition_probabilities[ab, j] * adm1.initial_probabilities[a, j] / divisor
            else:
                t[transform[ab], k-j-2] = 0.0
    return adm(t,i)

# Get all nucleotide k-mers
def get_kmers(k):
    assert k > 0
    result=[]
    A=['A']*k
    v=[0]*k
    v[k-1] = -1
    for current_string in range(pow(4,k)):
        i = k-1
        while v[i] == 3:
            v[i]=0
            A[i] = nucs[v[i]]
            i -= 1
        v[i] += 1
        A[i] = nucs[v[i]]
        result.append("".join(A))
      
    return result


def adm_probability(a, s):
    k=len(s)
    assert a.k >= k
    current=to_int[s[0]]
    prob = a.initial_probabilities[current, 0]
    for i in range(k-1):
        next=to_int[s[i+1]]
        prob *= a.transition_probabilities[4*current+next, i]
        current=next
    return prob

def product_adm_probability(a1, a2, s):
    k=len(s)
    assert a1.k >= k
    assert a2.k >= k
    current=to_int[s[0]]
    prob = a1.initial_probabilities[current, 0] * a2.initial_probabilities[current, 0]
    for i in range(k-1):
        next=to_int[s[i+1]]
        prob *= a1.transition_probabilities[4*current+next, i] * a2.transition_probabilities[4*current+next, i]

        current=next
    return prob

def probability_of_equal_paths(a1, a2, k):
    kmers=get_kmers(k)
    prob = 0.0
    for s in kmers:
        prob += product_adm_probability(a1, a2, s)
    return prob

def conditional_product_probability(a1, a2, s, k):
    assert len(s) <= k
    return product_adm_probability(a1, a2, s) / probability_of_equal_paths(a1, a2, k)


class equal_future(object):

    def helper(self, h, c, d):
        assert c in nucs
        assert d in nucs
        assert h < self.k-1
        return self.a1.transition(h, c, d) * self.a2.transition(h, c, d)

    def __init__(self, a1, a2):

        self.a1 = a1
        self.a2 = a2
        self.k = k = a1.k
        self.t = numpy.empty((4, k))
        for i in range(k-1):
            length = k - i - 1
            kmers = get_kmers(length)
            for a in nucs:
                for s in kmers:
                    prob = 1.0
                    prev=a
                    j = i
                    for x in s:
                        prob *= self.helper(j, prev, x)
                        prev = x
                        j += 1
                    self.t[to_int[a], i] += prob
        for a in nucs:
            self.t[to_int[a], k-1] = 1.0
            
    def get(self, i, c):
        return self.t[c, i]
        
def force_adms_equal(a1, a2):
    k = a1.k
    assert k == a2.k
    
    trans = numpy.empty((16, k-1))
    init = numpy.empty((4, k))
    temp = [0.0]*4
    ef = equal_future(a1, a2)
    for i in range(4):
        temp[i] = a1.initial_probabilities[i,0] * a2.initial_probabilities[i,0]
        temp[i] *= ef.get(0, i)
    init[:,0] = normalize(temp)
    for j in range(k-1):
        for a in range(4):
            temp = a1.transition_probabilities[4*a:4*(a+1),j] * a2.transition_probabilities[4*a:4*(a+1),j]
            for b in range(4):
                temp[b] = temp[b] * ef.get(j+1, b) / ef.get(j, a)
            #temp = normalize(temp)
            trans[4*a:4*(a+1),j] = temp

    a=adm(trans, init)

    return generate_all_initial_probabilities(a)


# The list 'l' contains sequences, the key'th element is these sequences are used for comparison    
def my_argmax(l, key):
    max_arg=0
    max_value=l[0][key]
    for i, x in enumerate(l):
        if x[key] > max_value:
            max_value = x[key]
            max_arg = i
    return max_arg


def underline(s):
    print(s)
    print("="*len(s))



def generate_all_initial_probabilities(myadm):
    k = myadm.k
    init=numpy.zeros((4,k))
    for a in range(4):
        init[a, 0] = myadm.initial_probabilities[a,0]
    for i in range(1, k):
        for b in range(4):
            for a in range(4):
                init[b, i] += init[a, i-1] * myadm.get_transition(i-1, a, b)
    return adm(myadm.transition_probabilities, init)

def generate_all_initial_probabilities2(transition):   # This version takes as input the normalized 16xk sized transition matrix
    rows,cols = transition.shape
    k=cols
    init=numpy.zeros((4,k))
    for a in range(4):
        init[a, 0] = transition[a,0]
    for i in range(1, k):
        for b in range(4):
            for a in range(4):
                init[b, i] += init[a, i-1] * transition[4*a+b, i]
    return init

def max_string_for_adm(myadm):
    # The state is a list of size four with each element a pair of previous char and current probability
    current_state = [ (0, myadm.get_initial(x)) for x in "ACGT" ]   # 0 is a dummy because no previous char exists
    path = []
    k = myadm.k
    for current_pos in range(k-1):
        new_state = []
        for j in range(4):
            temp=[ current_state[i][1] * myadm.get_transition(current_pos, i, j) for i in range(4) ]
            c=numpy.argmax(temp)
            p=max(temp)
            new_state.append((c, p))
        i = my_argmax(new_state, 1)
        path.append(nucs[new_state[i][0]])  # Previous char on the path to the winner
        current_state=new_state
    i=my_argmax(new_state, 1)
    path.append(nucs[i])
    return "".join(path), new_state[i][1]  # Return the maximum path and its probability

def read_adm_from_count_file(filename, pseudo_count=0.0):
    nucs="ACGT"
    with open(filename, "r") as f:
        lines = f.readlines()
    split_lines=[x.split() for x in lines]
    a=np.array(split_lines)
    a = a.astype(float) + pseudo_count
    return adm(a)

def read_adm_from_list_of_lines(lines):
    nucs="ACGT"
    if verbose:
        for line in lines:
            print(line.rstrip("\n"))
    #print len(lines)
    assert len(lines) == 20
    split_lines=[x.split() for x in lines]
    for i in range(16):
        assert split_lines[i][-2] == "ADM_DI"
        assert split_lines[i][-1] == nucs[i//4] + nucs[i%4]
    for i in range(16, 20):
        assert split_lines[i][-1] == "ADM_MONO_%s" %  nucs[i-16]

    t_probs = []
    for i in range(16):
        t_probs.append(list(map(float, split_lines[i][:-2])))

    i_probs = []
    for i in range(16,20):
        i_probs.append(list(map(float, split_lines[i][:-1])))

    return adm(numpy.array(t_probs), numpy.array(i_probs))

def read_adm_from_file(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    return read_adm_from_list_of_lines(lines)

def normalize(v):
    return [ x / sum(v) for x in v ]

def disturbe_adm(myadm, seed):
    s=[to_int[x] for x in seed]
    k=myadm.k
    assert k == len(seed)
    newadm = adm(k)
    new_init = [ myadm.initial(i)*myadm.transition(0, i, seed[1]) for i in range(4) ]
    new_init = normalize(new_init)
    for i in range(4):
        newadm.set_initial(i, new_init[i])
    for pos in range(k-2):
        for i in range(4):
            new_trans = [ myadm.transition(pos, i, j)*myadm.transition(pos+1, j, seed[pos+2]) for j in range(4) ]
            new_trans = normalize(new_trans)
            for j in range(4):
                newadm.set_transition(pos, i, j, new_init[j])
    for i in range(4):
        for j in range(4):
            newadm.set_transition(k-2, i, j, myadm.transition(k-2, i, j))
    return newadm

# Maximum elementwise distance between corresponding initial and transition matrices
def adm_distance(a1, a2):
#    if verbose:
#        print "Differences between two models:"
    x=abs(a1.initial_probabilities[:,0] - a2.initial_probabilities[:,0])
#    if verbose:
        #print x
        #    base.printmatrix(x)
        #print numpy.max(x)

    y=abs(a1.transition_probabilities - a2.transition_probabilities)
#    if verbose:
        #base.printmatrix(y)
        #print numpy.max(y)

    all_deviations = numpy.concatenate((x.flatten(), y.flatten()))
#    if verbose:
#        print "Average:", numpy.average(all_deviations)
#        print "Std:", numpy.std(all_deviations)
#    return max(all_deviations), numpy.average(all_deviations)
    return max(numpy.max(x), numpy.max(y))

class myaccumulate(object):
    def __init__(self, f, init):
        self.f = f
        self.v = init

    def get(self):
        return self.v

    def add(self, x):
        self.v = self.f(self.v, x)
        return self.v
    
# Maximum elementwise distance between corresponding initial and transition matrices
# This version doesn't compare transition probabilities if the corresponding initial
# probability is zero.
def adm_distance2(a1, a2):
    assert a1.k == a2.k
    k=a1.k
    x=abs(a1.initial_probabilities[:,0] - a2.initial_probabilities[:,0])
#    y=abs(a1.transition_probabilities - a2.transition_probabilities)
#    base.printmatrix(y, format=myf)
    #max_dist=0.0
    ac = myaccumulate(max, 0)
    for i in range(4):
        for j in range(k-1):
            if a1.initial_probabilities[i, j] <= 0.01 or a2.initial_probabilities[i, j] <= 0.01:
                continue
            for l in range(4):
                ac.add(abs(a1.transition_probabilities[4*i+l, j] - a2.transition_probabilities[4*i+l, j]))
                
    return max(numpy.max(x), ac.get())

# Maximum elementwise distance between corresponding initial and transition matrices
# This version weights the distance between transition probabilities by the corresponding initial
# probabilities.

def weighted_adm_distance(a1, a2):
    assert a1.k == a2.k
    k=a1.k
    x=abs(a1.initial_probabilities[:,0] - a2.initial_probabilities[:,0])
#    y=abs(a1.transition_probabilities - a2.transition_probabilities)
#    base.printmatrix(y, format=myf)
    #max_dist=0.0
    ac = myaccumulate(max, 0)
    for i in range(4):        # Compares probabilities of moving from i to l when in position j
        for j in range(k-1):
            p1 = a1.initial_probabilities[i, j]
            p2 = a2.initial_probabilities[i, j]
            for l in range(4):
                ac.add(abs(p1*a1.transition_probabilities[4*i+l, j] - p2*a2.transition_probabilities[4*i+l, j]))
                
    return max(numpy.max(x), ac.get())


def random_distribution():
    a=random.uniform(0,1)
    b=random.uniform(0,1)
    c=random.uniform(0,1)
    l=[a, b, c]
    l=sorted(l)
    a,b,c = l
    return [a, b-a, c-b, 1-c]

def random_adm(k):
    init=np.zeros((4,k))
    trans=np.zeros((16,k-1))
    init[:, 0] = random_distribution()
    for i in range(k-1):
        for a in range(4):
            trans[4*a:4*(a+1), i] = random_distribution()
    a=adm(trans, init)
    a=generate_all_initial_probabilities(a)
    return a

# use numpy.random.choice instead
def random_value(p):
    cp=np.cumsum(p)
    x=random.uniform(0,1)
    for i in range(len(p)):
        if x <= cp[i]:
            return i
    return len(p)-1



def generate(n, a):
    "Generate 'n' sequences using adm model 'a'"
    
    k=a.k
    data=[]
    for i in range(n):
        c=random_value(a.initial_probabilities[:,0])
        seq=[nucs[c]]
        for j in range(k-1):
            cnew=random_value(a.transition_probabilities[4*c:4*(c+1), j])
            seq.append(nucs[cnew])
            c=cnew
        data.append("".join(seq))
    return data



def generate2(n, a, seed):
    """Generated sequences using adm model 'a' until we have accumulated 'n' sequences
that are within Hamming distance 2 from the seed"""
    
    k=a.k
    data=[]
    i = 0
    while len(data) < n:
        hd=0
        c=random_value(a.initial_probabilities[:,0])
        seq=[nucs[c]]
        if seq[0] != seed[0]:
            hd += 1
        for j in range(k-1):
            cnew=random_value(a.transition_probabilities[4*c:4*(c+1), j])
            seq.append(nucs[cnew])
            if seq[j+1] != seed[j+1]:
                hd += 1
                if hd > 2:
                    break
            c=cnew
        if hd <= 2:
            data.append("".join(seq))
    return data



def generate3(n, a, seed):
    """Generated sequences using adm model 'a' until we have accumulated 'n' sequences
that are within Hamming distance 2 from the seed. If the distance is two, then
the mismatches must be consequent."""
    
    k=a.k
    data=[]
    i = 0
    while len(data) < n:
        hd=0
        c=random_value(a.initial_probabilities[:,0])
        seq=[nucs[c]]
        error_in_previous_pos=False
        if seq[0] != seed[0]:
            hd += 1
            error_in_previous_pos=True
        for j in range(k-1):
            cnew=random_value(a.transition_probabilities[4*c:4*(c+1), j])
            seq.append(nucs[cnew])
            if seq[j+1] != seed[j+1]:
                if hd==1 and not error_in_previous_pos:
                    hd=3
                    break
                hd += 1
                error_in_previous_pos=True
                if hd > 2:
                    break
            else:
                error_in_previous_pos=False
            c=cnew
        if hd <= 2:
            data.append("".join(seq))
    return data



def normalize_transition_matrix(t):
    rows, cols = t.shape
    assert rows == 16
    for i in range(cols):
        for j in range(4):
            s=sum(t[4*j:4*(j+1), i])
            if s > 0.0:
                for l in range(4):
                    t[4*j+l, i] /= s
    return t



def is_match(s, pattern):
    k=len(s)
    assert k == len(pattern)
    for i in range(k):
        if s[i] != pattern[i] and pattern[i] != 'N':
            return False
    return True



def get_mismatches(s, pattern):
    k=len(s)
    assert k == len(pattern)
    mismatches=[]
    for i in range(k):
        if s[i] != pattern[i]:
            mismatches.append(i)
    return mismatches



def learn_adm(data):
    k=len(data[0])
    init=np.zeros((4,k))
    trans=np.zeros((16,k-1))
    for seq in data:
        for i, c in enumerate(seq):
            init[to_int[c], i] += 1
        for i in range(k-1):
            trans[to_int[seq[i]]*4 + to_int[seq[i+1]], i] += 1

        
    init += pseudo_count  # Add pseudo count
    init = base.normalize(init)
    trans += pseudo_count  # Add pseudo count
    trans = normalize_transition_matrix(trans)
    learned = adm(trans, init)
    if verbose:
        print(learned)
    return learned



def learn_multinomial_adm(data, seed):
    k=len(data[0])
    init=np.zeros((4,k))
    trans=np.zeros((16,k-1))
    number_of_sites=0
    for seq in data:
        mm = get_mismatches(seq, seed)
        if len(mm) > 2:
            continue
        if len(mm) == 2 and mm[0] + 1 != mm[1]:  # Not consecutive
            continue
        number_of_sites += 1
        if len(mm) == 0:
            for i in range(k):
                init[to_int[seq[i]], i] += 1
            for i in range(k-1):
                trans[to_int[seq[i]]*4 + to_int[seq[i+1]], i] += 1
        elif len(mm) == 2:
            i = mm[0]
#            init[to_int[seq[i]], i] += 1
#            init[to_int[seq[i+1]], i+1] += 1
            trans[to_int[seq[i]]*4 + to_int[seq[i+1]], i] += 1
        elif len(mm) == 1:
            i = mm[0]
            init[to_int[seq[i]], i] += 1
            if i == 0:
#                init[to_int[seq[i+1]], i+1] += 1
                trans[to_int[seq[i]]*4 + to_int[seq[i+1]], i] += 1
            elif i == k-1:
#                init[to_int[seq[i-1]], i-1] += 1
                trans[to_int[seq[i-1]]*4 + to_int[seq[i]], i-1] += 1
            else:
#                init[to_int[seq[i-1]], i-1] += 1
#                init[to_int[seq[i+1]], i+1] += 1
                trans[to_int[seq[i-1]]*4 + to_int[seq[i]], i-1] += 1
                trans[to_int[seq[i]]*4 + to_int[seq[i+1]], i] += 1
                
            
    init += pseudo_count  # Add pseudo count
    init = base.normalize(init)
    trans += pseudo_count  # Add pseudo count
    trans = normalize_transition_matrix(trans)
    learned = adm(trans, init)
    if verbose:
        print(learned)
        print("The number of sites used: %i" % number_of_sites)
    return learned, number_of_sites



def correct_multinomial_adm(myadm, seed):
    t=copy.deepcopy(myadm.transition_probabilities)
    init=copy.deepcopy(myadm.initial_probabilities)
    k = myadm.k
    for i in range(k-3, -1, -1):
        s = to_int[seed[i+2]]
        for a in range(4):   # Starting character
            temp=[0.0] * 4
            for b in range(4): # Ending character
                temp[b] = t[4*a+b, i] / t[4*b+s,i+1]
            temp=normalize(temp)
            for b in range(4):
                t[4*a+b, i] = temp[b]

    # Normalize initial probabilities                
    temp=[0.0] * 4
    s = to_int[seed[1]]
    for b in range(4): # Ending character
        temp[b] = init[b, 0] / t[4*b+s, 0]
    temp=normalize(temp)
    for b in range(4):
        init[b, 0] = temp[b]

    myadm2=adm(t, init)
    init3=generate_all_initial_probabilities(myadm2)
#    return adm(t, copy.deepcopy(myadm.initial_probabilities))
    return adm(t, init3)

def is_consistent(a1):
    a2=generate_all_initial_probabilities(a1)
    x=abs(a1.initial_probabilities - a2.initial_probabilities)
    y=abs(a1.transition_probabilities - a2.transition_probabilities)
    dist = max(numpy.max(x), numpy.max(y))
    print("distance is %f" % dist)
    return dist < 0.01

usage=""

if __name__ == "__main__":   # Are we importing as module
    optlist, args = getopt.getopt(sys.argv[1:], 'vrsc')
    optdict=dict(optlist)
    args = [sys.argv[0]] + args

    try:
        optdict['-v']
        verbose=True
    except KeyError:
        verbose=False

    try:
        optdict['-c']    # Input matrix contains counts of dinucleotide probabilities (of size 16xk)
        counts=True
    except KeyError:
        counts=False

    pwmfilename=args[1]
    if counts:
        orig = read_adm_from_count_file(pwmfilename)
    else:
        orig = read_adm_from_file(pwmfilename)
    recreated_adm=generate_all_initial_probabilities(orig)
    
    try:
        optdict['-r']
        print(recreated_adm.str2(), end='')
        sys.exit(0)
    except KeyError:
        pass

    try:
        optdict['-s']
        if is_consistent(orig):
            print("is consistent")
        else:
            print("not consistent")
        sys.exit(0)
    except KeyError:
        pass
    
    print("Model length is %i" % orig.k)
    print(orig)

    # a=disturbe_adm(myadm, path)
    # if verbose:
    #     print "Disturbed ADM"
    #     print a

    # d,ad = adm_distance(myadm, a)
    # if verbose:
    #     print "Distance between original and disturbed model is %f" % d
    # print d, ad

    print("Recreated ADM")
    print(recreated_adm.str2())
    if verbose:
        underline("Recreate initial probabilities:")
        base.printmatrix(recreated_adm.initial_probabilities, format=myf)
        print("Distance between original and recreated original: %.4f" % adm_distance(recreated_adm, orig))

    path, prob =  max_string_for_adm(recreated_adm)
    print("Path %s gives maximum probability %f" %(path, prob))

    sys.exit(0)

    k = recreated_adm.shape[1]
    kmers = get_kmers(k)
    pairs = [ (kmer, adm_probability(recreated_adm, kmer)) for kmer in kmers]
    pairs.sort(key = lambda kmer_prob : kmer_prob[1])
    for kmer, prob in pairs:
        print("%s\t%f" % (kmer, prob))
        
    count=int(args[2])


    #data=generate(count, recreated_adm)
    #data=generate2(count, recreated_adm, path)
    data=generate3(count, recreated_adm, path)

    underline("\nNormal learning method:")
    a=learn_adm(data)
    #print "Distance between recreated original and traditional aligned: %.4f %.4f" % adm_distance2(recreated_adm, a)
    #print "Weighted distance between recreated original and traditional aligned: %.4f %.4f" % weighted_adm_distance(recreated_adm, a)
    print("Distance between recreated original and traditional aligned: %.4f" % max(adm_distance2(recreated_adm, a)))
    print("Weighted distance between recreated original and traditional aligned: %.4f" % max(weighted_adm_distance(recreated_adm, a)))

    underline("\nMultinomial learning method:")
    b,sites=learn_multinomial_adm(data, path)
    print("Number of sites used: %i" % sites)
    #print "Distance between recreated original and multinomial learned: %.4f %.4f" % adm_distance2(recreated_adm, b)
    #print "Weighted distance between recreated original and multinomial learned: %.4f %.4f" % weighted_adm_distance(recreated_adm, b)
    print("Distance between recreated original and multinomial learned: %.4f" % max(adm_distance2(recreated_adm, b)))
    print("Weighted distance between recreated original and multinomial learned: %.4f" % max(weighted_adm_distance(recreated_adm, b)))

    underline("\nCorrected multinomial learning method:")
    c=correct_multinomial_adm(b, path)
    print()
    #print "Distance between recreated original and corrected multinomial learned: %.4f %.4f" % adm_distance2(recreated_adm, c)
    #print "Weighted distance between recreated original and corrected multinomial learned: %.4f %.4f" % weighted_adm_distance(recreated_adm, c)
    print("Distance between recreated original and corrected multinomial learned: %.4f" % max(adm_distance2(recreated_adm, c)))
    print("Weighted distance between recreated original and corrected multinomial learned: %.4f" % max(weighted_adm_distance(recreated_adm, c)))


    if verbose:
        print("Orig:\n", orig)
        print("Recreated:\n", recreated_adm)
        print("b:\n", b)
        print("c:\n", c)
