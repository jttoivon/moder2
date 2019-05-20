#!/usr/bin/env python3

import sys
import io
import string
import numpy
import math
import tempfile
import subprocess
import re
import subprocess

def hms_to_seconds(s):
    hours_given = s.count(":") == 2
    if hours_given:
        hours, minutes, seconds = list(map(float, s.split(":")))
    else:
        minutes, seconds = list(map(float, s.split(":")))
        hours=0
    return hours*3600 + minutes*60 + seconds

iupac_complement=dict(list(zip(
    "ACGTRYSWKMBDHVN.acgtryswkmbdhvn",
    "TGCAYRWSMKVHDBN.tgcayrwsmkvhdbn")))

iupac_codes={
    "A":"A",
    "C":"C",
    "G":"G",
    "T":"T",
    "R":"AG",
    "Y":"CT",
    "S":"CG",
    "W":"AT",
    "K":"GT",
    "M":"AC",
    "B":"CGT",
    "D":"AGT",
    "H":"ACT",
    "V":"ACG",
    "N":"ACGT",
#    ".":"ACGT"
}

iupac_probabilities = {
    "A":[1.0, 0, 0, 0],
    "C":[0, 1.0, 0, 0],
    "G":[0, 0, 1.0, 0],
    "T":[0, 0, 0, 1.0],
    "R":[0.5, 0, 0.5, 0],
    "Y":[0, 0.5, 0, 0.5],
    "S":[0, 0.5, 0.5, 0],
    "W":[0.5, 0, 0, 0.5],
    "K":[0, 0, 0.5, 0.5],
    "M":[0.5, 0.5, 0, 0],
    "B":[0, 1.0/3, 1.0/3, 1.0/3],
    "D":[1.0/3, 0, 1.0/3, 1.0/3],
    "H":[1.0/3, 1.0/3, 0, 1.0/3],
    "V":[1.0/3, 1.0/3, 1.0/3, 0],
    "N":[0.25, 0.25, 0.25, 0.25],
}

to_int = dict(list(zip("ACGT", list(range(4)))))

def dna_to_number(s):
    code = 0
    for c in s:
        code = (code << 2) + to_int[c]
    return code

def number_to_dna(code, k):
    s=[]
    nucs="ACGT"
    mask = 3
    for i in range(k):
        s.append(nucs[(code >> (i*2)) & mask])
    s.reverse()
    return "".join(s)
        
# "str" is the starting tag.
# pos is relative to the starting tag
# return count lines
def find_lines(x, str, pos, count):
    resultfile= io.StringIO(string.join(x,""))
    resultlist=[]
    # read matrix header
    while True:
        line=resultfile.readline()
        if re.match(str, line):
            #line = resultfile.readline()
            if pos == 0:
                resultlist.append(line)
                count -= 1
            else:
                pos -= 1
            break
        elif line=="":
            raise AttributeError("Tag '%s' not found" % str)
    while pos > 0 and line != "":
        line = resultfile.readline()
        pos -= 1
    while count > 0:
        count -= 1
        line = resultfile.readline()
        resultlist.append(line)
    return resultlist

# x is a list of lines
# The first line defines the dimensions: e.g. 4x10
# Subsequent lines define the element separated by tabs
def readmatrix(x):
    result=[]
    try:
        rows,cols=x[0].split("x")
        first=1
    except ValueError:
        first=0   # no header line
    for line in x[first:]:
        line=line.strip()
#        tmp=map(float,line.split('\t'))
        tmp=list(map(float,line.split()))
        result.append(tmp)
    return numpy.array(result)


def read_matrix_from_file(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    m = readmatrix(lines)
    return m

def readmatrixfile(filename= ""):
    if filename == "":
        lines=sys.stdin.readlines()
    else:
        with open(filename) as f:
            lines=f.readlines()

    m=readmatrix(lines)
    return m

def writematrixfile(x, filename):
    with open(filename, "w") as f:
        printmatrix(x, f)




def printmatrix(x, file=sys.stdout, headers=[], colheaders=[], format="%f", sep="\t"):
    rows, cols = x.shape
    printheaders =   len(headers) != 0
    printcolheaders =   len(colheaders) != 0
    assert(printheaders==False or len(headers) == rows) 
    if printcolheaders:
        if printheaders:
            file.write(sep)
        for j in range(cols-1):
            file.write("%s" % colheaders[j])
            file.write(sep)
        file.write("%s" % colheaders[cols-1])
        file.write("\n")

    for i in range(rows):
        if printheaders:
            file.write("%s%s" % (headers[i], sep))
        file.write(format % x[i,0])
        for j in range(1,cols):
            file.write(sep)
            file.write(format % x[i,j])
        file.write("\n")

def printintegermatrix(x, file, headers=[], colheaders=[]):
    printmatrix(x, file, headers, colheaders, format="%i")

# Normalize a PFM
def normalize(m):
    for i in range(0,m.shape[1]):
        if sum(m[:,i]) != 0:
            m[:,i] /= sum(m[:,i])
    return m

def is_integer_matrix(m):
    for c in range(0,m.shape[1]):
        for r in range(0,m.shape[0]):
            if not float.is_integer(m[r,c]):
                return False
    return True


# Entropy of a probability distribution 'l'
def entropy(l):
    sum=0
    for f in l:
        if f != 0:
            try:
                sum+=f*math.log(f,2)
            except ValueError:
                print(l)
                raise
    return -sum;

def information_content(l):
    return 2-entropy(l)

def matrix_information_content(m):
    columns=m.transpose().tolist()  # l is list of columns
    total_ic = 0.0
    for column in columns:
        total_ic += information_content(column) 
    return total_ic

# scale probabilities by information content
def logo_form(m):
    rows, cols = m.shape
    res=numpy.matrix(m)
    heights=[]
    # get information content for columns
    l=res.transpose().tolist()  # l is list of columns
    for c in range(cols):
        res[:,c] *= information_content(l[c]) 
    return res

# scale probabilities by 2. This is for background matrices
def logo_form2(m):
    rows, cols = m.shape
    res=numpy.matrix(m)
    for c in range(cols):
        for r in range(rows):
            res[r,c] *= 2
    return res


def complement(c):
    return iupac_complement[c]

def reverse_complement(x):
    x=x[::-1]
    result=list(map(complement, x))
    return "".join(result)

def jaccard_index(A, B):
    """Computes the Jaccard index for sets A and B"""
    return len(A & B) / len(A | B)


def hamming_distance(x, y):
    assert len(x) == len(y)
    s = 0
    for i in range(len(x)):
        c = (1 - jaccard_index(set(iupac_codes[x[i]]), set(iupac_codes[y[i]])))
        #print(x[i], y[i], c)
        s += c
    if float.is_integer(s):
        return int(s)
    else:
        return s

def is_palindrome(x):
    return x == reverse_complement(x)



def reverse_complement_pwm(m):
    result=m.copy()
    for c in range(0,m.shape[1]):
        for r in range(0,m.shape[0]):
            result[m.shape[0]-r-1,m.shape[1]-c-1] = m[r,c]
    return result


def mycommand(s):
    p = subprocess.Popen(s, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (output, error) = p.communicate()
    return (p.returncode, output, error)

def KL_distance(p, q):
    assert len(p) == len(q)
    sum=0
    for i in range(len(p)):
        sum += p[i]*math.log(p[i]/q[i], 2)
    return sum

def symmetric_KL_distance(p,q):
    return 0.5*KL_distance(p, q) + 0.5*KL_distance(q, p)

def matrix_symmetric_KL_distance(m1, m2):
    pseudo_count=0.00001
    x1=normalize(m1+pseudo_count)
    x2=normalize(m2+pseudo_count)
    return symmetric_KL_distance(x1.flatten(), x2.flatten())

def matrix_KL_distance(m1, m2):
    pseudo_count=0.00001
    x1=normalize(m1+pseudo_count)
    x2=normalize(m2+pseudo_count)
    return KL_distance(x1.flatten(), x2.flatten())


def mybinom(k, h):
    from scipy.special import binom
    if k==0 and h > 0:
        return 0   # scipy.special.binom returns -0 when k=0 and d > 0 is even
    else:
        return binom(k, h)
    
# The neighbourhood includes the strings with Hamming distance exactly the given radius or below
# k is the string length
def hamming_neighbourhood_size(k, radius):
    return sum([ mybinom(k,h) * 3**h for h in range(radius+1)])
    
# The neighbourhood includes the strings with Hamming distance exactly the given radius.
# k is the string length
def hamming_border_size(k, radius):
    h = radius
    return mybinom(k,h) * 3**h

    
