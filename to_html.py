#!/usr/bin/env python3

import adm
import getopt
import sys
import re
import subprocess
import os
import numpy as np
import string
import math
import io
import matplotlib.pyplot as plt
import matplotlib
import heatmap
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

css="""

h1 {
    margin-left: 100px;
}


ul {
    list-style-type: square;

}

#info {
    overflow: auto;
    background-color: lightblue;
    float: left; 
    border: solid green;
}

#programinfo {
    overflow: auto;
    background-color: lightblue;
    float: left; 
    clear: left;
    border: solid green;
    margin-top: 20px;
    margin-bottom: 20px;
    padding-right: 20px;
}

#programinfo ul {
    max-width: 60em;
}

div#lambdaTable th{
    text-align: left;
}


div#factors {
#    clear: left;
    float: left;
    text-align: center;
}
div#factors img{
#    clear: left;
#    float: left;
    padding: 20px;
}

div#flankfactors {
    clear: both;
}


#cobs {
    clear:left;
    float: left;
    text-align: center;
//    color: red;
}

div#iterations {
    float: left;
    #text-align: center;
    color: red;
}

div#iterations img{
    width: 200;
    height: 200;
}

div#background {
    overflow: auto;
    clear: left;
}

table.ppmtable th{
    background-color: lightblue;
}

table.ppmtable td,
table.ppmtable th
{
    padding: 10px;
    border: 1px solid black;
    text-align: center;
}


table.ppmtable  {
    border-collapse: collapse;
    border-spacing: 0px;
}

h3.tableheading {
    margin-top: 40px;
    margin-bottom: 10px;
}

div.rc {
    border: 2px solid white;
}
div.normal {
    border: 2px solid red;
}
div.logoContainer {vertical-align: middle;}
div.logoButtonContainer {display: inline-block; vertical-align: middle;}
div.normal { text-align: center; cursor: pointer;}
div.rc { text-align: center; cursor: pointer;}
div.logoImageContainer {display: inline-block; vertical-align: middle;}
"""


dna_orients=["HT", "HH", "TT", "TH"]
dna_orient_dict = {"HT" : 0, "HH" : 1, "TT" : 2, "TH" : 3}

rna_orients=["HT", "TH"]
rna_orient_dict = {"HT" : 0, "TH" : 1}

max_logo_width = 500 # This is defined in myspacek40

dont_create_visualizations=True

# Entropy of a probability distribution 'l'
def entropy(l):
    assert abs(sum(l) - 1.0) < 0.001, "The distribution must sum to 1.0, got %e" % sum(l)
    assert 0.0 <= min(l) and max(l) <= 1.0, "The values in the distribution must be between 0.0 and 1.0"
    s=0
    for f in l:
        if f != 0:
            try:
                s += f*math.log(f,2)
            except ValueError:
                print(l)
                raise
    return -s;

def information_content(l):
    return 2-entropy(l)

def adm_information_content(adm):
    cols = adm.k
    result=[]
    result.append(information_content(adm.initial_probabilities[:,0]))
    for c in range(cols-1):
        temp = 0.0
        for a in range(4):
            init = adm.initial_probabilities[a,c]
            if init > 0.0:
                try:
                    temp += init * information_content(adm.transition_probabilities[4*a:4*(a+1),c])
                except AssertionError as err:
                    print("Initial probability is %f, %s" % (init, err), file=sys.stderr)
                    #raise
                result.append(temp)
    return result

def matrix_information_content(m):
    columns=m.transpose().tolist()  # l is list of columns
    total_ic = 0.0
    for column in columns:
        total_ic += information_content(column) 
    return total_ic

def model_information_content(m):
    if is_adm(m):
        ics = adm_information_content(m)
        return sum(ics)
    else:
        return matrix_information_content(m)
    
# Normalize a PFM
def normalize(m):
    for i in range(0,m.shape[1]):
        if sum(m[:,i]) != 0:
            m[:,i] /= sum(m[:,i])
    return m

def is_adm(m):
    return type(m) == adm.adm

transform = [15, 11, 7, 3, 14, 10, 6, 2,
             13, 9, 5, 1, 12, 8, 4, 0]


def reverse_complement_adm(m):
    k = m.k
    i=reverse_complement_pwm(m.initial_probabilities)
    t = np.zeros((16, k-1))
    for j in range(k-1):
        for ab in range(16):
            a = ab // 4
            b = ab % 4
            divisor = m.initial_probabilities[b, j+1]
            if divisor > 0.0:
                t[transform[ab], k-j-2] = m.transition_probabilities[ab, j] * m.initial_probabilities[a, j] / divisor
    return adm.adm(t,i)

def reverse_complement_pwm(m):
    result=m.copy()
    for c in range(0,m.shape[1]):
        for r in range(0,m.shape[0]):
            result[m.shape[0]-r-1,m.shape[1]-c-1] = m[r,c]
    return result

def reverse_complement_model(m):
    if is_adm(m):
        return reverse_complement_adm(m)
    else:
        return reverse_complement_pwm(m)
    
def matrices_in_orientation(o, p1, p2):
    if o == "HT":
        return (p1, p2)
    elif o == "HH":
        return (p1, reverse_complement_model(p2))
    elif o == "TT":
        return (reverse_complement_model(p1), p2)
    elif o == "TH":
        return (reverse_complement_model(p1), reverse_complement_model(p2))

def find_lines(x, str, pos, count):
    resultfile= io.StringIO("".join(x))
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

def readarray(lines):
    result=[]
    for line in lines:
        line = line.rstrip('\n')
        tmp=line.split("\t")
        result.append(tmp)
    return np.array(result)

def readmodel(x):
    if len(x)==20:
        return adm.read_adm_from_list_of_lines(x)
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
    return np.array(result)

def right_extend_adm(m, extension):
    assert extension >= 0
    orig_k = m.shape[1]
    k = orig_k + extension
    result = np.zeros((16,k))
    result[:, 0:orig_k] = m.representation()
    result[:, orig_k:] = 0.25
    return adm.adm(result)

def left_extend_adm(m, extension):
    assert extension >= 0
    orig_k = m.shape[1]
    k = orig_k + extension
    result = np.zeros((16,k))
    result[:, extension:] = m.representation()
    result[4:8, extension] = result[8:12, extension] = result[12:16, extension] = result[0:4, extension]
    result[:, :extension] = 0.25
    return adm.adm(result)

def force_adms_equal(adm1, adm2):
    k = adm1.shape[1]
    assert k == adm2.shape[1]

    product = adm1.representation() * adm2.representation()

    r = np.zeros((4, k+1))
    for a in range(4):
        r[a, k] = 1.0
    
    for j in range(k-1, -1, -1):
        amax = 1 if j==0 else 4
        for a in range(amax):
            for b in range(4):
                r[a, j] += product[a*4+b, j] * r[b, j+1]

    result = np.zeros((16, k))
    for j in range(k):
        amax = 1 if j==0 else 4
        for a in range(amax):
            for b in range(4):
                if r[a, j] > 0.0:
                    result[4*a+b, j] = product[a*4+b, j] * r[b, j+1] / r[a, j]
                else:
                    assert product[a*4+b, j] * r[b, j+1] == 0.0
    return adm.adm(result)
        
def compute_expected_adm(adm1, adm2, o, d):
    adm1,adm2 = matrices_in_orientation(o, adm1, adm2)
    k1 = adm1.shape[1]
    k2 = adm2.shape[1]
    dimer_len = k1 + k2 + d
    a1 = right_extend_adm(adm1, dimer_len - k1)
    a2 = left_extend_adm(adm2, dimer_len - k2)
    return force_adms_equal(a1, a2)
    
def compute_expected_pwm(pwm1, pwm2, o, d):
    m1,m2 = matrices_in_orientation(o, pwm1, pwm2)
    k1 = m1.shape[1]
    k2 = m2.shape[1]
    dimer_len = k1+k2+d
    fill=1.00

    # Left occurrence
    result1 = np.empty((4, dimer_len))
    result1.fill(fill)
    for pos in range(0, k1):
        result1[:,pos] = m1[:,pos]

    # Right occurrence
    result2 = np.empty((4, dimer_len))
    result2.fill(fill)
    for pos in range(k1+d, dimer_len):
        result2[:,pos] = m2[:,pos-(k1+d)]

    expected = normalize(result1 * result2)
    return expected

def compute_expected(pwm1, pwm2, o, d):
    if is_adm(pwm1):
        return compute_expected_adm(pwm1, pwm2, o, d)
    else:
        return compute_expected_pwm(pwm1, pwm2, o, d)
    
# Write results for a cob case
def write_results(cob, o, d, pwm1, pwm2, observed, expected, deviation, last_iteration_output, get_flanks, motif_ending):
    k1 = pwm1.shape[1]
    k2 = pwm2.shape[1]
    if use_adm:
        rows=20
        start=1
    else:
        rows=4
        start=1

    if get_flanks:
        try:
            flank=readmodel(find_lines(last_iteration_output, "Flank dimer case matrix %s %s %i:" % (cob, o, d), start, rows))
        except AttributeError:
            flank=np.zeros(expected.shape)

    oname="observed.%s.%s.%i.%s" % (cob, o, d, motif_ending)
    ename="expected.%s.%s.%i.%s" % (cob, o, d, motif_ending)
    dname="deviation.%s.%s.%i.dev" % (cob, o, d)
    fname="flank.%s.%s.%i.%s" % (cob, o, d, motif_ending)

    writematrixfile(observed, oname)
    writematrixfile(expected, ename)
    writematrixfile(deviation, dname)
    if get_flanks:
        writematrixfile(flank, fname)



    oname_rc="observed.%s.%s.%i-rc.%s" % (cob, o, d, motif_ending)
    ename_rc="expected.%s.%s.%i-rc.%s" % (cob, o, d, motif_ending)
    dname_rc="deviation.%s.%s.%i-rc.dev" % (cob, o, d)
    fname_rc="flank.%s.%s.%i-rc.%s" % (cob, o, d, motif_ending)

    observed_rc=reverse_complement_model(observed)
    expected_rc=reverse_complement_model(expected)
    deviation_rc=reverse_complement_model(deviation)
    if get_flanks:
        flank_rc=reverse_complement_model(flank)

    writematrixfile(observed_rc, oname_rc)
    writematrixfile(expected_rc, ename_rc)
    writematrixfile(deviation_rc, dname_rc)
    if get_flanks:
        writematrixfile(flank_rc, fname_rc)



    if not dont_create_visualizations:
        # Forward direction
        myrun("myspacek40 %s --logo %s %s" % (myspacek_flags, oname, oname.replace(".%s"%motif_ending, ".svg")))
        myrun("myspacek40 %s --logo %s %s" % (myspacek_flags, ename, ename.replace(".%s"%motif_ending, ".svg")))
        #if not use_adm:         # This does not work for adm models
        myrun("myspacek40 %s --difflogo %s %s %s" % (myspacek_flags, oname, ename, dname.replace(".dev", ".svg")))          # Deviation logo
        if get_flanks:
            myrun("myspacek40 %s -core=%i,%i,%i --logo %s %s" % (myspacek_flags, k1, k2, d, fname, fname.replace(".%s"%motif_ending, ".svg")))

        # Reverse complement
        myrun("myspacek40 %s --logo %s %s" % (myspacek_flags, oname_rc, oname_rc.replace(".%s"%motif_ending, ".svg")))
        myrun("myspacek40 %s --logo %s %s" % (myspacek_flags, ename_rc, ename_rc.replace(".%s"%motif_ending, ".svg")))
    #    if not use_adm:         # This does not work for adm models
        myrun("myspacek40 %s --difflogo %s %s %s" % (myspacek_flags, oname_rc, ename_rc, dname_rc.replace(".dev", ".svg")))      # Deviation logo
        if get_flanks:
            myrun("myspacek40 %s -core=%i,%i,%i --logo %s %s" % (myspacek_flags, k2, k1, d, fname_rc, fname_rc.replace(".%s"%motif_ending, ".svg")))

    for rc in ["", "-rc"]:
        with open("three.%s.%s.%i%s.html" % (cob, o, d, rc), "w") as f:
            oname = "observed.%s.%s.%i%s" % (cob, o, d, rc)
            ename = "expected.%s.%s.%i%s" % (cob, o, d, rc)
            dname = "deviation.%s.%s.%i%s" % (cob, o, d, rc)
            #if use_adm:
            #    myrun("mv %s.adm_minus_%s.adm.svg %s.svg" % (oname, ename, dname))
            #else:
            #    myrun("mv %s.pfm_minus_%s.pfm.svg %s.svg" % (oname, ename, dname))
            f.write('<h1>%s %s %i</h1>' % (cob, o, d))
            f.write('<figure><figcaption>Observed:</figcaption><a href="%s.%s"><img src="%s.svg"\></a></figure>' % (oname, motif_ending, oname)) 
            f.write('<figure><figcaption>Expected:</figcaption><a href="%s.%s"><img src="%s.svg"\></a></figure>' % (ename,motif_ending,ename)) 
            f.write('<figure><figcaption>Deviation:</figcaption><a href="%s.dev"><img src="%s.svg"\></a></figure>' % (dname,dname)) 




def get_cob_case(cob, o, d, pwm1, pwm2, last_iteration_output, get_flanks, motif_ending):        
    expected = compute_expected(pwm1, pwm2, o, d)
    if use_adm:
        rows=16
        start=1
    else:
        rows=4
        start=1
    try:
        deviation=readmodel(find_lines(last_iteration_output, "Deviation matrix %s %s %i:" % (cob, o, d), start, rows))
    except AttributeError:
        deviation=np.zeros(expected.shape)
    if use_adm:
        observed = expected.representation() + deviation
    else:
        observed = expected + deviation
    g = np.vectorize(lambda x : max(x,0)) # This cuts negative values to 0
    observed = normalize(g(observed)) # Because precision of (about) 6 digits is used, some elements can be slightly negative
    if use_adm:
        observed = adm.adm(observed)
        
    write_results(cob, o, d, pwm1, pwm2, observed, expected, deviation, last_iteration_output, get_flanks, motif_ending)
    return observed


def seconds_to_hms(seconds, use_si_units = True):
    fraction=seconds-int(seconds)
    seconds=int(seconds)
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    if use_si_units:
        if h != 0:
            result = "%dh %02dm %02.2fs" % (h, m, s+fraction)
        elif m != 0:
            result = "%dm %02.2fs" % (m, s+fraction)
        else:
            result = "%.2fs" % (s+fraction)
    else:
        if h != 0:
            result = "%d:%02d:%02.2f" % (h, m, s+fraction)
        elif m != 0:
            result = "%d:%02.2f" % (m, s+fraction)
        else:
            result = "%.2f" % (s+fraction)
    return result

def extract(query, list):
    m = re.search(query, list)
    try:
        value = m.group(1)
    except:
        value = None
    return value

def extract_list(query, list):
    return re.findall(query, list)

def logo_container(anchor, image, ending, title=""):
    if len(title) > 0:
        attr='title="%s"' % title
    else:
        attr=''
    return """<div class="logoContainer">
                   <div class="logoButtonContainer">
	             <div class="normal"  onclick='myclick(event, "%s")'>+</div>
	             <div class="rc" onclick='myclick(event, "%s")'>-</div>
                   </div>
                   <div class="logoImageContainer" >
	             <a href="%s"><img class="image" src="%s" %s></a>
                   </div>
                 </div>""" % (ending, ending, anchor, image, attr)

# f ends with .svg, and is the filename of the logo
#def print_logo_container(f, logo, ending, title=""):
#    f.write(logo_container(logo.replace(".svg", ending), logo, ending, title))

def mycommand(s):
    p = subprocess.Popen(s, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
    (output, error) = p.communicate()
    return (p.returncode, output, error)

def make_table_h(f, files, headers, motif_ending, titles=[]):
    f.write("<table>")
    if len(headers) > 0:
        f.write("<tr>")
        for h in headers:
            f.write("<th>")
            f.write(h)
            f.write("</th>")
        f.write("</tr>")
    f.write("<tr>")
    if len(titles) != len(files):
        titles = [""]*len(files)
    for logo, title in zip(files, titles):
        link=logo.replace(".svg", ".%s"%motif_ending)
        f.write("<td>")
        f.write(logo_container(link, logo, ".%s"%motif_ending, title))
        f.write("</td>")
    f.write("</tr>")
    f.write("</table>")

def make_table_v(files, f, header=""):
    f.write("<table>")
    # if len(header) > 0:
    #     print "<tr>"
    #     print "<th>"
    #     print f
    #     print "</th>"
    #     print "</tr>"
    for plot in files:
        f.write("<tr>")
        f.write("<td>")
        f.write('<a href="%s"><img src="%s"\></a>' % (plot,plot))
        f.write("</td>")
        f.write("</tr>")
    f.write("</table>")

def make_table_v2(f, files, headers=[], links=[], titles=[]):
    if len(headers) != len(files):
        headers = [""]*len(files)
    if len(links) != len(files):
        links = [""]*len(files)
    if len(titles) != len(files):
        titles = [""]*len(files)
    f.write("<table>")
    for x,h,l in zip(files, headers, links):
        f.write("<tr>")
        if len(h) > 0:
            f.write("<th>%s</th>" % h)
        f.write("<td>")
        if l == "":
            l=x
        f.write('<a href="%s"><img src="%s"\></a>' % (l,x))
        f.write("</td>")
        f.write("</tr>")
    f.write("</table>")

def make_table_v3(f, files, headers, links, ending, titles=[]):
    if len(titles) != len(files):
        titles = [""]*len(files)
    f.write("<table>")
    for x,h,link, title in zip(files, headers, links, titles):
        f.write("<tr>")
        f.write("<th>%s</th>" % h)
        f.write("<td>")
        f.write(logo_container(link, x, ending, title))
#f.write('<a href="%s"><img src="%s"\></a>' % (l,f))
        f.write("</td>")
        f.write("</tr>")
    f.write("</table>")

def make_table(table, headers, f):
    f.write("<table>")
    if len(headers) > 0:
        f.write("<tr>")
        for h in headers:
            f.write("<th>")
            f.write(h)
            f.write("</th>")
        f.write("</tr>")
    for row in table:
        f.write("<tr>")
        for c in row:
            f.write("<td>")
            f.write(str(c))
            f.write("</td>")
        f.write("</tr>")
    f.write("</table>")

def make_better_table(f, table, headers=[], row_headers=[], htmlclass=""):
    if len(htmlclass) > 0:
        f.write('<table class="%s">' % htmlclass)
    else:
        f.write("<table>")
    if len(headers) > 0:
        f.write("<tr>")
        for x in headers:
            f.write("<th>%s</th>" % x)
        f.write("</tr>")
    for i, row in enumerate(table):
        f.write("<tr>")
        if len(row_headers) > 0:
            f.write("<th>%s</th>" % row_headers[i])
        for c in row:
            f.write("<td>")
            f.write(c)
            f.write("</td>")
        f.write("</tr>")
    f.write("</table>")
    

def get_best_cob_cases(cob_codes):
    best_cases=[]                                                          # ["best.observed.0-0.svg", ...]
    best_cases_links=[]                                                    # ["best.0-0.html", ...]
    best_cases_headers=[]                                                  # ["best.observed.0-0.svg", ...]
    for c in cob_codes:
        best_cases.append("best.observed.%s.svg" % c)
        best_cases_links.append("best.%s.html" % c)
        try:
            best_cases_headers.append(os.readlink("best.observed.%s.svg" % c))
        except OSError:
            best_cases_headers.append("best.observed.%s.svg(nonexisting)"%c)
    return best_cases, best_cases_links, best_cases_headers


def get_lambda_table(results_output, factors, cobs, cob_codes):
    lambda_headers = factors + cobs + ["Background", "Sum"]
    # lambda_table = [ [factors[i], 0.0] for i in xrange(number_of_factors) ]  + \
    #     [ [cobs[i],   0.0]  for i in xrange(number_of_cobs) ] +\
    #     [ ["Background", 0.0], ["Sum", 0.0] ]

    bg_lambda = float(extract(r"Background lambda is (.*)", results_output))
    temp = extract(r"Monomer lambdas are (.*)", results_output) 
    monomer_lambdas = eval(temp)
    # for i,dummy in enumerate(factors):
    #     lambda_table[i][1] = monomer_lambdas[i]

    number_of_cobs = len(cobs)
    cob_lambdas = [0] * number_of_cobs
    for i, c in enumerate(cob_codes):
        cob_lambdas[i] = float(extract(r"Sum of dimer lambdas of cob table %s is (.*)" % c, results_output))

    lambdas = monomer_lambdas + cob_lambdas + [bg_lambda]
    lambda_sum = sum(lambdas)
    lambdas.append(lambda_sum)
    lambda_table2 = list(zip(lambda_headers, lambdas))
    
    return lambda_table2

def get_info(results_output, full_output, cob_codes):
    maxiter = int(extract(r"Maximum number of iterations is (.*)", full_output))
    iterations =  int(extract(r"EM-algorithm took (.*) iterations", results_output))
    Lmin =  int(extract(r"Minimum sequence length is (.*)", full_output))
    Lmax =  int(extract(r"Maximum sequence length is (.*)", full_output))
    try:
        lines = int(extract(r"Using (.*) sequences", full_output))
    except TypeError:
        lines = int(extract(r"Read (.*) lines from file",full_output))

    command = extract(r"Command line was: (.*)", full_output)
    start_time = extract(r"Starting program at (.*)", full_output)
    version = extract(r"MODER version (.*)", full_output)
    hostname = extract(r"Running on host: (.*)", full_output)
    threads = extract(r"Using (.*) openmp threads", full_output)
    epsilon = float(extract(r"Epsilon is (.*)", full_output))
    excluded = 0
    for c in cob_codes:
        excluded += len(extract(r"Cob %s excluded cases: \[(.*)\]" % c, results_output).split(", "))
    return iterations, maxiter, Lmin, Lmax, lines, epsilon, excluded, command, start_time, version, hostname, threads


def get_seeds(full_output, number_of_factors):
    temp=extract_list(r"Monomer seeds are \[(.+)\]", full_output)
    temp2=[x.split(", ") for x in temp]
    seeds_begin = extract(r"Initial monomer seeds are \[(.+)\]", full_output).split(", ")
    #seeds_begin = ["GACCGGAAGCG", "CACCTG"]

    try:
        seeds_end=temp2[-1]
    except IndexError:
        seeds_end=seeds_begin

    with open("seeds.txt", "w") as f:
        f.write("Initial %s\n" % " ".join(seeds_begin))
        try:
            prev=temp2[0]
        except IndexError:
            pass
        for i, t in enumerate(temp2): 
            ch=[" "] * number_of_factors
            for j in range(number_of_factors): 
                ch[j] = str(j+1) if prev[j] != t[j] else " "
            f.write("Round %02i %s\t%s\n" % (i, " ".join(t), " ".join(ch))) # Third field contains a number for each factor that had its seed changes
                                                                            # compared to the one from the previous iteration

            prev = t
    return seeds_begin, seeds_end



def get_run_time(results_output):
    temp=extract(r"Whole program took (.+) seconds wall-time", results_output)
    try:
        s=seconds_to_hms(float(temp))
    except TypeError:
        s="unknown"
    return s


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

def writematrixfile(x, filename):
    with open(filename, "w") as f:
        if is_adm(x):
            f.write(x.str2())
        else:
            printmatrix(x, f)


def myrun(command):
    result=os.system("%s 2>&1 > /dev/null" % command)
    if result == 0:
        print("SUCCESS %s" % command)
    else:
        print("FAILURE %s" % command)


def myargmax(a):
    return np.unravel_index(np.argmax(a), a.shape)

def add_pseudo_counts(model, pseudo_count):
    if use_adm:
        a = model.representation()
        k = model.shape[1]
        for col in range(k):
            for row in range(4):
                idx = np.s_[row*4:(row+1)*4, col]  # this could be used instead of repeating the slicing notation
                s = np.sum(a[row*4:(row+1)*4, col])
                if s > 0.0:
                    a[row*4:(row+1)*4, col] += pseudo_count
                    s = s + 4*pseudo_count
                    a[row*4:(row+1)*4, col] /= s
        model = adm.adm(a)
    else:
        model += pseudo_count
        model /= model.sum(axis=0)
    return model

def get_monomers(factors, results_output, last_iteration_output, model_rows, motif_ending, get_flanks):
    factor_lengths = [0]*len(factors)
    factor_ics = [0]*len(factors)
    factor_models = [0]*len(factors)
#    start = 1 if use_adm else 2
    start = 1
    for i, factor in enumerate(factors):
        lines=find_lines(results_output, "Monomer matrix %i:" % i, start, model_rows)
#        with open("%s.pfm" % factor, "w") as f:
        with open("monomer.%i.%s" % (i,motif_ending), "w") as f:
            f.writelines(lines)
        model = readmodel(lines)
        model = add_pseudo_counts(model, 0.000001) # Because float were stored it the file using only 6 decimals, some entries maybe zero, this fixes it
        factor_models[i] = model
        model_rc=reverse_complement_model(model)
#        writematrixfile(pwm_rc, "%s-rc.pfm" % factor)
        writematrixfile(model_rc, "monomer.%i-rc.%s" % (i,motif_ending))
        factor_lengths[i] = model.shape[1]
        factor_ics[i] = model_information_content(model)
        if get_flanks:
            lines=find_lines(last_iteration_output, "Flank monomer matrix %i:" % i, start, model_rows)
            with open("flank-%i.%s" % (i,motif_ending), "w") as f:
                f.writelines(lines)
            flank_model=readmodel(lines)
            flank_model_rc=reverse_complement_model(flank_model)
            writematrixfile(flank_model_rc, "flank-%i-rc.%s" % (i,motif_ending))
    return factor_lengths, factor_ics, factor_models

def get_dimer_cases(results_output, iterations, last_iteration_output, cob_factors, use_rna, cob_codes, dmin, dmax, cob_tables, cob_ic_tables,
                    cob_length_tables, orients, orient_dict, monomer_pwms, get_flanks, motif_ending, cob_titles, best_cases_titles):
    for i, cob_factor in enumerate(cob_factors):
        number_of_orientations = 1 if use_rna else 3
        
        if cob_factor[0] != cob_factor[1]:
            number_of_orientations += 1
            
        tf1, tf2 = cob_factor
        # Find the cob table
        lines=find_lines(results_output, "Dimer lambdas %s:" % cob_codes[i], 1, number_of_orientations + 1)
        temp = readarray(lines)
        dmin[i] = int(temp[0,1])
        dmax[i] = int(temp[0,-1])
        cob_tables[i] = temp[1:,1:].astype(float)  # Remove column and row names
        #print cob_tables[i]
        cob_ic_tables[i]=np.zeros(cob_tables[i].shape)      # Information contents of corresponding dimeric PWM
        cob_length_tables[i]=np.zeros(cob_tables[i].shape)  # Lengths of corresponding dimeric PWM
        #print len(cob_ic_tables)
        #print cob_ic_tables[i].shape
        #overlap_table = cob_tables[i][:,0:-dmin[i]] # Only the overlap area of the cob table
        try:
            oi,di = myargmax(cob_tables[i])
            empty=False
            best_o = orients[oi]  # Maximum lambda cob case is (o,d)
            best_d = dmin[i]+di
        except ValueError:
            empty=True

        excluded=extract(r"Cob %s excluded cases: \[(.*)\]" % cob_codes[i], results_output)
        if excluded:
            excluded=excluded.split(", ")
        # for s in excluded:
        #     o2,d2 = s.split(" ")
        #     temp[1+orient_dict[o2],1+int(d2)-dmin[i]] = -0.0002   # Gets gray colour in heatmap

        # This is for the best case in this cob table
        if float(temp[1+orient_dict[best_o],1+int(best_d)-dmin[i]]) > 0.0:
            os.system("ln -f -s observed.%s.%s.%i.svg best.observed.%s.svg" % (cob_codes[i], best_o, best_d, cob_codes[i]))
            os.system("ln -f -s observed.%s.%s.%i-rc.svg best.observed.%s-rc.svg" % (cob_codes[i], best_o, best_d, cob_codes[i]))
    #        os.system("ln -f -s expected.%s.%s.%i.svg best.expected.%s.svg" % (cob_codes[i], best_o, best_d, cob_codes[i]))
    #        os.system("ln -f -s deviation.%s.%s.%i.svg best.deviation.%s.svg" % (cob_codes[i], best_o, best_d, cob_codes[i]))
            os.system("ln -f -s three.%s.%s.%s.html best.%s.html" % (cob_codes[i], best_o, best_d, cob_codes[i]))
            os.system("ln -f -s three.%s.%s.%s-rc.html best.%s-rc.html" % (cob_codes[i], best_o, best_d, cob_codes[i]))

        # Compute IC and length of each dimeric PPM model related to cob_tables[i].
        # They are used in the hovering tool tip in the resulting html page.
        for row in range(0, number_of_orientations):
            for d in range(dmin[i], dmax[i]+1):
                if float(temp[1+row,1+d-dmin[i]]) > 0.00000:
                    #command="get_cob_case.py %s %i %s %s %i %s" % ("-f" if get_flanks else "", iterations-1, cob_codes[i], orients[row], d, inputfile)
                    #myrun(command)
                    try:
                        # The following return the observed model, but also computes and stores the expectations and deviations
                        dimer_pwm=get_cob_case(cob_codes[i], orients[row], d, monomer_pwms[tf1], monomer_pwms[tf2],
                                               results_output, get_flanks, motif_ending)
                        ic = model_information_content(dimer_pwm)
                    except AssertionError as error:
                        print("Error with o=%s d=%i: %s!" % (orients[row], d, error), file=sys.stderr)
                        raise
                    cob_ic_tables[i][row, d-dmin[i]] = ic
                    cob_length_tables[i][row, d-dmin[i]] = dimer_pwm.shape[1]

        # Write the cob table to a file
        with open("cob.%s.cob" % cob_codes[i], "w") as f:
            for row in range(temp.shape[0]):
                f.write("\t".join(temp[row,:]))
                f.write("\n")
        g = np.vectorize(lambda x,y,z : "Lambda %f&#010;IC: %f&#010;Length: %i" % (x,y,z)) # The hex number 0x10 is the new line character
        cob_titles[i] = g(cob_tables[i], cob_ic_tables[i], cob_length_tables[i])
        best_cases_titles[i]=cob_titles[i][oi, di]

def create_monomer_logos(factors, factor_lengths, motif_ending, get_flanks):
    # Monomers
    for i, f in enumerate(factors):
        #os.system("to_logo.sh -n -t %s %s.pfm" % (f, f))
#        myrun("myspacek40 -noname -paths --logo %s.pfm %s.svg" % (f, f))
#        myrun("myspacek40 -noname -paths --logo %s-rc.pfm %s-rc.svg" % (f, f))
        if not dont_create_visualizations:
            myrun("myspacek40 %s --logo monomer.%i.%s monomer.%i.svg" % (myspacek_flags, i, motif_ending, i))
            myrun("myspacek40 %s --logo monomer.%i-rc.%s monomer.%i-rc.svg" % (myspacek_flags, i, motif_ending, i))
            if get_flanks:
                g="flank-%i" % i
                myrun("myspacek40 %s -core=%i --logo %s.%s %s.svg" % (myspacek_flags, factor_lengths[i], g, motif_ending, g))
                myrun("myspacek40 %s -core=%i --logo %s-rc.%s %s-rc.svg" % (myspacek_flags, factor_lengths[i], g, motif_ending, g))




            
def visualize_cobs(cobs, cob_codes, cob_tables, dmin, dmax, orients):            
    for i, (cob, code) in enumerate(zip(cobs, cob_codes)):
        f = "cob.%s" % code
        #    myrun('heatmap.R -z 8 -c -s -f "%%.3f" -t %s %s.cob' % (f, f))
        data=cob_tables[i]
        vfunc = np.vectorize(lambda x: x if x > 0.0 else -0.0002)
        data=vfunc(data)
        drange = list(range(dmin[i], dmax[i]+1))
        if not dont_create_visualizations:
            print("Creating heatmap for %s" % cob)
            heatmap.make_heatmap(data, drange, orients, "svg", cob, "%s.svg" % f, fontsize=20.0, cell_labels=True)
#        myrun('heatmap.R -z 12 -c -s -i -t %s %s.cob 2> /dev/null > /dev/null' % (cob, f))
#        myrun("sed -i '/page/d' %s.svg" % f)  # R or its pheatmap package make svg files corrupt. This fixes it.
    #    myrun('heatmap.R -z 8 -c -s -f "%%.3f" -t %s %s.cob' % (f, f))


###################################################################################
#
#  Create distance, information content, lambda, log likelihood and parameter plots
#
###################################################################################

def myplot(data, title="", xlab="", ylab="", ymax=None, headers=[], outputfile=""):
    linewidth=2.0
    fontsize=32.0
    labelfontsize=fontsize*0.6
    tickfontsize=fontsize*0.6
    fig = plt.figure()
    ax = plt.subplot(111)
    if title:
        plt.title(title, fontsize=fontsize, fontname='sans-serif')
#    if ylab=="mll":
#        plt.ylabel(ylab, fontsize=labelfontsize, labelpad=20)
#    else:
    plt.ylabel(ylab, fontsize=labelfontsize)
    plt.xlabel(xlab, fontsize=labelfontsize)
    plt.xticks(fontsize=tickfontsize)
    plt.yticks(fontsize=tickfontsize)
    if not ymax is None:
        a = ymax*0.05            # Add margin to y-axis
        plt.ylim(0.0-a, ymax+a)
    ax.plot(data, linewidth=linewidth)
    ax.margins(x=0.05)           # Add margin to x-axis
    plt.grid(True)
    if headers:
        # Shrink current axis by 30%
        box = ax.get_position()
        if ylab=="mll":           # Make some more room for ylabel because of longer ytick labels
            ax.set_position([box.x0+0.10*box.width, box.y0, box.width * 0.7, box.height])
        else:
            ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
        # Put a legend to the right of the current axis
        plt.legend(headers, title="Models", loc='center left', bbox_to_anchor=(1, 0.5))
    if outputfile:
        plt.savefig(outputfile, format="svg")
    else:
        plt.show()

def create_graphs(full_output, factors, cobs, cob_codes, name):
    temp=extract_list(r"Log likelihood is (.+)", full_output)
    mll_header=["Log likelihood"]
    mll_data=[]
    with open("mll.txt", "w") as f:
        f.write("%s\n" % ("\t".join(mll_header)))
        for t in temp:
            mll_data.append(float(t))
            f.write("%s\n" % t)

    temp=extract_list(r"Total number of parameters is (.+)", full_output)
    parameters_header=["Number of parameters"]
    parameters_data=[]
    with open("parameters.txt", "w") as f:
        f.write("%s\n" % ("\t".join(parameters_header)))
        for t in temp:
            parameters_data.append(float(t))
            f.write("%s\n" % t)

    temp=extract_list(r"Intermediate average information content of monomer models: \[(.+)\]", full_output)
    ics_header=factors
    ics_data=[]
    with open("ics.txt", "w") as f:
        f.write("%s\n" % ("\t".join(ics_header)))
        for t2 in temp:
            t = t2.split(", ")
            ics_data.append(list(map(float, t)))
            f.write("%s\n" % "\t".join(t))

    distances_header=factors+cobs
    distances_data=[]
    with open("distances.txt", "w") as f:
        cols=len(factors+cobs)
        temp=extract_list(r"Monomer distances are \[(.+)\]", full_output)
        rows=len(temp)
        a=np.empty((rows,cols))
        for r,t in enumerate(temp):
            for c,x in enumerate(t.split(", ")):
                a[r,c]=float(x)
        for c, cob in enumerate(cob_codes):
            query=r"Max distance in deviation %s is (.+)" % cob
            #print query
            temp=extract_list(query, full_output)
            #print temp
            for r,t in enumerate(temp):
                a[r,len(factors)+c] = float(t)
        f.write("%s\n" % ("\t".join(distances_header)))
        #printmatrix(a, sys.stdout, headers=[], colheaders=[], format="%f", sep="\t")
        printmatrix(a, f, headers=[], colheaders=[], format="%f", sep="\t")
        distances_data = a.tolist()

    lambdas_header=factors+cobs+["bg"]
    lambdas_data=[]
    with open("lambdas.txt", "w") as f:
        temp=extract_list(r"Intermediate monomer lambdas are \[(.+)\]", full_output)
        atemp=[t.split(", ") for t in temp]
        #print atemp
        a=np.array(atemp).astype(float)
        btemp=[]
        for cob in cob_codes:
            #print "Cob is %s" % cob
            temp=extract_list(r"Intermediate sum of dimer lambdas %s is (.+)" % cob, full_output)
            btemp.append(temp)
    #    print btemp
        b=np.transpose(np.array(btemp)).astype(float)
        
        temp=extract_list(r"Intermediate background lambda is (.+)", full_output)
        c=np.transpose(np.array([temp])).astype(float)
        f.write("%s\n" % ("\t".join(lambdas_header)))

        print(a.shape, b.shape, c.shape)
        if len(cob_codes) > 0:
            d=np.concatenate((a,b,c), 1)
        else:
            d=np.concatenate((a,c), 1)
        #printmatrix(d, format="%s")
        printmatrix(d, f, headers=[], colheaders=[], format="%s", sep="\t")
        lambdas_data=d.tolist()

    myplot(mll_data,        "%s log likelihood" % name,       "iterations", "mll",                headers=lambdas_header,    outputfile="mll.svg")
    myplot(parameters_data, "%s number of parameters" % name, "iterations", "params",             headers=parameters_header, outputfile="parameters.svg")
    myplot(ics_data,        "%s information content" % name,  "iterations", "IC", ymax=2.0,       headers=ics_header,        outputfile="ics.svg")
    myplot(distances_data,  "%s convergence" % name,          "iterations", "distance", ymax=1.0, headers=distances_header,  outputfile="distances.svg")
    myplot(lambdas_data,    "%s lambdas" % name,              "iterations", "lambda", ymax=1.0,   headers=lambdas_header,    outputfile="lambdas.svg")


def get_start_positions(last_iteration_output, factors):
    lines = last_iteration_output.split('\n')
    for i, line in enumerate(lines):
        if "Monomer 0 start position distribution (forward):" in line:
            break
    if i < len(lines)-1:
        global start_positions
        start_positions = True
        temp = lines[i+1:i+len(factors)*4:2]
        temp = [x.strip('[]') for x in temp]
        temp = [x.replace(", ", "\t") for x in temp]
        with open("start_positions.tsv", "w") as f:
            for line in temp:
                f.write("%s\n" % line)
        temp = [x.split() for x in temp]
        data = np.array(temp).astype(float)
        ylabels=[]
        for i, f in enumerate(factors):
            ylabels.append("Monomer %i (forward)" % i)
            ylabels.append("Monomer %i (backward)" % i)
        xlabels = list(range(0, data.shape[1]))
        heatmap.make_heatmap(data, xlabels, ylabels, "svg", "Monomer start positions", "start_positions.svg", fontsize=20.0)

def get_monomer_modularity(full_output):
    temp=extract_list(r"Is monomer pwm learnt purely modularly: \[(.+)\]", full_output)
    temp2=[x.split(", ") for x in temp]
    return temp2[-1]

def main():        
    # This is to locate helper scripts in the same directory as this file
    execdir=os.path.abspath(os.path.dirname(sys.argv[0]))
    #print "execdir is %s" % execdir
    path=os.getenv("PATH")
    os.putenv("PATH", path+":"+execdir)

    usage="""Usage:
    \tto_html.py tf1name,tf2name,... moderoutputfile [ moderreportdir ]

    Parses the output created by MODER and converts it to a graphical
    html page. The first parameter is a comma separated list
    of factor names. The second parameter is the name of the
    MODER output file.

    A directory 'moderoutputfile.report' is created, and
    the report will be in file 'moderoutputfile.report/index.html'
    """

    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'hd', ["help", "debug"])
    except getopt.GetoptError as e:
        print(e)
        sys.stderr.write(usage)
        sys.exit(1)

    optdict = dict(optlist)
    args = [sys.argv[0]]+ args
    debug=False
    start_positions = False

    #print optdict
    for o, a in optlist:
            if o in ("-h", "--help"):
                print(usage)
                sys.exit(0)
            elif o in ("-d", "--debug"):
                debug=True
                print("Debugging on")
            else:
                sys.stderr.write("Unknown option: %s\n" % o)
                sys.stderr.write(usage)
                sys.exit(1)



    if len(args) == 1:
        print(usage)
        sys.exit(0)
    elif len(args) < 3:
        sys.stderr.write("Error, give at least two parameters.\n")
        sys.stderr.write(usage)
        sys.exit(1)

    name=args[1]
    orig=inputfile=args[2]

    if len(args) == 4:   # name of report directory given as third parameter
        mydir=args[3]
    else:
        if inputfile[1:].count(".") > 0:                         # Contains a file extension
            mydir=re.sub("\.[^.]*?$", ".report", inputfile)      # Strip the extension, and use .report instead
        else:
            mydir="%s.report" % inputfile
    reportfile="%s/index.html" % mydir

    inputfile = "../%s" % os.path.basename(inputfile)

    try:
        os.mkdir(mydir)
    except OSError:
        pass

    os.chdir(mydir)

    factors=name.split(',')  # For example ['FLI1a', 'FLI1b']

    try:
        with open(inputfile) as f:
            full_output = "".join(f.readlines())
    except IOError:
        sys.stderr.write("Could not read file %s. Exiting!\n" % orig)
        sys.exit(1)

    binding_model = extract(r"Using binding model: (.*)", full_output)
    global use_adm
    if binding_model == "ppm":
        use_adm=False
        motif_ending="pfm"
        model_rows=4
    else:
        use_adm=True
        motif_ending="adm"
        model_rows=20

    use_rna = extract(r"Use RNA alphabet: (.*)", full_output)
    global myspacek_flags
    if use_rna == "yes":
        use_rna = True
        orients = rna_orients
        orient_dict = rna_orient_dict
        myspacek_flags="-paths -noname -rna"
    else:
        use_rna = False
        orients = dna_orients
        orient_dict = dna_orient_dict
        myspacek_flags="-paths -noname"

    cob_factors=extract(r"Cob combinations are ([0-9,-]*)", full_output)
    if cob_factors:
        cob_factors=cob_factors.split(",")
    else:
        cob_factors=[]
    cob_factors=[list(map(int, x.split("-"))) for x in cob_factors]              # cob_factors=[[0,0], [1,1], [0,1]]
    number_of_cobs=len(cob_factors)
    factor_codes=set()
    for x,y in cob_factors:
        factor_codes.add(x)
        factor_codes.add(y)
    # In the if clause below, if only one name e.g. HNF4A is given and four codes 0,1,2,3 are used in cob types, then form names
    # HNF4Aa,HNF4Ab,HNF4Ac,HNF4Ad
    if len(factors) == 1 and len(factor_codes) > 1:
        new_factors=[factors[0]+chr(ord('a')+x) for x in range(max(factor_codes)+1) ]
        factors = new_factors

    number_of_factors=len(factors)



    cobs=['-'.join([factors[x[0]], factors[x[1]]]) for x in cob_factors ]  # cobs=["TEAD4-TEAD4", "ERG-ERG",  "TEAD4-ERG"]
    cob_codes=[ "-".join(map(str,x)) for x in cob_factors ]                # cob_codes=["0-0", "1-1", "0-1"]


    command="sed -n '/Results/,$p' %s" % inputfile                         # KORJAA TÄMÄ!!!!!!!!!! MUUALLAKIN KÄYTETÄÄN SEDIÄ, POISTA NE
    (ret_val, results_output, error) = mycommand(command)
    assert ret_val == 0


    ######################################################################################################
    #
    # to_logos
    #


    cob_tables=[0]*number_of_cobs
    cob_ic_tables=[0]*number_of_cobs
    cob_length_tables=[0]*number_of_cobs
    cob_titles=[0]*number_of_cobs
    dmin=[0]*number_of_cobs
    dmax=[0]*number_of_cobs

    get_flanks = extract(r"Get full flanks: (.*)", full_output) == "yes"

    iterations =  int(extract(r"EM-algorithm took (.*) iterations", results_output))

    command="sed -n '/Round %i/,/^-+$/p' %s" % (iterations-1,inputfile)   # Output from last iteration onwards
    (ret_val, last_iteration_output, error) = mycommand(command)

#    monomer_lengths, monomer_ics, monomer_pwms, monomer_flanked_pwms = get_monomers(factors, results_output, last_iteration_output, model_rows, motif_ending, get_flanks)
    monomer_lengths, monomer_ics, monomer_pwms = get_monomers(factors, results_output, last_iteration_output, model_rows, motif_ending, get_flanks)

    best_cases_titles=[0]*number_of_cobs
    get_dimer_cases(results_output, iterations, last_iteration_output, cob_factors, use_rna, cob_codes, dmin, dmax, cob_tables, cob_ic_tables,
                    cob_length_tables, orients, orient_dict, monomer_pwms, get_flanks, motif_ending, cob_titles, best_cases_titles)
    create_monomer_logos(factors, monomer_lengths, motif_ending, get_flanks)
    visualize_cobs(cobs, cob_codes, cob_tables, dmin, dmax, orients)
    if debug:
        create_graphs(full_output, factors, cobs, cob_codes, name)
        get_start_positions(last_iteration_output, factors)

    # temp=extract_list(r"Background distribution \(intermed\): \[(.+)\]", full_output)
    # bg=[x.split(", ") for x in temp]
    # with open("background.txt", "w") as f:
    #     bg_t=np.array(bg).transpose()
    #     for t in bg_t: 
    #         f.write("%s\n" % ("\t".join(t)))
    # myrun('myspacek40 --logo -paths background.txt background.svg')



    ######################################################################################################
    #
    # Print html
    #



    #logo_files = [ s+".svg" for s in factors ]                             # logo_files=["TEAD4.svg", "ERG.svg"]
    logo_files = [ "monomer.%i.svg" % i for i in range(number_of_factors) ]   # logo_files=["monomer.0.svg", "monomer.1.svg"]
    #cob_files = [ s+".svg" for s in cobs ]                                 # cob_files=["TEAD4-TEAD4.svg", "ERG-ERG.svg",  "TEAD4-ERG.svg"]
    cob_files = [ "cob.%s.svg" % s for s in cob_codes ]                       # cob_files=["cob.0-0.svg", "cob.1-1.svg",  "cob.0-1.svg"]
    cob_links = [ "cob.%s.array.html" % s for s in cob_codes ]             # ["cob.0-0.array.html", "cob.1-1.array.html", ... ]


    best_cases, best_cases_links, best_cases_headers = get_best_cob_cases(cob_codes)

    lambda_table = get_lambda_table(results_output, factors, cobs, cob_codes)



    iterations, maxiter, Lmin, Lmax, lines, epsilon, excluded, command, start_time, version, hostname, threads = get_info(results_output, full_output, cob_codes)

    #Background distribution: [0.30201, 0.294127, 0.202849, 0.201013]
    bg_dist = list(map(float, extract(r"Background distribution: \[(.*)\]", results_output).split(', ')))

    seeds_begin, seeds_end = get_seeds(full_output, number_of_factors)

    monomer_modularity = get_monomer_modularity(full_output)


    runtime = get_run_time(results_output)




    f = open('index.html', 'w')

    f.write("<html>\n")

    f.write("<head>\n")
    javascript='''
      <script type="text/javascript">
        function myclick(e, ending)
        {
        t=e.target;
        t.style.borderColor="red";
        logo_container = t.parentNode.parentNode;
        image_node = logo_container.getElementsByClassName("image")[0];
        anchor_node = logo_container.getElementsByTagName("a")[0];
        if (t.className=="normal") {
          var other=logo_container.getElementsByClassName("rc")[0];
          image_node.src = image_node.src.replace("-rc.svg", ".svg");
          anchor_node.href= anchor_node.href.replace("-rc".concat(ending), ending);
        } else {
          var other=logo_container.getElementsByClassName("normal")[0];
          if (! image_node.src.endsWith("-rc.svg")) {
            image_node.src= image_node.src.replace(".svg", "-rc.svg");
            anchor_node.href= anchor_node.href.replace(ending, "-rc".concat(ending));
          }
        }
        other.style.borderColor="white";
        }
      </script>
    '''
    f.write(javascript)

    f.write('<link rel="stylesheet" href="style.css" type="text/css" />\n')
    f.write("<title>%s - %s</title>\n" % (name, re.sub("^../", "", inputfile)))

    f.write("</head>\n")

    f.write("<body>\n")

    f.write("<h1>MODER - MOtif DEtectoR</h1>\n")
    f.write("<h2>%s - %s</h2>\n" % (name, re.sub("^../", "", inputfile)))


    ###################
    #
    # Print the infobox

    f.write('<div id="info">\n')

    f.write('<div style="float: left; padding: 20px;">\n')
    f.write("<ul>\n")
    #print """<li>Result file: <a href="%s" onclick="window.open('%s', 'newwindow', 'width=300, height=250'); return false;">%s</a></li>""" % (inputfile, inputfile, inputfile)
    f.write("""<li>Result file: <a href="%s">%s</a></li>""" % (inputfile, re.sub("^../", "", inputfile)))

    if Lmin == Lmax:
        L="%s" % Lmin
    else:
        L="%s-%s" % (Lmin, Lmax)
    f.write("<li>Data contains %i sequences of length %s </li>" % (lines, L))
    f.write("<li>Running time was (wall-clock) %s</li>" % runtime)
    if iterations == maxiter:
        f.write("<li>EM-algorithm took <span style='color: red;'>%i iterations</span> (max-iter=%i)</li>" % (iterations, maxiter))
    else:
        f.write("<li>EM-algorithm took %i iterations (max-iter=%i)</li>" % (iterations, maxiter))
    f.write("<li>Convergence criterion cutoff is %g</li>" % epsilon)
    f.write("<li>Excluded cob cases: %i</li>" % excluded)
    f.write("<li>Are monomers learnt modularly: %s</li>" % " ".join(monomer_modularity))
    f.write('<li>Initial and final <a href="seeds.txt">consensus sequences</a> of lengths %s:</li>' % (" ".join([str(len(x)) for x in seeds_begin])))
    f.write("<ul>")
    f.write('<li style="font-family: monospace;">%s</s>' % (" ".join(seeds_begin)))
    f.write('<li style="font-family: monospace;">%s</s>' % (" ".join(seeds_end)))
    f.write("</ul>")
    f.write("<li>Bg: %s</li>" % (" ".join(["%.2f" % x for x in bg_dist])))
    f.write("</ul>")
    f.write("</div>\n")
    f.write('<div id="lambdaTable" style="float: left; padding: 20px;">\n')
    make_table(lambda_table, ["Model", "Lambda"], f)
    f.write("</div>\n")
    f.write("</div>\n")


    ###################
    #
    # Print the factors

    monomer_lambdas=[ y for x,y in lambda_table][0:number_of_factors]
    # These are for the title attribute of the images
    monomer_titles=[ "Lambda: %f&#010;IC: %.2f&#010;Length: %i" % (l,i, length) for l, i, length in zip(monomer_lambdas, monomer_ics, monomer_lengths)]
    f.write('<div id="factors">')
    f.write('<h2>Monomer motifs</h2>')
    make_table_h(f, logo_files, factors, motif_ending, monomer_titles)
    f.write("</div>")


    ###################
    #
    # Print the cob tables and the best case from each cob table

    if number_of_cobs > 0:
        f.write('<div id="cobs">')
        f.write('<h2>COB tables</h2>')
        make_table_v2(f, cob_files, [""]*number_of_cobs, cob_links)
        f.write('<h2>Strongest dimeric case from each cob table</h2>')
        make_table_v3(f, best_cases, best_cases_headers, best_cases_links, ".html", best_cases_titles)
        f.write("</div>")

    ###############################
    #
    # Print the factors with flanks

    if get_flanks:
        f.write('<div id="flankfactors">')
        f.write('<h2>Monomer motifs with flanks</h2>')
        flank_logo_files=["flank-%i.svg" % i for i in range(number_of_factors)]
        make_table_v3(f, flank_logo_files, factors, [x.replace(".svg", ".%s"%motif_ending) for x in flank_logo_files], ".%s"%motif_ending, monomer_titles)
        f.write("</div>")

    ###################
    #
    # Print behaviour as function of iterations

    if debug:
        f.write('<div id="iterations">\n')
        #(files, headers=[], links=[], f, titles=[])
        make_table_v(["distances.svg", "ics.svg", "lambdas.svg", "mll.svg", "parameters.svg"], f=f)
        f.write("</div>\n")

        if start_positions:
            f.write('<div id="startpositions">\n')
            f.write('<img src="start_positions.svg" />\n')
            f.write("</div>\n")

    #print '<div id="background">'
    #print '<p><em>Note.</em> Background model is a multinomial distribution for mononucleotides. In the below logo each position gives the background distribution of the corresponding EM-iteration.</p>'
    #print '<img src="background.svg"\>'
    #print '</div>'

    citation="""Jarkko Toivonen, Teemu Kivioja, Arttu Jolma, Yimeng Yin, Jussi Taipale, Esko Ukkonen (2018) Modular discovery of monomeric and dimeric
    transcription factor binding motifs for large data sets, <i>Nucleic Acids Research</i>,
    Volume 46, Issue 8, 4 May 2018, Pages e44."""

    bibtex=""
    moder_doi="https://doi.org/10.1093/nar/gky027"

    f.write('<div id="programinfo">\n')
    f.write("<ul>\n")
    f.write("<li>Command line was: %s</li>\n" % command)
    f.write("<li>Program started at: %s</li>\n" % start_time)
    f.write("<li>MODER version: %s</li>\n" % version)
    f.write("<li>MODER was run on host: %s</li>\n" % hostname)
    f.write("<li>Number of simultaneous threads: %s</li>\n" % threads)
    f.write("<li>MODER is available from <a href='https://github.com/jttoivon/MODER'>GitHub</a></li>")
    f.write("<li>If you use MODER in your research, please cite: %s <a href='%s'>Link to article.</a>\n</li>" % (citation,moder_doi))
    f.write("</ul>\n")
    f.write("</div>")

    f.write("</body>")

    f.write("</html>")





    ##################################################################################################################
    #
    # Print the cob.x-y.array.html files that contain the observed, expected and deviation logos for all dimeric cases
    #

    # dmin=[0]*number_of_cobs
    # dmax=[0]*number_of_cobs
    # cob_tables=[0]*number_of_cobs

    for i, cob_factor in enumerate(cob_factors):
        number_of_orientations = 1 if use_rna else 3
        if cob_factor[0] != cob_factor[1]:
            number_of_orientations += 1


        link_table = []
        link_expected_table = []
        link_deviation_table = []
        link_flank_table = []
        for row in range(0, number_of_orientations):
            temp_list=[]
            temp_list2=[]
            temp_list3=[]
            temp_list4=[]
            for d in range(dmin[i], dmax[i]+1):
                if cob_tables[i][row,d-dmin[i]] > 0.00000:
                    title=cob_titles[i][row, d-dmin[i]]
                    ending="%s.%s.%i" % (cob_codes[i], orients[row], d)
                    html="three.%s.html" % ending
                    # temp_list.append('<a href="three.%s.html"><img src="observed.%s.svg"\></a>' % (ending, ending))
                    temp_list.append(logo_container(html, "observed.%s.svg" % ending, ".html", title))
                    temp_list2.append(logo_container(html, "expected.%s.svg" % ending, ".html", title))
                    temp_list3.append(logo_container(html, "deviation.%s.svg" % ending, ".html", title))
                    temp_list4.append(logo_container(html, "flank.%s.svg" % ending, ".html", title))
                else:
                    temp_list.append('<p>-</p>')
                    temp_list2.append('<p>-</p>')
                    temp_list3.append('<p>-</p>')
                    temp_list4.append('<p>-</p>')
            link_table.append(temp_list)
            link_expected_table.append(temp_list2)
            link_deviation_table.append(temp_list3)
            link_flank_table.append(temp_list4)
        with open("cob.%s.array.html" % cob_codes[i], "w") as f:
            f.write("<title>%s PPMs - %s</title>\n" % (name, inputfile))
            f.write('<link rel="stylesheet" href="style.css" type="text/css" />\n')
            f.write(javascript)
            f.write('<a href="cob.%s.cob" type="text/plain"><img src="cob.%s.svg"\></a>\n' % (cob_codes[i], cob_codes[i]))

            # Display the logos involved in the cob table
            f.write('<div id="factors">')
            index_set=set(cob_factor)  # If both indices are the same, then return just one
            make_table_h(f, [logo_files[i1] for i1 in index_set], [factors[i2] for i2 in index_set], motif_ending, [monomer_titles[i3] for i3 in index_set])
            f.write("</div>")

            f.write('<h3 class="tableheading">Observed Matrices</h3>\n')
            make_better_table(f, link_table, [""]+list(range(dmin[i], dmax[i]+1)), orients, htmlclass="ppmtable")
            f.write('<h3 class="tableheading">Expected Matrices</h3>\n')
            make_better_table(f, link_expected_table, [""]+list(range(dmin[i], dmax[i]+1)), orients, htmlclass="ppmtable")
            f.write('<h3 class="tableheading">Deviation Matrices</h3>\n')
            make_better_table(f, link_deviation_table, [""]+list(range(dmin[i], dmax[i]+1)), orients, htmlclass="ppmtable")
            if get_flanks:
                for j in range(len(factors)):
                    f.write('<h3 class="tableheading">Flanked Monomer Matrix %i</h3>\n' % j)

                    f.write(logo_container("flank-%i.%s" % (j,motif_ending), "flank-%i.svg" % j, ".%s"%motif_ending, monomer_titles[j]))
                    #f.write('<img src="flank-%i.svg"\>\n' % j)
                f.write('<h3 class="tableheading">Flanked Dimer Matrices</h3>\n')
                make_better_table(f, link_flank_table, [""]+list(range(dmin[i], dmax[i]+1)), orients, htmlclass="ppmtable")

    #myrun("cp %s/style.css ." % execdir)

    with open("monomer_weights.txt", "w") as f:
        f.write("%s\n" % (",".join(factors)))
        f.write("%s\n" % (",".join(map(str, monomer_lambdas))))

    with open("style.css", "w") as f:
        f.write(css)

    print("The report is in file %s" % reportfile)

if __name__ == '__main__':
    main()
