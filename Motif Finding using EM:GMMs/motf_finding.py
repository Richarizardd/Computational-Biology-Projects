# # HW4 Extra Credit: Motif Finding
# #### EN.600.438 Comptutational Genomics: Data Analysis
# #### Richard Chen
# #### Worked with: Steven Chen, Yunfan Fan

# ### 0. Dependencies
import numpy as np
import re
import math


# ### 1. E Step: Estimate y(Zij) from P
def E_Step(W, read_list, read_length, PFM, nucleotides_dict):
    for pos in range(read_length):
        all_curr_pos = "".join([line[pos] for line in read_list])
        for base in range(len(nucleotides_dict.keys())):
            PFM[base][pos] = all_curr_pos.count(nucleotides_dict.keys()[base])
    PFM = PFM/len(read_list)
    PFM_init = PFM[0:4,:W]
    return PFM_init


# ### 2. M Step: Estimate Pck and Bc using y(Z)
def M_Step(W, read_list, read_length, PFM, nucleotides_dict, nucleotides_prob, Pck):
    gamma = np.zeros((len(read_list),read_length))
    for ind in range(len(read_list)):
        for pos in range(read_length-W+1):
            curr_motif = list(read_list[ind][pos:pos+W]); else_seq = list(read_list[ind][:0] + read_list[ind][W:])
            curr_motif_prob = reduce(lambda x,y: x*y, [Pck[:,i][[nucleotides_dict[base] for base in list(curr_motif)][i]] for i in range(W)])
            else_seq_prob = reduce(lambda x,y: x*y, [nucleotides_prob[base] for base in else_seq])
            gamma[ind,pos] = curr_motif_prob*else_seq_prob
    rowsums = np.sum(gamma, axis = 1)
    gamma = gamma/rowsums[:,None]

    nck = np.zeros((4,W))
    for i in range(len(nucleotides_dict.keys())):
        for ind in range(len(read_list)):
            indices = [m.start() for m in re.finditer(nucleotides_dict.keys()[i], read_list[ind][:read_length-W+1])]
            for j in range(W):
                nck[i][j] = nck[i][j] + sum([gamma[ind][index-j] for index in indices])
    Pck = (nck + 1)/(np.sum(nck, axis = 0) + 4)

    B = []
    for ind in range(len(nucleotides_dict.keys())):
        mc = 0
        for read in read_list:
            mc += read.count(nucleotides_dict.keys()[ind])
        gc = mc - sum(nck[ind]) + 1
        B.append(gc)
    denom = sum(B)
    B = B/denom
    
    return gamma, Pck, B


# ### 3. Compute log likelihood of the data
def loglikelihood(gamma):
    return sum([math.log10(p) for p in np.sum(gamma, axis = 1)])


# ### 4a. Implement a function findmotif(file name, motif width, iterations) using the EM algorithm described above.
def findmotif(file_name, motif_width, iterations):
    read_list = [line.rstrip() for line in open(file_name).readlines()]
    read_length = len(read_list[0]); W = motif_width; iterations = 100; loglik_prev = 0; itr = 1; conv = False
    
    PFM = np.zeros((4,read_length))
    nucleotides_dict = {'A':0, 'C':1, 'T':2, 'G':3}
    nucleotides_prob = {'A':0.25, 'C':0.25, 'T':0.25, 'G':0.25}
    
    Pck = E_Step(W, read_list, read_length, PFM, nucleotides_dict)

    while (itr < 100 and conv == False):
        gamma, Pck, B = M_Step(W, read_list, read_length, PFM, nucleotides_dict, nucleotides_prob, Pck)
        
        nucleotides_prob['A'] = B[0]
        nucleotides_prob['C'] = B[1]
        nucleotides_prob['T'] = B[2]
        nucleotides_prob['G'] = B[3]

        loglik_curr = loglikelihood(gamma)
        if (abs(loglik_curr) - loglik_prev < 0.001):
            conv = True
        loglik_prev = loglik_curr
        itr += 1
    return Pck


# ### 4b. Test your function using sequences in file shortmotif.txt.
Pck = findmotif('shortmotif.txt', 5, 100)
np.savetxt("shortmotif_PWM.csv", Pck, delimiter="\t")
print(Pck)
max_indices = [list(Pck[:,i]).index(max(Pck[:,i])) for i in range(5)]
motif = "".join([['A', 'C', 'T', 'G'][i] for i in max_indices])
print(motif)


# ### 4c. Test your function using sequences in file longmotif.txt.
Pck = findmotif('longmotif.txt', 8, 50)
np.savetxt("longmotif_PWM.csv", Pck, delimiter="\t")
print(Pck)
max_indices = [list(Pck[:,i]).index(max(Pck[:,i])) for i in range(8)]
motif = "".join([['A', 'C', 'T', 'G'][i] for i in max_indices])
print(motif)

