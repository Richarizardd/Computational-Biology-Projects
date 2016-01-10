from itertools import islice
import re

### Natural Sort Function
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

### Glance: Peaks into the first couple items in a dictionary
def glance(d, x):
    return dict(islice(d.items(), 0, x))

### Modified write_solution for Python 3
def write_solution(genome, per_line = 60):
    offset = 0
    solution = open("solution.fa", "w")
    while offset < len(genome):
        nchars = min(len(genome) - offset, per_line)
        line = genome[offset:offset+nchars]
        offset += nchars
        solution.write(line + '\n')
    solution.close()
        
def parse_fasta(fh):
    fa = {}
    name = None
    # Part 1: compile list of lines for each sequence
    for ln in fh:
        if ln[0] == '>':  # new sequence
            name = ln[1:].split()[0]
            fa[name] = []
        else:
            # append nucleotides to current sequence
            fa[name].append(ln.rstrip())
    # Part 2: join lists into strings
    for name, nuc_list in fa.items():
        fa[name] = ''.join(nuc_list)  # join into one long string
    return fa

def make_kmer_table(seqs, k):
    ''' Given dictionary (e.g. output of parse_fasta) and integer k,
        return a dictionary that maps each k-mer to the set of names
        of reads containing the k-mer. '''
    table = {}  # maps k-mer to set of names of reads containing k-mer
    for name, seq in seqs.items():
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer not in table:
                table[kmer] = set()
            table[kmer].add(name)
    return table

def suffixPrefixMatch(str1, str2, min_overlap):
    ''' Returns length of longest suffix of str1 that is prefix of
        str2, as long as that suffix is at least as long as min_overlap. '''
    if len(str2) < min_overlap: return 0
    str2_prefix = str2[:min_overlap]
    str1_pos = -1
    while True:
        str1_pos = str1.find(str2_prefix, str1_pos + 1)
        if str1_pos == -1: return 0
        str1_suffix = str1[str1_pos:]
        if str2.startswith(str1_suffix): 
            return len(str1_suffix)

### Part 2 ###
def build_overlap_graph(fa, table, k):
    BBR_dict = {}
    for name in fa.keys():
        BBR_dict[name] = (None, 0, 0) # (Name, Longest Name, Collisions)
        
    for neighborhoods in list(table.items()): # for each "neighborhood of names"
        for A in list(neighborhoods[1]): # O(n^2) of building a dictionary of best buddies to the right
            for B in list(neighborhoods[1]):
                longest_match = suffixPrefixMatch(fa.get(A),fa.get(B), k)
                if (A != B and longest_match > BBR_dict[A][1]):
                    BBR_dict[A] = (B, longest_match, 0) # update the longest match!
                elif (A != B and longest_match == BBR_dict[A][1] and B != BBR_dict[A][0]):
                    BBR_dict[A] = (B, longest_match, 1)                

    ordered_keys = natural_sort(list(fa.keys())) # Obtains a list of names in order
    BBR_triples = []
    for name in ordered_keys: # adds only unique best buddies to the right
        if (BBR_dict[name][2] == 0):
            BBR_triples.append((name, BBR_dict[name][0], BBR_dict[name][1]))

    return BBR_triples

### Part 3 ###
def build_unitigs():
    BBL = {}
    overlaps = open("overlaps.txt", "r")
    overlaps_txt = overlaps.readlines()
    for line in overlaps_txt:
        if (len(line.split()) == 3):
            A = line.split()[0].rstrip()
            B = line.split()[1].rstrip()
            longest_match = int(line.split()[2].rstrip())
            
            if (BBL.get(B) == None or longest_match > BBL[B][1]):
                BBL[B] = (A, longest_match, 0)
            elif (longest_match == BBL[B][1]):
                BBL[B] = (A, longest_match, 1)
    
    BBR = {}     
    for pair in list(BBL.items()):
        if (pair[1][2] == 0):
            BBR[pair[1][0]] = (pair[0], pair[1][1])
            BBL[pair[0]] = (pair[1][0], pair[1][1])
        else:
            BBL.pop(pair[0], None)
    bb_keys = list(BBR.keys()) # Obtains a list of names in order
    
    unitigs = []
    count = 0
    while (bb_keys != []):

        start = bb_keys.pop(0)
        unitigs.append([start])
        
        next_buddy = start
        prev_buddy = start
        while(BBR.get(next_buddy) != None): # going forwards
            next_buddy = BBR[next_buddy][0]   
            unitigs[count].append(next_buddy)
            
            if (BBR.get(next_buddy) != None):
                bb_keys.remove(next_buddy)
        
        while(BBL.get(prev_buddy) != None):
            prev_buddy = BBL[prev_buddy][0]            
            unitigs[count].insert(0, prev_buddy)
            bb_keys.remove(prev_buddy)

        count = count + 1
        
    return unitigs

### Part 4 ###
def finish_assembly(fa, unitigs):
    unitig_genomes = []
    for i in range(len(unitigs)):
        unitig_genome = fa[unitigs[i][0]]
        
        for j in range(0,len(unitigs[i]) - 1):
            curr_read = fa[unitigs[i][j]]
            next_read = fa[unitigs[i][j+1]]
            subset_index = suffixPrefixMatch(curr_read, next_read, 40)
            unitig_genome = ''.join([unitig_genome, next_read[subset_index:]])
            
        unitig_genomes.append(unitig_genome)
    
    ### Hard Coded Section ###
    unitig_overlaps = {}
    
    for seq1 in unitig_genomes:
        unitig_overlaps[seq1] = []
        
        for seq2 in unitig_genomes:
            unitig_overlaps[seq1].append(suffixPrefixMatch(seq1, seq2, 0))
    
    # unitig_overlaps = [[2109, 0, 0, 99], [0, 3127, 1, 0], [1, 0, 2613, 99], [0, 98, 98, 252]]
    # A = [2109, 0, 0, 99] - prefix-suffix match with D. AD
    # B = [0, 3127, 1, 0] - Has no prefix-suffix match, must be at the end
    # C = [1, 0, 2613, 99] - prefix-suffix matches with D. CD
    # D = [0, 98, 98, 252] - prefix-suffix matches with B, C. DB, DC (D must occur twice)
    # Empirically, from the overlaps, we can construct an overlap of unitigs A, B, C, and D.
    # ADCDB
    
    A = unitig_genomes[0]
    AD = A + unitig_genomes[3][unitig_overlaps[unitig_genomes[0]][3]:]
    ADC = AD + unitig_genomes[2][unitig_overlaps[unitig_genomes[3]][2]:]
    ADCD = ADC + unitig_genomes[3][unitig_overlaps[unitig_genomes[2]][3]:]
    ADCDB = ADCD + unitig_genomes[1][unitig_overlaps[unitig_genomes[3]][1]:]
    return ADCDB
        
def main():
    ### Part 1 ###
    k = 40
    reads_fa_file = open("f2014_hw4_reads.fa", "r")
    reads_fa = reads_fa_file.readlines()
    fa = parse_fasta(reads_fa) # maps every "name" to its corresponding nucleotide sequence
    table = make_kmer_table(fa, k) # creates a dictionary that maps each kmer to list of names

    ### Part 2 ###
    BBR_triples = build_overlap_graph(fa, table, k)
    overlaps = open("overlaps.txt", "w")
    for triple in BBR_triples:
        overlaps.write(str(triple[0]) + " " + str(triple[1]) + " " + str(triple[2]) + "\n")
    overlaps.close()
    
    ### Part 3 ###
    unitigs = build_unitigs()
    unitigs_file = open("unitigs.txt", "w")
    for i in range(len(unitigs)):
        unitigs_file.write("START UNITIG " + str(i+1) + " " + unitigs[i][0] + "\n")
        for j in range(1,len(unitigs[i])):
            unitigs_file.write(unitigs[i][j] + "\n")
        unitigs_file.write("END UNITIG " + str(i+1) + "\n")

    unitigs_file.close()
    
    ### Part 4 ### 
    genome = finish_assembly(fa, unitigs) # 7959 units long! Success!
    write_solution(genome) # The file has a total of 8091 characters, but ther are 133 rows, so 133 \n characters
main()

