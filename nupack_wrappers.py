"""
Wrapper functions around nupack 
"""


import numpy as np
import subprocess
import os
import re
from glob import glob


#############################
### Single strand folding ###
#############################

def pfunc_single(seq, material='dna', temp=37, sodium=0.15, magnesium=0.003):
    """
    Compute the partition function for a single nucleic acid strand
    Return both the ensemble energy and the partition function
    Inputs:
        seq = nucleic acid sequence
        material = nupack code for parameter set (default is 'dna')
        temp = temperature in C
        sodium = sodium concentration in M
        magnesium = magnesium concentration in M
    """
    # First, write a temporary file containing the necessary information for nupack
    temp_prefix = '_temp'
    with open(temp_prefix + '.in', 'w') as f:
        f.write(seq)
        f.write('\n')
    
    # Now, run nupack pfunc and collect output
    pfunc_out = subprocess.check_output(['pfunc', '-material', material, '-sodium', str(sodium), '-magnesium', str(magnesium), 
                     '-T', str(temp), temp_prefix])
    
    # Remove temporary files
    for f in glob(temp_prefix + '*'):
        os.remove(f)
    
    # Parse output to collect partition function and ensemble energy
    lines = pfunc_out.split('\n')
    results = []
    for l in lines:
        if (len(l) == 0) or (l[0] == '%'):
            continue
        results.append(float(l))
    # Energy (kcal / mol) reported first, then partition function
    return results[0], results[1]



def mfe_single(seq, material='dna', temp=37, sodium=0.15, magnesium=0.003):
    """
    Compute the mfe energy of a list of nucleic acid strands
    Inputs:
        seq = sequence to calculate
        material = nupack code for parameter set (default is 'dna')
        temp = temperature in C
        sodium = sodium concentration in M
        magnesium = magnesium concentration in M
    """
    # First, write a temporary file containing the necessary information for nupack
    temp_prefix = '_temp'
    with open(temp_prefix + '.in', 'w') as f:
        f.write(seq)
        f.write('\n')
    
    # Now, run nupack pfunc and collect output
    subprocess.call(['mfe', '-material', material, '-sodium', str(sodium), '-magnesium', str(magnesium), 
                     '-T', str(temp), temp_prefix])
    
    
    # First, just read file into list of lines
    with open(temp_prefix + '.mfe', 'r') as f:
        lines = [line.strip() for line in f]
            
    # Remove temporary files
    for f in glob(temp_prefix + '*'):
        os.remove(f)
    
    # mfe output file is formatted as such:
    # First several lines are comments denoting settings (all begin with a single '%')
    # blank line following header (no '%')
    # result bookend ('% %+ %')
    # Each result contains the same first three lines:
    # 1 : number of bases in ordered complex
    # 2 : minimum free energy (kcal/mol)
    # 3 : dot-parens-dot representation of MFE structure
    # 4 - n: MFE structure in pair list notation

    # Parse mfe output file to mfe energy and dot-bracket notation
    div_pat = re.compile(r'\% \%+ \%')
        
    # Get indicees of result bookends
    idxs = [i for i, line in enumerate(lines) if div_pat.match(line)]
    result = lines[idxs[0]+1:idxs[1]]
    mfe = float(result[1])
    dot_bracket = result[2]
        
    return mfe, dot_bracket



################################
### two-strand hybridization ###
################################


def pfunc_multi(seq_list, material='dna', temp=37, sodium=0.15, magnesium=0.003):
    """
    Compute the partition function for a list of nucleic acid strands.
    Return both the ensemble energy and the partition function
    Inputs:
        seq_list = list of sequences (must be > 1)
        material = nupack code for parameter set (default is 'dna')
        temp = temperature in C
        sodium = sodium concentration in M
        magnesium = magnesium concentration in M
    """
    if len(seq_list) <= 1:
        print "Error: must have >= 1 sequences in seq_list"
        return None
    
    # First, write a temporary file containing the necessary information for nupack
    temp_prefix = '_temp'
    with open(temp_prefix + '.in', 'w') as f:
        f.write("{}\n".format(len(seq_list)))
        for s in seq_list:
            f.write(s + '\n')
        f.write(' '.join([str(x) for x in range(1, len(seq_list) + 1, 1)]))
    
    # Now, run nupack pfunc and collect output
    pfunc_out = subprocess.check_output(['pfunc', '-multi', '-material', material, '-sodium', str(sodium), '-magnesium', str(magnesium), 
                     '-T', str(temp), temp_prefix])
    
    # Remove temporary files
    for f in glob(temp_prefix + '*'):
        os.remove(f)
    
    # Parse output to collect partition function and ensemble energy
    lines = pfunc_out.split('\n')
    results = []
    for l in lines:
        if (len(l) == 0) or (l[0] == '%'):
            continue
        results.append(float(l))
    
    # Energy (kcal / mol) reported first, then partition function
    return results[0], results[1]




def mfe_multi(seq_list, material='dna', temp=37, sodium=0.15, magnesium=0.003):
    """
    Compute the mfe energy of a list of nucleic acid strands
    Inputs:
        seq_list = list of sequences (must be > 1)
        material = nupack code for parameter set (default is 'dna')
        temp = temperature in C
        sodium = sodium concentration in M
        magnesium = magnesium concentration in M
    """
    if len(seq_list) <= 1:
        print "Error: must have >= 1 sequences in seq_list"
        return None
    
    # First, write a temporary file containing the necessary information for nupack
    temp_prefix = '_temp'
    with open(temp_prefix + '.in', 'w') as f:
        f.write("{}\n".format(len(seq_list)))
        for s in seq_list:
            f.write(s + '\n')
        f.write(' '.join([str(x) for x in range(1, len(seq_list) + 1, 1)]))
    
    # Now, run nupack pfunc and collect output
    subprocess.call(['mfe', '-multi', '-material', material, '-sodium', str(sodium), '-magnesium', str(magnesium), 
                     '-T', str(temp), temp_prefix])
    
    
    # First, just read file into list of lines
    with open(temp_prefix + '.mfe', 'r') as f:
        lines = [line.strip() for line in f]
            
    # Remove temporary files
    for f in glob(temp_prefix + '*'):
        os.remove(f)
    
    # mfe output file is formatted as such:
    # First several lines are comments denoting settings (all begin with a single '%')
    # blank line following header (no '%')
    # result bookend ('% %+ %')
    # Each result contains the same first three lines:
    # 1 : number of bases in ordered complex
    # 2 : minimum free energy (kcal/mol)
    # 3 : dot-parens-dot representation of MFE structure
    # 4 - n: MFE structure in pair list notation

    # Parse mfe output file to mfe energy and dot-bracket notation
    div_pat = re.compile(r'\% \%+ \%')
        
    # Get indicees of result bookends
    idxs = [i for i, line in enumerate(lines) if div_pat.match(line)]
    result = lines[idxs[0]+1:idxs[1]]
    mfe = float(result[1])
    dot_bracket = result[2]
        
    return mfe, dot_bracket

