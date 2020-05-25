# -*- coding: utf-8 -*-
"""
Script to convert a Position specific scoring matrices (PSSMs) file
from RegulonDB to the minimal MEME format with position weight matrix
(PWM) data.
"""
import re
import pandas as pd

# Open PSSM file downloaded from RegulonDB
PSSM_file = open('pwm.txt', mode = 'r')
PSSM_data = PSSM_file.readlines()
PSSM_file.close()

# Splits list into list of lists, each TF getting its own list
split_PSSM = []
for item in PSSM_data:
    if re.match("Transcription Factor ID", item) or len(split_PSSM) == 0:
        split_PSSM.append([item])
    else:
        split_PSSM[-1].append(item)
        
split_PSSM = split_PSSM[1:len(PSSM_data)]

# Gets PCM data from each TF sublist as panda dataframe
pcm_data = []
for item in split_PSSM:
    cur_pcm_data = item[-7:-3]
    cur_pcm_lol = []
    for line in cur_pcm_data:
        cur_line = line.split(sep = '\t')
        cur_line = cur_line[1:len(cur_line)]
        cur_pcm_lol.append(cur_line)
    cur_pcm_dataframe = pd.DataFrame(cur_pcm_lol).transpose()
    cur_pcm_dataframe.columns = ['A', 'C', 'G', 'T']
    cur_pcm_dataframe = cur_pcm_dataframe.astype('float')
    pcm_data.append(cur_pcm_dataframe)

# Convert Position Count Data to Position Weight Data
pwm_data = []
pcm_copy = pcm_data.copy()
for item in pcm_copy:
    item['counts'] = item.sum(axis = 1)
    item['A'] = item['A']/item['counts']
    item['C'] = item['C']/item['counts']
    item['T'] = item['T']/item['counts']
    item['G'] = item['G']/item['counts']
    pwm_data.append(item[['A', 'C', 'G', 'T']])
    
# Generates MEME compatable file from PWM
header_text = 'MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n'
string_out = header_text
for i in range(0,len(split_PSSM)):
    TF_name = split_PSSM[i][1].split(':')[1]
    cur_pwm = pwm_data[i].to_string(index = False, header= False)
    string_out = string_out + 'MOTIF' + TF_name + 'letter-probability matrix:' + '\n' + cur_pwm + '\n\n'
    

text_file = open("Minimal_MEME.txt", "w")
text_file.write(string_out)
text_file.close()
