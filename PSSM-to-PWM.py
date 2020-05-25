# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import re
import numpy
import pandas as pd

# Open PWM file downloaded from RegulonDB
pwm_file = open('pwm.txt', mode = 'r')
pwm_data = pwm_file.readlines()

# Splits list into list of lists, each TF getting its own list
split_pwms = []
for item in pwm_data:
    if re.match("Transcription Factor ID", item) or len(split_pwms) == 0:
        split_pwms.append([item])
    else:
        split_pwms[-1].append(item)
        
split_pwms = split_pwms[1:len(pwm_data)]

# Gets PCM data from each TF sublist as panda dataframe
pcm_data = []
for item in split_pwms:
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
