#!/usr/bin/python

import sys

input_file_name = sys.argv[1]

times = dict()
algos = []
file_names = []

def addTime(times, algo, procs, fname, time):
    if fname not in file_names: 
        file_names.append(fname)
    if algo not in times:
        times[algo] = dict()
        algos.append(algo)
    if fname not in times[algo]:
        times[algo][fname] = dict()
    if procs not in times[algo][fname]:
        times[algo][fname][procs] = time 
    if float(time) < float(times[algo][fname][procs]): 
        times[algo][fname][procs] = time

def getMin(ls):
    res = ls[0]
    for item in ls:
        if float(item) < float(res):
            res = item
    return res
def makeShort(algo):
    return algo.replace("parallel", "par").replace("serial","ser").replace("Divsufsort", "Dss")

with  open(input_file_name, "r") as ins:
    for line in ins:
        tokens = line.split(" ")
        if len(tokens) == 2:
            try:
                float(tokens[1])
                addTime(times, cur_algorithm, cur_num_procs, cur_file, tokens[1].strip())
            except:
                cur_algorithm = tokens[0]
                cur_num_procs = 1
                cur_file = tokens[1].strip()
        if len(tokens) == 3:
            cur_algorithm = tokens[0]
            cur_num_procs = int(tokens[1])
            cur_file = tokens[2].strip()
            
line = "file & "
serial_algo = ["serialDivsufsort", "serialKS"]
algos.remove("parallelRangeLight")
for algo in algos:
    if algo in serial_algo:
        line += makeShort(algo) + " & "
    else:
        line += makeShort(algo) + "($T_1$) & " + makeShort(algo) + "($T_{64}$) & "
print line[:-2] + "\\\\"
print "\hline"
for fname in file_names:
    line = fname + " & "
    for algo in algos:
        if algo in serial_algo:
            line += times[algo][fname][1] + " & "
        else:
            line += times[algo][fname][1] + " & " + getMin(times[algo][fname].values()) + " & "
    print line[:-2] + "\\\\"
