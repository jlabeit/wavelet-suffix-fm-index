#!/usr/bin/python

import sys

input_file_name = sys.argv[1]

times = dict()
memory = dict()
algos = []
file_names = []

def addMemory(memory, algo, procs, fname, mem):
    if fname not in file_names: 
        file_names.append(fname)
    if algo not in memory:
        memory[algo] = dict()
        if algo not in algos:
            algos.append(algo)
    if fname not in memory[algo]:
        memory[algo][fname] = dict()
    if procs not in memory[algo][fname]:
        memory[algo][fname][procs] = mem 
    if float(mem) < float(memory[algo][fname][procs]): 
        memory[algo][fname][procs] = mem

   
def addTime(times, algo, procs, fname, time):
    if fname not in file_names: 
        file_names.append(fname)
    if algo not in times:
        times[algo] = dict()
        if algo not in algos:
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
                if tokens[0] == "Peak-memory:":
                    addMemory(memory, cur_algorithm, cur_num_procs, cur_file, tokens[1].strip())
                else:
                    addTime(times, cur_algorithm, cur_num_procs, cur_file, tokens[1].strip())
            except:
                cur_algorithm = tokens[0]
                cur_num_procs = 1
                cur_file = tokens[1].strip()
        if len(tokens) == 3:
            cur_algorithm = tokens[0]
            cur_num_procs = int(tokens[1])
            cur_file = tokens[2].strip()
            
write_memory = True
line = "file & "
serial_algo = ["serialDivsufsort", "serialKS"]
algos.remove("parallelRangeLight")
for algo in algos:
    if algo in serial_algo:
        line += makeShort(algo) + " & "
    else:
        if write_memory:
            line += makeShort(algo) + " & "
        else:
            line += makeShort(algo) + "($T_1$) & " + makeShort(algo) + "($T_{40}$) & "
print line[:-2] + "\\\\"
print "\\midrule"
for fname in file_names:
    line = fname + " & "
    for algo in algos:
        if algo in serial_algo:
            if write_memory:
                line += memory[algo][fname][1] + " & "
            else:
                line += times[algo][fname][1] + " & "
        else:
            if write_memory:
                line += getMin(memory[algo][fname].values()) + " & "
            else:
                line += times[algo][fname][1] + " & " + getMin(times[algo][fname].values()) + " & "
    print line[:-2] + "\\\\"
print "\\bottomrule"
