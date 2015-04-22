#!/usr/bin/python

import sys
import collections

def read_file(f):
    times = collections.OrderedDict()
    memory = collections.OrderedDict()
    content = f.readlines()
    for line in content:
        if line[0] == '#': #ignore marked lines
            continue
        tokens = map(str.strip,  line.split(":"))
        times[tokens[1]] = tokens[3].replace("'","")
        memory[tokens[1]] = tokens[4].replace(".","")
    return (f.name, times,memory)
files = [open(file_name,"r") for file_name in sys.argv[2:]]
data = map(read_file,files)
#first parameter m/t/a (memory, times, all)
mem = tim = False
if sys.argv[1] == 'm':
    mem = True
if sys.argv[1] == 't':
    tim = True
if sys.argv[1] == 'a':
    tim = mem = True

# print first line
sys.stdout.write("File name")
for (n,x,y) in data:
    if tim: 
        sys.stdout.write("," +  n + "(t)")
    if mem:
        sys.stdout.write("," + n + "(m)")
sys.stdout.write("\n")
#print the following lines
(n,x,y) = data[0]
for name in x.keys():
    line = ""
    for (n,times,memory) in data:
        if tim:
            line += "," + times[name] 
        if mem:
            line += ","  + memory[name] 
    sys.stdout.write(name + line + "\n")
