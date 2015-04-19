#!/usr/bin/python

import sys

def read_file(f):
    times = dict()
    memory = dict()
    content = f.readlines()
    for line in content:
        if line[0] == '#': #ignore marked lines
            continue
        tokens = map(str.strip,  line.split(":"))
        times[tokens[1]] = tokens[3].replace("'","")
        memory[tokens[1]] = tokens[4].replace(".","")
    return (f.name, times,memory)
files = [open(file_name,"r") for file_name in sys.argv[1:]]
data = map(read_file,files)

# print first line
sys.stdout.write("File name")
for (n,x,y) in data:
    sys.stdout.write("," +  n + "(t)" + "," + n + "(m)")
sys.stdout.write("\n")
#print the following lines
(n,x,y) = data[0]
for name in x.keys():
    line = ""
    for (n,times,memory) in data:
        line += "," + times[name] + ","  + memory[name] 
    sys.stdout.write(name + line + "\n")
