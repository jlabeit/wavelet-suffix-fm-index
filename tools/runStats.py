import subprocess
import os

artificial = ('artificial repetitive', ['fib41', 'rs.13', 'tm29'])
real = ('real repetitive', ['Escherichia_Coli', 'cere', 'coreutils', 'einstein.de.txt',
        'influenza', 'kernel', 'para', 'world_leaders'])
pseudo= ('pseudo repetitive', ['dblp.xml.00001.1', 'dblp.xml.00001.2', 'dblp.xml.0001.1',
    'dblp.xml.0001.2', 'dna.001.1', 'english.001.2', 'proteins.001.1', 'sources.001.2'])
classic = ('non-repetitive', ['sources', 'pitches', 'proteins', 'dna', 'english.1024MB', 'dblp.xml'])
special = ('special cases', ['aaa', 'abab', 'aabbaabb'])

inputs = [classic, real, artificial, pseudo, special]
result_file = open('stats.tex', 'w')

def get_command(f):
    executable = './STATS'
    return [executable, f, '../input/']

def write_output(s, f):
    result_file.write(s)
    result_file.flush()
def write_label(l):
    result_file.write('\\multicolumn{2}{r}{\\textit{' + label + '}} \\\\\n')
    
def run_stat(f):
    commands = get_command(f)
    print 'Running: ', commands
    proc = subprocess.Popen(commands, stdout=subprocess.PIPE)
    write_output(proc.communicate()[0], f)

result_file.write('\\begin{tabular}{llll}\n')
result_file.write('\\toprule\n')
result_file.write('Collection & $n$ & $\sigma$ & avg. LCP \\[2ex]\n')
result_file.write('\hline\n')
for (label, files) in inputs:
    write_label(label)
    for f in files:
        run_stat(f)
result_file.write('\\bottomruler\n')
result_file.write('\\end{tabular}\n')
