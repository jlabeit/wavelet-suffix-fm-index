import subprocess
import os

integer = ['rnd-8', 'rnd-12', 'rnd-16', 'rnd-20']
artificial = ['fib41', 'rs.13', 'tm29']
real = ['Escherichia_Coli', 'cere', 'coreutils', 'einstein.de.txt',
        'influenza', 'kernel', 'para', 'world_leaders']
pseudo= ['dblp.xml.00001.1', 'dblp.xml.00001.2', 'dblp.xml.0001.1', 'dblp.xml.0001.2',
        'dna.001.1', 'english.001.2', 'proteins.001.1', 'sources.001.2']
classic = ['sources', 'pitches', 'proteins', 'dna', 'english.1024MB', 'dblp.xml']
special = ['aaa', 'abab', 'aabbaabb']

INT = True
if INT:
    inputs = [integer]
    result_file = open('results.int', 'w')
else:
    inputs = [special, real, pseudo, artificial, classic]
    result_file = open('results', 'w')
threads = [1, 2, 4, 8, 12, 16, 24, 32, 40, 48, 56, 64]

salgos = ['serialWT', 'sdslWT']
palgos = ['ddWT', 'levelWT', 'recWT']

def get_command(algo, f, p):
    executable = './' + algo + '/WT'
    filedir = '../input/' + f + '.bwt'
    if (p == 1):
        return [executable, filedir]
    else:
        return ['numactl', '--interleave=all', executable, filedir]

def get_env(algo, p):
    result = os.environ.copy()
    result["CILK_NWORKERS"] = str(p)
    return result

def parse_output(s):
    times = ''
    mems = ''
    lines = s.split('\n')
    for l in lines:
        tokens = l.split(': ')
        if len(tokens) > 1:
            if tokens[0].strip() == 'PBBS-time':
                times = tokens[1] 
            elif tokens[0].strip() == 'Peak-memory':
                mems = tokens[1]
            else:
                print "error ", tokens[0]
    return (times, mems)

def write_output(s, algo, f, p):
    (time,mem) = parse_output(s)
    result_file.write(algo + ", " + f + ", " + str(p) + ", " + mem + ", " + time + '\n')
    result_file.flush()

    
def run_test(algo, f, p):
    commands = get_command(algo, f, p)
    cur_env = get_env(algo, p)
    print 'Running: ', commands
    proc = subprocess.Popen(commands, stdout=subprocess.PIPE, env=cur_env)
    write_output(proc.communicate()[0], algo, f, p)
        

for files in inputs:
    for f in files:
        # Parallel
        for algo in palgos:
            for p in threads:
                run_test(algo, f, p)
        # Sequential
        for algo in salgos:
            run_test(algo, f, 1)
