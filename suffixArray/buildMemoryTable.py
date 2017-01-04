import csv

sizes = {'abab': 104857600, 'kernel': 257961544, 'dblp.xml.0001.1': 104857600, 'influenza': 154808555, 'dblp.xml.0001.2': 104857600, 'sources.001.2': 104857600, 'proteins.001.1': 104857600, 'sources': 210866603, 'einstein.de.txt': 92758441, 'english.001.2': 104857600, 'dna': 403927746, 'Escherichia_Coli': 112689515, 'world_leaders': 46968181, 'coreutils': 205281760, 'fib41': 267914296, 'english.1024MB': 1073741816, 'dblp.xml': 296135874, 'para': 429265758, 'aabbaabb': 104857600, 'pitches': 55814376, 'dna.001.1': 104857600, 'tm29': 268435456, 'rs.13': 216747218, 'cere': 461286644, 'aaa': 104857600, 'proteins': 1184051855, 'dblp.xml.00001.2': 104857600, 'dblp.xml.00001.1': 104857600}

# Read result data.

def parseRow(row):
    return [row[0].strip(), row[1].strip(), int(row[2]), float(row[3]), float(row[4])]

data = []
with open('results', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in spamreader:
        data.append(parseRow(row))

artificial = ('artificial repetitive', ['fib41', 'rs.13', 'tm29'])
real = ('real repetitive', ['Escherichia_Coli', 'cere', 'coreutils', 'einstein.de.txt',
        'influenza', 'kernel', 'para', 'world_leaders'])
pseudo= ('pseudo repetitive', ['dblp.xml.00001.1', 'dblp.xml.00001.2', 'dblp.xml.0001.1',
    'dblp.xml.0001.2', 'dna.001.1', 'english.001.2', 'proteins.001.1', 'sources.001.2'])
classic = ('non-repetitive', ['sources', 'pitches', 'proteins', 'dna', 'english.1024MB', 'dblp.xml'])
special = ('special cases', ['aaa', 'abab', 'aabbaabb'])

inputs = [classic, real, artificial, pseudo, special]
#threads = [1, 2, 4, 8, 12, 16, 24, 32, 40, 48, 56, 64]
threads = [1, 8, 16, 32, 64]

salgos = ['divsufsort']
palgos = ['parallelKS', 'parallelRange', 'parallelDivsufsort', 'pScan']


result_file = open('memorySuffix.tex', 'w')

def getMemory(a, b, c):
    size = 1.0
    for row in data:
        if row[0] == a and row[1] == b and row[2] == c:
            size = row[3]
    return size * 1024 * 1024 / sizes[b];

def write_label(l):
    result_file.write('\\multicolumn{2}{l}{\\textit{' + label + '}} \\\\\n')

def getStrMemory(f, algo, smallest):
    if (algo in salgos):
        return '%.2f '% getMemory(algo, f, 1)
    else:
        t64 = getMemory(algo, f, 64)
        if (smallest):
            return '\\textbf{%.2f} &' % t64 
        return '%.2f & ' % t64 

def getSmallest(f):
    mini = 100000000.0;
    result = ''
    for algo in palgos:
        cur = getMemory(algo, f, 64)
        if mini > cur:
            mini = cur
            result = algo
    return result

def write_result(f):
    result_file.write(f.replace('_', '\\_') + ' & ')
    smallest = getSmallest(f)
    for algo in palgos:
        result_file.write(getStrMemory(f, algo, smallest == algo))
    for algo in salgos:
        result_file.write(getStrMemory(f, algo, smallest == algo))
    result_file.write('\\\\\n')
    

result_file.write('\\begin{tabular}{l | r r r r | r}\n')

result_file.write('Input')
for a in ['KS', 'Range', 'parDss', 'Scan']:
    result_file.write(' & {%s}' % a)
result_file.write(' & serDss \\\\\n')

result_file.write('\hline\n')
for (label, files) in inputs:
    write_label(label)
    for f in files:
        write_result(f)
result_file.write('\\bottomrule\n')
result_file.write('\\end{tabular}\n')
