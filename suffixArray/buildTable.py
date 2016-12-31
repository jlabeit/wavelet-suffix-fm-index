import csv

# Read result data.

def parseRow(row):
    return [row[0].strip(), row[1].strip(), int(row[2]), int(row[3]), float(row[4])]

data = []
with open('results', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in spamreader:
        data.append(parseRow(row))
print data

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
palgos = ['parallelDivsufsort', 'parallelKS', 'parallelRange', 'pScan']


result_file = open('resultSuffix.tex', 'w')

def getTime(a, b, c):
    for row in data:
        if row[0] == a and row[1] == b and row[2] == c:
            return row[4]
    return 1.0;

def write_label(l):
    result_file.write('\\multicolumn{2}{l}{\\textit{' + label + '}} \\\\\n')

def getTiming(f, algo):
    if (algo in salgos):
        return '%.2f '% getTime(algo, f, 1)
    else:
        t1 = getTime(algo, f, 1)
        t64 = getTime(algo, f, 64)
        return '%.2f & %.2f & %.2f &' % (t1, t64, t1 / t64) 

def write_result(f):
    result_file.write(f + ' & ')
    for algo in palgos:
        result_file.write(getTiming(f, algo))
    for algo in salgos:
        result_file.write(getTiming(f, algo))
    result_file.write('\\\\\n')
    

result_file.write('\\begin{tabular}{l | r r r | r r r | r r r | r r r | r}\n')
result_file.write('& \\multicolumn{3}{c}{parDss}  & \\multicolumn{3}{c}{KS} & \\multicolumn{3}{c}{Range} & \multicolumn{3}{c}{Scan} & serDss \\\\\n')
result_file.write('Input & & & & & & & & & & & & & \\\\\n')
result_file.write('\hline\n')
for (label, files) in inputs:
    write_label(label)
    for f in files:
        write_result(f)
result_file.write('\\bottomrule\n')
result_file.write('\\end{tabular}\n')
