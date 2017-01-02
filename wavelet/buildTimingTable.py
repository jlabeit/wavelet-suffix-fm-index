import csv

# Read result data.

def parseRow(row):
    return [row[0].strip(), row[1].strip(), int(row[2]), int(row[3]), float(row[4])]

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
integer = ('integer', ['rnd-8', 'rnd-12', 'rnd-16', 'rnd-20'])

inputs = [classic, real, artificial, pseudo, special, integer]
#threads = [1, 2, 4, 8, 12, 16, 24, 32, 40, 48, 56, 64]
threads = [1, 8, 16, 32, 64]

salgos = ['serialWT', 'sdslWT']
palgos = ['levelWT', 'recWT', 'ddWT']

result_file = open('resultsWavelet.tex', 'w')

def getTime(a, b, c):
    for row in data:
        if row[0] == a and row[1] == b and row[2] == c:
            return row[4]
    return 1.0;

def write_label(l):
    result_file.write('\\multicolumn{2}{l}{\\textit{' + label + '}} \\\\\n')

def getTiming(f, algo, fastest):
    if (algo in salgos):
        if algo == 'sdslWT':
            return '%.2f '% getTime(algo, f, 1)
        else:
            return '%.2f & '% getTime(algo, f, 1)
    else:
        t1 = getTime(algo, f, 1)
        t64 = getTime(algo, f, 64)
        if (fastest):
            return '%.2f & \\textbf{%.2f} & %.2f &' % (t1, t64, t1 / t64) 
        return '%.2f & %.2f & %.2f &' % (t1, t64, t1 / t64) 

def getFastest(f):
    mini = 100000000.0;
    result = ''
    for algo in palgos:
        cur = getTime(algo, f, 64)
        if mini > cur:
            mini = cur
            result = algo
    return result

def write_result(f):
    result_file.write(f.replace('_', '\\_') + ' & ')
    fastest = getFastest(f)
    for algo in palgos:
        result_file.write(getTiming(f, algo, fastest == algo))
    for algo in salgos:
        result_file.write(getTiming(f, algo, fastest == algo))
    result_file.write('\\\\\n')
    

result_file.write('\\begin{tabular}{l | r r r | r r r | r r r | r r}\n')

for a in ['levelWT', 'recWT', 'ddWT']:
    result_file.write('& \\multicolumn{3}{c|}{%s} ' % a)
result_file.write(' & serWT & sdslWT \\\\\n')

result_file.write('Input & ')
for a in ['levelWT', 'recursiveWT', 'ddWT']:
    result_file.write('$T_1$ & $T_{64}$ & $\\frac{T_1}{T_{64}}$ & ')
result_file.write('$T_1$ & $T_1$ \\\\\n')

result_file.write('\hline\n')
for (label, files) in inputs:
    write_label(label)
    for f in files:
        write_result(f)
result_file.write('\\bottomrule\n')
result_file.write('\\end{tabular}\n')
