import numpy as np
import matplotlib.pyplot as plt
import csv

# Read result data.
def parseRow(row):
    return [row[0].strip(), row[1].strip(), int(row[2]), int(row[3]), float(row[4])]

data = []
with open('results', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in spamreader:
        data.append(parseRow(row))

salgos = ['serialWT', 'sdslWT']
palgos = [('b', 'o', 'levelWT'),
          ('g', '', 'recWT'),
          ('r', 'x', 'ddWT')]

def getTime(algo, f, p):
    for row in data:
        if row[0] == algo and row[1] == f and row[2] == p:
            return row[4]
    return 1.0;



def getSpeedup(algo, f):
    return [getTime('serialWT', f, 1) / getTime(algo, f, p) for p in threads]
lines = dict()

def makePlot(plot, filename):
    for (c, m, algo) in palgos:
        speedup = getSpeedup(algo, filename)
        lines[algo], = plot.plot(threads, speedup, c, marker=m, label=algo)
    plot.set_xticks([1, 4, 8, 12, 16, 24, 32, 40, 48, 56, 64])
    plot.set_title(filename.replace('_', '\\_'), fontsize=16)

# Example data
threads = [1, 2, 4, 8, 12, 16, 24, 32, 40, 48, 56, 64]

f, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharex='col', sharey='row')

makePlot(ax1, 'english.1024MB')
makePlot(ax2, 'dna')
makePlot(ax3, 'coreutils')
makePlot(ax4, 'world_leaders')
makePlot(ax5, 'fib41')
makePlot(ax6, 'rs.13')
makePlot(ax7, 'dblp.xml.00001.1')
makePlot(ax8, 'dblp.xml.0001.2')
ax1.set_ylabel(r'absolute speedup', fontsize=16)
ax3.set_ylabel(r'absolute speedup', fontsize=16)
ax5.set_ylabel(r'absolute speedup', fontsize=16)
ax7.set_ylabel(r'absolute speedup', fontsize=16)
ax7.set_xlabel(r'number of threads', fontsize=16)
ax8.set_xlabel(r'number of threads', fontsize=16)

f.legend([lines['levelWT'], lines['recWT'], lines['ddWT']], 
        ['levelWT', 'recWT', 'ddWT'], 'lower center', ncol=3)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.savefig('waveletPlot.eps', format='eps')
plt.show()
