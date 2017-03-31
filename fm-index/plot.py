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

salgos = ['ser']
palgos = [('black', '|' , 'par')]

def getTime(algo, f, p):
    for row in data:
        if row[0] == algo and row[1] == f and row[2] == p:
            return row[4]
    return 1.0;



def getSpeedup(algo, f):
    return [getTime('ser', f, 1) / getTime(algo, f, p) for p in threads]
lines = dict()

def makePlot(plot, filename):
    for (c, m, algo) in palgos:
        speedup = getSpeedup(algo, filename)
        lines[algo], = plot.plot(threads, speedup, c, marker=m, label=algo)
    plot.set_xticks([1, 4, 8, 12, 16, 24, 32, 40, 48, 56, 64])
    plot.set_title(filename.replace('_', '\\_'), fontsize=16)
    plot.set_ylim([0,5])

# Example data
threads = [1, 2, 4, 8, 12, 16, 24, 32, 40, 48, 56, 64]

f, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8), (ax9, ax10)) = plt.subplots(5, 2, sharex='col')

makePlot(ax1, 'dna')
makePlot(ax2, 'sources')
makePlot(ax3, 'cere')
makePlot(ax4, 'world_leaders')
makePlot(ax5, 'tm29')
makePlot(ax6, 'fib41')
makePlot(ax7, 'proteins.001.1')
makePlot(ax8, 'dna.001.1')
makePlot(ax9, 'proteins.001.1')
makePlot(ax10, 'dna.001.1')
ax1.set_ylabel(r'absolute speedup', fontsize=16)
ax3.set_ylabel(r'absolute speedup', fontsize=16)
ax5.set_ylabel(r'absolute speedup', fontsize=16)
ax7.set_ylabel(r'absolute speedup', fontsize=16)
ax9.set_ylabel(r'absolute speedup', fontsize=16)
ax9.set_xlabel(r'number of threads', fontsize=16)
ax10.set_xlabel(r'number of threads', fontsize=16)


f.legend([lines['par']], ['par'], 'lower center', ncol=4)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.savefig('fmindexPlot.eps', format='eps')
plt.show()
