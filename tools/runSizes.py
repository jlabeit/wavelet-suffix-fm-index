import os

artificial = ('artificial repetitive', ['fib41', 'rs.13', 'tm29'])
real = ('real repetitive', ['Escherichia_Coli', 'cere', 'coreutils', 'einstein.de.txt',
        'influenza', 'kernel', 'para', 'world_leaders'])
pseudo= ('pseudo repetitive', ['dblp.xml.00001.1', 'dblp.xml.00001.2', 'dblp.xml.0001.1',
    'dblp.xml.0001.2', 'dna.001.1', 'english.001.2', 'proteins.001.1', 'sources.001.2'])
classic = ('non-repetitive', ['sources', 'pitches', 'proteins', 'dna', 'english.1024MB', 'dblp.xml'])
special = ('special cases', ['aaa', 'abab', 'aabbaabb'])

inputs = [classic, real, artificial, pseudo, special]
result_file = open('sizes', 'w')

def getSize(f):
    return os.path.getsize('../input/' + f)
result = dict()
for (label, files) in inputs:
    for f in files:
        result[f] = getSize(f)
result_file.write(str(result))
