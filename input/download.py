import wget
import os

def download(url, files):
    for f in files:
        print url + f + '.gz'
        if not os.path.isfile(f):
            wget.download(url + f + '.gz')
            os.system('gzip -d -f ' + f + '.gz')
            os.system('../suffixArray/pScan/tools/delete-bytes-255/delete255 ' + f + " > " +  f + '.zero')
            os.system('mv ' + f +'.zero ' + f)
        if not os.path.isfile(f + '.bwt'):
            os.system('../tools/BWT ' + f + ' ' + f + '.bwt')

def writeRep(pat, count, filename):
    if not os.path.isfile(filename):
        f = open(filename, 'w')
        for i in range(0, count):
            f.write(pat)
        f.close()
    if not os.path.isfile(filename + '.bwt'):
        os.system('../tools/BWT ' + filename + ' ' + filename + '.bwt')

def writeRnd(sigma, n):
    filename = 'rnd-%d.bwt' % sigma
    if not os.path.isfile(filename):
        f = open(filename, 'w')
        os.system('../sequenceData/randomSeq -r %d -t int %d %s' %
                (2**sigma, n, filename))

sigmas = [8, 12, 16, 20]
for s in sigmas:
    writeRnd(s, 1024*1024*100)

# Writing 100 MB of some pattern.
writeRep('a', 1024*1024*100, 'aaa')
writeRep(('a'*512*1024)+('b'*512*1024) ,100, 'aabbaabb')
writeRep('ab',1024*1024*50, 'abab')

artificial = ['fib41', 'rs.13', 'tm29']
artificial_url = 'http://pizzachili.dcc.uchile.cl/repcorpus/artificial/'
print 'Downloading artificial repetitive collection'
download(artificial_url, artificial)

real = ['Escherichia_Coli', 'cere', 'coreutils', 'einstein.de.txt',
        'influenza', 'kernel', 'para', 'world_leaders']
real_url = 'http://pizzachili.dcc.uchile.cl/repcorpus/real/'
print 'Downloading real repetitive collection'
download(real_url, real)

pseudo= ['dblp.xml.00001.1', 'dblp.xml.00001.2', 'dblp.xml.0001.1', 'dblp.xml.0001.2',
        'dna.001.1', 'english.001.2', 'proteins.001.1', 'sources.001.2']
pseudo_url = 'http://pizzachili.dcc.uchile.cl/repcorpus/pseudo-real/'
print 'Downloading pseudo-real repetitive collection'
download(pseudo_url, pseudo)

print 'Downloading classic collection'
classic = ['sources', 'pitches', 'proteins', 'dna', 'english.1024MB', 'dblp.xml']
download('http://pizzachili.dcc.uchile.cl/texts/code/', ['sources'])
download('http://pizzachili.dcc.uchile.cl/texts/music/', ['pitches'])
download('http://pizzachili.dcc.uchile.cl/texts/protein/', ['proteins'])
download('http://pizzachili.dcc.uchile.cl/texts/dna/', ['dna'])
download('http://pizzachili.dcc.uchile.cl/texts/nlang/', ['english.1024MB'])
download('http://pizzachili.dcc.uchile.cl/texts/xml/', ['dblp.xml'])


