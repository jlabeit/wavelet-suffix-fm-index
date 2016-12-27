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
classic = ['sources', 'pitches', 'proteins', 'dna', 'english', 'dblp.xml']
download('http://pizzachili.dcc.uchile.cl/texts/code/', ['sources'])
download('http://pizzachili.dcc.uchile.cl/texts/music/', ['pitches'])
download('http://pizzachili.dcc.uchile.cl/texts/protein/', ['proteins'])
download('http://pizzachili.dcc.uchile.cl/texts/dna/', ['dna'])
download('http://pizzachili.dcc.uchile.cl/texts/nlang/', ['english'])
download('http://pizzachili.dcc.uchile.cl/texts/xml/', ['dblp.xml'])
