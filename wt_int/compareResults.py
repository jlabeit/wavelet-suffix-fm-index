from subprocess import call

testProgram = "./WT"
answerProgram = "./WTans"

inputDir = "../sequenceData/data/"
inputFiles = [ "trigramString_10000000", "chr22.dna", "etext99", "randomSeq_", "wikisamp.xml" ]


for f in inputFiles:
    # Run both programs
    call([testProgram, "-o",  "temp.out" , inputDir + f])
    call([answerProgram, "-o",  "temp.ans" , inputDir + f])
    call(["cmp", "temp.out", "temp.ans"])

