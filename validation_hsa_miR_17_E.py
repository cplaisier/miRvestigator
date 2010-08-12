#################################################################
# @Program: validation_hsa_miR_17.py                            #
# @Version: 1                                                   #
# @Author: Chris Plaisier                                       #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 1441 North 34th Street                                        #
# Seattle, Washington  98103-8904                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyrighted by Chris Plaisier  6/18/2010                      #
#################################################################

# Libraries needed to run
import sys, re, os, math, shutil
from subprocess import *
from copy import deepcopy
from random import sample
import cPickle

# Libraries for plotting
import numpy, corebio                     # http://numpy.scipy.org and http://code.google.com/p/corebio/
from numpy import array, float64, log10   # http://numpy.scipy.org
from weblogolib import *                  # http://code.google.com/p/weblogo/

# Custom library

# Plot a PSSM using weblogo
def plotPssm(pssm, fileName):
    dist = numpy.array( pssm.getMatrix(), numpy.float64 ) 
    data = LogoData.from_counts(corebio.seq.unambiguous_dna_alphabet, dist*100)
    options = LogoOptions()
    options.color_scheme = colorscheme.nucleotide
    format = LogoFormat(data, options)
    fout = open(fileName, 'w')
    png_formatter(data, format, fout)
    fout.close()

# Run weeder and parse its output
# First weederTFBS -W 6 -e 1, then weederTFBS -W 8 -e 2, and finally adviser
# The weeder program can be found at:  http://159.149.109.9/modtools/
# I modified the C code and recompiled to make Weeder look for the FreqFiles
# folder in /local/FreqFiles. Then I made symbolic links in my PATH so that
# weeder could be run from the command line as weederlauncher. You will also
# have to add weederTFBS.out and adviser.out to the PATH in order to run.
def weeder(seqFile=None, percTargets=50, revComp=False):
    if not os.path.exists('tmp/weeder'):
        os.makedirs('tmp/weeder')
    
    # First run weederTFBS for 6bp motifs
    weederArgs = ' '+str(seqFile)+' HS small T50'
    if revComp==True:
        weederArgs += ' -S'
    errOut = open('tmp/weeder/stderr.out','w')
    weederProc = Popen("weederlauncher " + weederArgs, shell=True,stdout=PIPE,stderr=errOut)
    output = weederProc.communicate()
    
    """# First run weederTFBS for 6bp motifs
    weederArgs = '-f '+str(seqFile)+' -W 6 -e 1 -O HS -R '+str(percTargets)
    if revComp==True:
        weederArgs += ' -S'
    errOut = open('tmp/weeder/stderr.out','w')
    weederProc = Popen("weeder " + weederArgs, shell=True,stdout=PIPE,stderr=errOut)
    output = weederProc.communicate()
    
    # Second run weederTFBS for 8bp motifs
    weederArgs = '-f '+str(seqFile)+' -W 8 -e 2 -O HS -R '+str(percTargets)
    if revComp==True:
        weederArgs += ' -S'
    weederProc = Popen("weeder " + weederArgs, shell=True,stdout=PIPE,stderr=errOut)
    output = weederProc.communicate()
    
    # Finally run adviser
    weederArgs = str(seqFile)
    weederProc = Popen("adviser " + weederArgs, shell=True,stdout=PIPE,stderr=errOut)
    output = weederProc.communicate()
    errOut.close()
    """
    # Now parse output from weeder
    PSSMs = []
    output = open(str(seqFile)+'.wee','r')
    outLines = [line for line in output.readlines() if line.strip()]
    hitBp = {}
    # Get top hit of 6bp look for "1)"
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('1) ') == -1:
            break
    hitBp[6] = outLine.strip().split(' ')[1:]

    # Scroll to where the 8bp reads wll be
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Searching for motifs of length 8') == -1:
            break

    # Get top hit of 8bp look for "1)"
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('1) ') == -1:
            break
    hitBp[8] = outLine.strip().split(' ')[1:]

    # Scroll to where the 8bp reads wll be
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Your sequences:') == -1:
            break
    
    # Get into the highest ranking motifs
    seqDict = {}
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('**** MY ADVICE ****') == -1:
            break
        splitUp = outLine.strip().split(' ')
        seqDict[splitUp[1]] = splitUp[3].lstrip('>')

    # Get into the highest ranking motifs
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Interesting motifs (highest-ranking)') == -1:
            break
    while 1:
        name = outLines.pop(0).strip() # Get match
        if not name.find('(not highest-ranking)') == -1:
            break
        # Get redundant motifs
        outLines.pop(0)
        redMotifs = [i for i in outLines.pop(0).strip().split(' ') if not i=='-']
        outLines.pop(0)
        outLines.pop(0)
        line = outLines.pop(0)
        instances = []
        while line.find('Frequency Matrix') == -1:
            splitUp = [i for i in line.strip().split(' ') if i]
            instances.append({'gene':seqDict[splitUp[0]], 'strand':splitUp[1], 'site':splitUp[2], 'start':splitUp[3], 'match':splitUp[4].lstrip('(').rstrip(')') })
            line = outLines.pop(0)
        # Read in Frequency Matrix
        outLines.pop(0)
        outLines.pop(0)
        matrix = []
        col = outLines.pop(0)
        while col.find('======') == -1:
            nums = [i for i in col.strip().split('\t')[1].split(' ') if i]
            colSum = 0
            for i in nums:
                colSum += int(i.strip())
            matrix += [[ float(nums[0])/float(colSum), float(nums[1])/float(colSum), float(nums[2])/float(colSum), float(nums[3])/float(colSum)]]
            col = outLines.pop(0)
        PSSMs += [pssm(biclusterName=name,nsites=instances,eValue=hitBp[len(matrix)][1],pssm=matrix,genes=redMotifs)]
    return PSSMs

# 1. Read in hsa-miR-17 targets: per target ['entrez IDs', 'affy-probes separated by spaces']
inFile = open('hsa_miR_17_E.csv','r')
targets = [i.strip() for i in inFile.readlines()]
inFile.close()

# 2. Read in sequences
seqFile = open('p3utrSeqs_Homo_sapiens.csv','r')
seqLines = seqFile.readlines()
ids = [i.strip().split(',')[0].upper() for i in seqLines]
sequences = [i.strip().split(',')[1] for i in seqLines]
seqs = dict(zip(ids,sequences))
seqFile.close()

# 3. Get sequences for each target
miR17Seqs = {}
for target in targets:
    if target in seqs:
        miR17Seqs[target] = seqs[target]
    else:
        print 'Did not find seq for',target

# 4. Make a FASTA file
if not os.path.exists('tmp/weeder/fasta'):
    os.makedirs('tmp/weeder/fasta')
fastaFile = open('tmp/weeder/fasta/hsa_miR_17.fasta','w')
for seq in miR17Seqs:
    fastaFile.write('>'+str(seq)+'\n'+str(miR17Seqs[seq])+'\n')
fastaFile.close()

# 5. Run Weeder
weederPSSMs1 = weeder(seqFile='tmp/weeder/fasta/hsa_miR_17.fasta', percTargets=50, revComp=False)
if len(weederPSSMs1)>0:
    print weederPSSMs1[0].getConsensusMotif(),weederPSSMs1[0].getEValue()
if len(weederPSSMs1)>1:
    print weederPSSMs1[1].getConsensusMotif(),weederPSSMs1[1].getEValue()

# 6. Look at PSSM (plot it preferably) but you can comment this out if you don't want to install
#    weblogolib
for i in range(len(weederPSSMs1)):
    plotPssm(weederPSSMs1[i],'pssm'+str(i)+'only.png')

# 7. Compare to miRDB using my program
from miRvestigator import miRvestigator
mV = miRvestigator(weederPSSMs1,seqs.values(),seedModel=[6,7,8],minor=True,p5=True,p3=True,wobble=True,wobbleCut=0.25)
print mV.getTopHit(weederPSSMs1[0].getName())

