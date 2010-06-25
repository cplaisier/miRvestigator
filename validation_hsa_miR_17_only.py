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

import sys, re, os, math, shutil
from pssm import pssm
from multiprocessing import Pool, cpu_count, Manager
from subprocess import *
from copy import deepcopy
from random import sample
import cPickle
from numpy import array, float64, log10
from weblogolib import *
import numpy, corebio

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

# Run meme and get the output into PSSMs
def meme(seqFile=None, bgFile=None, nMotifs=1, minMotifWidth=6, maxMotifWidth=12, revComp=False, seed=None):
    if not os.path.exists('tmp/meme'):
        os.makedirs('tmp/meme')
    # Arguments for tomtom
    memeArgs = str(seqFile)+' -bfile '+str(bgFile)+' -nostatus -text -time 600 -dna -maxsize 9999999 -evt 1e9 -mod zoops -minw ' + str(minMotifWidth) + ' -maxw ' + str(maxMotifWidth) + ' -nmotifs ' + str(nMotifs)
    #if revComp==True:
    #    memeArgs += ' -revcomp'
    if not seed==None:
        memeArgs += ' -cons ' + str(seed)
    print memeArgs
    #errOut = open('tmp/meme/stderr.out','w')
    memeProc = Popen("meme " + memeArgs, shell=True,stdout=PIPE) #,stderr=errOut)
    output = memeProc.communicate()[0].split('\n')
    #errOut.close()
    
    PSSMs = []
    # Now iterate through output and save data
    for i in range(len(output)):
        splitUp1 = output[i].strip().split(' ')
        if splitUp1[0]=='Motif' and splitUp1[2]=='position-specific' and splitUp1[3]=='probability':
            i += 2 # Skip the separator line, go to the summary line
            splitUp = output[i].strip().split(' ')
            width = int(splitUp[5])
            sites = splitUp[7]
            eValue = splitUp[9]
            matrix = []
            for j in range(width):
                i += 1
                matrix += [[float(let) for let in output[i].strip().split(' ') if let]]
            PSSMs.append(pssm(biclusterName=str(splitUp1[1]),nsites=sites,eValue=eValue,pssm=matrix,genes=[]))
    return PSSMs

# Function to run the meme function
def runMeme(i):
    #print clusterFileNames[i]
    meme(i,seqFile=clusterFileNames[i],bgFile=allVars['bgFile'],nMotifs=allVars['nMotifs'],minMotifWidth=allVars['minMotifWidth'], maxMotifWidth=allVars['maxMotifWidth'], revComp=allVars['revComp'], seed=seeds[i])

# Complement
def complement(seq):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    complseq = [complement[base] for base in seq]
    return complseq

# Reverse complement
def reverseComplement(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(complement(seq))

# Get a list of all possible nMers
def permute(depth, letters=['A','C','G','T'], seqs=[''],curdepth=0):
    newseqs = []
    for seq in seqs:
        for letter in letters:
            newseqs.append(seq + letter)
    if depth > curdepth:
        return(permute(depth,letters,newseqs,curdepth + 1))
    else:
        return(seqs)

# Create the background
def getBackground(allSeqs):
    outFile = open('tmp/meme/bgFile.meme','w')
    for w in range(1,3):
        allnmers = permute(w)
        nmersT = topNmers(w, allSeqs, True, True)
        nmersD = {}
        total = 0
        for nmer in allnmers:
            nmersD[nmer] = 1 #Pseudo count
            total = total + 1
        for nmer,count in nmersT[:]:
            try: 
                rc = reverseComplement(nmer)
                nmersD[nmer] = nmersD[nmer] + count
                nmersD[rc]   = nmersD[rc]   + count
                total = total + 2*count
            except KeyError:
                pass
        tmp1 = nmersD.keys()
        tmp1.sort()
        outFile.write('# '+str(w-1)+'th order Markov background model\n')
        for nmer in tmp1:
            outFile.write(str(nmer)+' '+str(float(nmersD[nmer])/total)+'\n')
    outFile.close()

# Assemble list of all nmers (kmers) with width n from supplied sequences.
# Option with_counts returns list of (kmer, count) tuples instead.
# Purge N's ignores kmers containing N's.
def topNmers(order, seqs, withCounts=True, purgeNs=True):
    nMers = {}
    for seq in seqs:
        for i in range(len(seq)-order+1):
            nMer = seq[i:i+order]
            if purgeNs:
                if nMer.find('N') >= 0: continue
            nMerRC = reverseComplement(nMer)
            tmp1 = [nMer, nMerRC]
            tmp1.sort()
            nMerKey = tmp1[0]        # _t used until here to get alphabetically first seq
            if nMers.has_key(nMerKey):
                nMers[nMerKey] = nMers[nMerKey] + 1
            else:
                nMers[nMerKey] = 1
    sorted = nMers.keys()
    sorted.sort(lambda x,y,D=nMers:cmp(D[y],D[x]) or cmp(x,y))
    if withCounts:
        return(zip(sorted,map(lambda x,N=nMers:N[x], sorted)))
    else:
        return(sorted)

# 1. Read in hsa-miR-17 targets: per target ['entrez IDs', 'affy-probes separated by spaces']
inFile = open('hsa_miR_17_targets_only.csv','r')
targets = {}
inLines = inFile.readlines()
for i in inLines:
    splitUp = i.strip().split(',')
    targets[splitUp[0]] = splitUp[1].split(' ')
inFile.close()

# 2. Read in sequences
# Could use conserved sequences 'p3utrSeqs_set3pUTR_Final_c.csv'
seqFile = open('p3utrSeqs_set3pUTR_Final.csv','r')
#seqFile = open('p3utrSeqs_set3pUTR_Final_c.csv','r')
seqLines = seqFile.readlines()
ids = [i.strip().split(',')[0].upper() for i in seqLines]
sequences = [i.strip().split(',')[1] for i in seqLines]
seqs = dict(zip(ids,sequences))
seqFile.close()

# 3. Get sequences for each target
miR17Seqs = {}
for target in targets:
    foundOne = 0
    for probe in targets[target]:
        if probe in seqs:
            miR17Seqs[target] = seqs[probe]
            foundOne = 1
            break
    if foundOne==0:
        print 'Did not find seq for',target

# 4. Make a FASTA file
if not os.path.exists('tmp/meme/fasta'):
    os.makedirs('tmp/meme/fasta')
fastaFile = open('tmp/meme/fasta/hsa_miR_17.fasta','w')
for seq in miR17Seqs:
    fastaFile.write('>'+str(seq)+'\n'+str(miR17Seqs[seq])+'\n')
fastaFile.close()

# 5. Make a background file
# Read in hsa-miR-17 targets: per target ['entrez IDs', 'affy-probes separated by spaces']
inFile = open('bkgdEntrezIds.csv','r')
bgs = {}
inLines = inFile.readlines()
for i in inLines:
    splitUp = i.strip().split(',')
    bgs[splitUp[0]] = splitUp[1].split(' ')
inFile.close()
# Get seqs for each bg
bgSeqs = {}
for bg in bgs:
    foundOne = 0
    for probe in bgs[bg]:
        if probe in seqs:
            bgSeqs[bg] = seqs[probe]
            foundOne = 1
            break
    if foundOne==0:
        print 'Did not find seq for',bg
# Make background distribution for meme runs
bgFile = 'tmp/meme/bgFile.meme'
if not os.path.exists(bgFile):
    print 'Making bacground distribution...'
    getBackground(bgSeqs.values())
    print 'Done.\n'
else:
    print 'Using previously computed background distribution file.'


# 5. Run MEME
memePSSMs1 = meme(seqFile='tmp/meme/fasta/hsa_miR_17.fasta', bgFile=bgFile, nMotifs=1, minMotifWidth=4, maxMotifWidth=9, revComp=False, seed=None)
print memePSSMs1[0].getConsensusMotif(),memePSSMs1[0].getEValue()

# 6. Look at PSSM (plot it preferably)
plotPssm(memePSSMs1[0],'pssm0only.png')

# 5. Run MEME
memePSSMs2 = meme(seqFile='tmp/meme/fasta/hsa_miR_17.fasta', bgFile=bgFile, nMotifs=1, minMotifWidth=4, maxMotifWidth=9, revComp=False, seed='GCACTTTG')
print memePSSMs2[0].getConsensusMotif(),memePSSMs2[0].getEValue()

# 6. Look at PSSM (plot it preferably)
plotPssm(memePSSMs2[0],'pssmSeedonly.png')

# 7. Compare to miRDB using my program
from miRNA2motifHMM import miRNA2motifHMM
m2m = miRNA2motifHMM(memePSSMs1,seqs.values(),seedModel=[6,7,8],minor=True,p5=True,p3=True)
print m2m.getTopHit('1')

