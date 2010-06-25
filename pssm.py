#################################################################
# @Program: pssm.py                                             #
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
# Copyrighted by Chris Plaisier  12/4/2009                      #
#################################################################

from copy import deepcopy
from math import log

# A class designed to hold a position specific scoring matrix
# and be able to output this in many different formats
#
# Variables:
# pssmName - name
# eValue - the significance of the motif
# nsites - number of genes that have the motif
# genes - the genes that have the motif
# pssmMatrix - a matrix of x by 4 [A, C, G, T]
#
# Functions:
# readPssm(pssmFileName) - Reads in the pssm from the specified file, file name should be as absolute as necessary.
# getMatrix() - returns the pssm matrix as is
# getMemeFormatted() - returns the pssm matrix as a meme 3 formatted string
# colConsensus() - gets the consensus letter for a column of the pssm matrix
# getConsensusMotif() - repeatedly calls colConsensus to get the complete consensus motif
# getName() - return the name of the pssm "bicluster#_motif#<optional _3pUTR>"
# getNSites() - return the number of genes with the motif
# getEValue() - return the significance of the motif for the bicluster
# getGenes() - return the names of the genes that have the motif
#
class pssm:
    # Initialize the pssm
    def __init__(self, pssmFileName=None, biclusterName=None, nsites=None, eValue=None, pssm=None, genes=None):
        if not pssmFileName==None:
            self.name = str(biclusterName)+'_'+(((pssmFileName.split('.'))[0]).split('/'))[-1]
            self.readPssm(pssmFileName)
        else:
            self.name = biclusterName
            self.nsites = nsites
            self.eValue = eValue
            self.matrix = pssm
            self.genes = genes

    # Read in the PSSM matrix
    def readPssm(self, pssmFileName):
        pssmFile = open(pssmFileName,'r')
        firstLine = pssmFile.readline().strip().split(',')
        self.eValue = firstLine[0]
        self.nsites = firstLine[1]
        self.genes = []
        for i in range(int(self.nsites)):
            self.genes.append(pssmFile.readline().strip().split(',')[0])
        self.matrix = []
        self.matrix += [[float(i) for i in line.strip().split(',')] for line in pssmFile.readlines() if line]
        pssmFile.close()

    # Returns the name of the PSSM
    def getName(self):
        return deepcopy(self.name)

    # Sets the name of the PSSM
    def setName(self,name):
        self.name = name

    # Returns the E-value of the PSSM
    def getEValue(self):
        return deepcopy(self.eValue)

    # Returns the number of sites for the PSSM
    def getNSites(self):
        return deepcopy(self.nsites)

    # Returns the number of genes of the PSSM
    def getNumGenes(self):
        return len(self.genes)

    # Returns the genes of the PSSM
    def getGenes(self):
        return deepcopy(self.genes)

    # Returns the matrix
    def getMatrix(self):
        return self.matrix

    # Pads the meme nucleotide frequencies with zeros
    def padMe(self,str1):
        if len(str1)<8:
            for i in range(8-len(str1)):
                str1 += '0'
        return str1

    # Retunrs a log-odds value
    def logOdds(self, p, f):
        p = float(p)
        f = float(f)
        if p==float(0):
            v1 = str(int(round(log(float(1E-300)/f,2),0)))
        else:
            v1 = str(int(round(log(p/f,2),0)))
        if len(v1)<6:
            for i in range(6-len(v1)):
                v1 = ' ' + v1
        return v1

    # Returns a meme 3 formatted string (letter-probability matrix)
    def getMemeFormatted(self,atFreq=0.25,cgFreq=0.25):
        memeFormatted = 'MOTIF '+self.name+'\n'
        memeFormatted += 'BL   MOTIF '+self.name+' width=0 seqs=0\n'
        memeFormatted += 'letter-probability matrix: alength= 4 w= '+str(len(self.matrix))+' nsites= '+self.nsites+' E= '+self.eValue
        for i in self.matrix:
            memeFormatted += '\n '+self.padMe(str(round(float(i[0]),6)))+'  '+self.padMe(str(round(float(i[1]),6)))+'  '+self.padMe(str(round(float(i[2]),6)))+'  '+self.padMe(str(round(float(i[3]),6)))
        return memeFormatted

    # Returns a mast formatted string (log-odds matrix)
    def getMastFormatted(self,atFreq=0.25,cgFreq=0.25):
        mastFormatted = 'log-odds matrix: alength= 4 w= '+str(len(self.matrix))
        for i in self.matrix:
            mastFormatted += '\n '+self.logOdds(i[0],atFreq)+'  '+self.logOdds(i[1],cgFreq)+'  '+self.logOdds(i[2],cgFreq)+'  '+self.logOdds(i[3],atFreq)
        return mastFormatted

   # Returns the consensus word for a motif
    def getConsensusMotif(self, lim1=0.6, lim2=0.8, three=0):
        consensus = ''
        for i in range(len(self.matrix)):
            consensus += self.colConsensus(self.matrix, i, lim1, lim2, three)
        return consensus

    # Returns the consensus letter for a motif column
    def colConsensus(self, pssm, i, lim1, lim2, three):
        two_base_l = ['Y','R','W','S','K','M']
        three_base_l = ['V','H','D','B']
        conLet = 'N'
        if float(pssm[i][0])>=lim1:
            conLet = 'A'
        elif float(pssm[i][1])>=lim1:
            conLet = 'C'
        elif float(pssm[i][2])>=lim1:
            conLet = 'G'
        elif float(pssm[i][3])>=lim1:
            conLet = 'T'
        else:
            two_base_c = [float(pssm[i][1])+float(pssm[i][3]), float(pssm[i][0])+float(pssm[i][2]), float(pssm[i][0])+float(pssm[i][3]), float(pssm[i][1])+float(pssm[i][2]), float(pssm[i][2])+float(pssm[i][3]), float(pssm[i][0])+float(pssm[i][1])]
            three_base_c = [float(pssm[i][0])+float(pssm[i][1])+float(pssm[i][2]), float(pssm[i][0])+float(pssm[i][1])+float(pssm[i][3]), float(pssm[i][0])+float(pssm[i][2])+float(pssm[i][3]), float(pssm[i][1])+float(pssm[i][2])+float(pssm[i][3])]
            pMax = 0
            for k in range(0,6):
                if two_base_c[k] > pMax:
                    pMax = two_base_c[k]
                    conLet = two_base_l[k]
            if not pMax>lim2 and three==1:
                for k in range(0,4):
                    if three_base_c[k] > pMax:
                        pMax = three_base_c[k]
                        conLet = three_base_l[k]
            if not pMax>lim2:
                conLet = 'N'
        return conLet

    
