__author__ = 'alex'

from math import log10

ANZAHL_STELLEN = 6


def toStr(x):
    n = int(log10(x)) + 1 
    zahl = str(round(x, ANZAHL_STELLEN - n))

    #mit 0en auffuellen
    nmb_digits = len(zahl) - 1
    if nmb_digits < ANZAHL_STELLEN:
        zahl += "0" * (ANZAHL_STELLEN - nmb_digits)
        
    return zahl

def mean(list):
    sum = 0
    for x in list:
        sum += x

    return sum / len(list)

def gmean(list):
    prod = 1
    for x in list:
        prod *= x

    return prod ** (1.0 / len(list))

class Stat(object):
    def __init__(self, ints, floats):
        self._Ints = ints
        self._Floats = floats

    def printString(self):
        line = 'Qdur: ' + toStr(self._Floats[3] / self._Floats[0]) + ' & ' + toStr(self._Floats[4] / self._Floats[1]) + ' & ' + toStr(self._Floats[5] / self._Floats[2]) + ' QavgD: ' + toStr(self._Floats[9] / self._Floats[6]) + ' & ' + toStr(self._Floats[10] / self._Floats[7]) + ' & ' + toStr(self._Floats[11] / self._Floats[8]) + ' QmaxC: ' + toStr(float(self._Ints[3]) / float(self._Ints[0])) + ' & ' + toStr(float(self._Ints[4]) / float(self._Ints[1])) + ' & ' + toStr(float(self._Ints[5]) / float(self._Ints[2])) + ' QmaxD: ' + toStr(float(self._Ints[9]) / float(self._Ints[6])) + ' & ' + toStr(float(self._Ints[10]) / float(self._Ints[7])) + ' & ' + toStr(float(self._Ints[11]) / float(self._Ints[8])) + ' QCut: ' + toStr(float(self._Ints[15]) / float(self._Ints[12])) + ' & ' + toStr(float(self._Ints[16]) / float(self._Ints[13])) + ' & ' + toStr(float(self._Ints[17]) / float(self._Ints[14])) + ' QTcv: ' + toStr(float(self._Ints[21]) / float(self._Ints[18])) + ' & ' + toStr(float(self._Ints[22]) / float(self._Ints[19])) + ' & ' + toStr(float(self._Ints[23]) / float(self._Ints[20]))

        print line


def analyze(list):
    dur1 = []
    dur2 = []
    aDil1 = []
    aDil2 = []
    mCon1 = []
    mCon2 = []
    mDil1 = []
    mDil2 = []
    cut1  = []
    cut2  = []
    tcv1  = []
    tcv2  = []

    for line in list:
        dur1.append(line._Duration1)
        dur2.append(line._Duration2)
        aDil1.append(line._AvgDil1)
        aDil2.append(line._AvgDil2)
        mCon1.append(line._MaxCon1)
        mCon2.append(line._MaxCon2)
        mDil1.append(line._MaxDil1)
        mDil2.append(line._MaxDil2)
        cut1.append(line._CUT1)
        cut2.append(line._CUT2)
        tcv1.append(line._TCV1)
        tcv2.append(line._TCV2)

    floats = []
    ints = []
	
    floats.append(min(dur1))
    floats.append(mean(dur1))
    floats.append(max(dur1))
    floats.append(min(dur2))
    floats.append(mean(dur2))
    floats.append(max(dur2))
    floats.append(min(aDil1))
    floats.append(mean(aDil1))
    floats.append(max(aDil1))
    floats.append(min(aDil2))
    floats.append(mean(aDil2))
    floats.append(max(aDil2))

    ints.append(min(mCon1))
    ints.append(mean(mCon1))
    ints.append(max(mCon1))
    ints.append(min(mCon2))
    ints.append(mean(mCon2))
    ints.append(max(mCon2))
    ints.append(min(mDil1))
    ints.append(mean(mDil1))
    ints.append(max(mDil1))
    ints.append(min(mDil2))
    ints.append(mean(mDil2))
    ints.append(max(mDil2))
    ints.append(min(cut1))
    ints.append(mean(cut1))
    ints.append(max(cut1))
    ints.append(min(cut2))
    ints.append(mean(cut2))
    ints.append(max(cut2))
    ints.append(min(tcv1))
    ints.append(mean(tcv1))
    ints.append(max(tcv1))
    ints.append(min(tcv2))
    ints.append(mean(tcv2))
    ints.append(max(tcv2))

    newS = Stat(ints, floats)

    return newS

def secondAnalyze(list):
    ints = []
    floats = []

    for x in range(0,12):
        floats.append(geomean([o._Floats[x] for o in list]))
    for x in range(0,24):
        ints.append(geomean([o._Ints[x] for o in list]))

    newS = Stat(ints, floats)
    return newS	

def geomean(list):
    product = 1.0
    for n in list:
        product *= (n ** (1.0/len(list)))
    return (product)
