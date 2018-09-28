__author__ = 'alex'
#manipulators Patrick Bisenius, Roland Glantz

from math import log10

ANZAHL_STELLEN = 5

def toStr(x):
    return "{:.5f}".format(x)

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

    def getPrintString(self):
        line = 'cut::: ' + toStr(self._Ints[0]) + ' & ' + toStr(self._Ints[1]) + ' & ' + toStr(self._Ints[2]) + ' mcv::: ' + toStr(self._Ints[3]) + ' & ' + toStr(self._Ints[4]) + ' & ' + toStr(self._Ints[5]) + ' duration::: '  + toStr(self._Floats[0]) + ' & ' + toStr(self._Floats[1]) + ' & ' + toStr(self._Floats[2]) + '"'
        return line

    def compare(self, comp):
        print 'Cut: ' + toStr(float(self._Ints[0]) / float(comp._Ints[0])) + ' & ' + toStr(float(self._Ints[1]) / float(comp._Ints[1])) + ' & ' + toStr(float(self._Ints[2]) / float(comp._Ints[2])) + ' MCV: ' + toStr(float(self._Ints[3]) / float(comp._Ints[3])) + ' & ' + toStr(float(self._Ints[4]) / float(comp._Ints[4])) + ' & ' + toStr(float(self._Ints[5]) / float(comp._Ints[5])) + ' Duration: ' + toStr(self._Floats[0] / comp._Floats[0]) + ' & ' + toStr(self._Floats[1] / comp._Floats[1]) + ' & ' + toStr(self._Floats[2] / comp._Floats[2])

def analyze(list):
    cuts = []
    mcv = []
    dur = []
    for line in list:
        cuts.append(line._CutSize)
        mcv.append(line._MCV)
        dur.append(line._Duration)
    ints = []
    floats = []
    
    floats.append(min(dur))
    floats.append(mean(dur))
    floats.append(max(dur))

    ints.append(min(cuts))
    ints.append(mean(cuts))
    ints.append(max(cuts))
    ints.append(min(mcv))
    ints.append(mean(mcv))
    ints.append(max(mcv))

    newS = Stat(ints, floats)

    return newS

def secondAnalyze(list):
    ints = []
    floats = []

    for x in range(0,3):
        floats.append(geomean([o._Floats[x] for o in list]))

    for x in range(0,6):
        ints.append(geomean([o._Ints[x] for o in list]))

    newS = Stat(ints, floats)
    return newS	

def geomean(list):
    product = 1.0
    for n in list:
        product *= (n ** (1.0/len(list)))
    return (product)
