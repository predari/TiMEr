# adapted from statistics.py by Alexander Noe
# adapted from statistics_resistance.py by Patrick Bisenius
# manipulator: Roland Glantz

import sys
import subprocess
import partition_helper
import os
from graphs import CommunicationGraph

class StatLine(object):
    def __init__(self, duration, cutSize, mcv):
        self._Duration = duration
        self._CutSize = cutSize
        self._MCV = mcv
    def toString(self):
        return ("Duration: {dur}, Cut size: {cutSize}, MCV: {mcv} ").format(dur=self._Duration, cutSize=self._CutSize, mcv=self._MCV)
                        
graphdir = os.path.expanduser('~/c+/KarMa/data/extra_social/')

exePath = os.path.expanduser('~/c+/KarMa/bin/partScenario')

social = [
    CommunicationGraph('/PGPgiantcompo', ''),
    CommunicationGraph('/as-22july06', ''),
    CommunicationGraph('/as-skitter', ''),
    CommunicationGraph('/citationCiteseer', ''),
    CommunicationGraph('/coAuthorsCiteseer', ''),
    CommunicationGraph('/coAuthorsDBLP', ''),
    CommunicationGraph('/coPapersCiteseer', ''),
    CommunicationGraph( '/coPapersDBLP', ''),
    CommunicationGraph( '/email-EuAll', ''),
    CommunicationGraph( '/loc-brightkite_edges', ''),
    CommunicationGraph('/loc-gowalla_edges', ''),
    CommunicationGraph('/p2p-Gnutella04', ''),
    CommunicationGraph('/soc-Slashdot0902', ''),
    CommunicationGraph('/web-Google', ''),
    CommunicationGraph('/wiki-Talk', '')
]

edge_ratings = [    #'weight', 
                    #'expansionstar2',
                    'zero',
                    'resistance'
]
        
statLines = []
stats = []

for currentPGraph in social:
    for k in [2]:
        compare = None
        for edge_rating in edge_ratings:
            #print('Using edge rating ' + edge_rating + " on " + currentPGraph._Graph + ' with k=' + str(k))
            for seed in range(1,31):
                graphFile = graphdir + currentPGraph._Graph + '.graph'
            
                process = subprocess.Popen([exePath, graphFile, '--k=' + str(k), '--preconfiguration=ecosocial', '--edge_rating=' + edge_rating, '--seed=' + str(seed)], stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
                out, err = process.communicate()
                outString = out.decode("UTF-8")

                if outString.startswith( 'Error' ):
                    print (outString)
                    sys.exit(0)

                keyword = outString.find("time spent for partitioning ")
                start = keyword + len("time spent for partitioning ")
                end = outString.find("\n", keyword)
                if (keyword == -1):
                    print(outString)
                    print(err)
                    sys.exit(0)
                duration = float(outString[start:end])

                keyword = outString.find("cut")
                start = keyword + 6
                end = outString.find("\n", keyword)
                if (keyword == -1):
                    print(outString)
                    print(err)
                    sys.exit(0)
                cutSize = int(outString[start:end])
                
                keyword = outString.find("max_comm_vol")
                start = keyword + 13
                end = outString.find("\n", keyword)
                if (keyword == -1):
                    print(outString)
                    print(err)
                    sys.exit(0)
                mcv = int(outString[start:end])

                line = StatLine(duration, cutSize, mcv)
                statLines.append(line)
                #print('cut, mcv, duration are ' + str(cutSize) + ' ' + str(mcv) + ' ' + str(duration))

            stat = partition_helper.analyze(statLines)
            statLines = []
            stats.append(stat)

            
            bigStat = partition_helper.secondAnalyze(stats)
            stats = []
            #print (currentPGraph._Graph + ',' + bigStat.getPrintString())

            if compare is None:
                compare = bigStat
            else:
                bigStat.compare(compare)
