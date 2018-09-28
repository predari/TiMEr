#_author__ = 'alex'
#manipulator Roland Glantz

import sys
import subprocess
import map_helper
import os
from graphs import ProcessorGraph

class StatLine(object):
    def __init__(self, duration1, duration2, avgDil1, avgDil2, maxCon1, maxCon2, maxDil1, maxDil2, cut1, cut2, tcv1, tcv2):
        self._Duration1 = duration1
        self._Duration2 = duration2
        self._AvgDil1 = avgDil1
        self._AvgDil2 = avgDil2
        self._MaxCon1 = maxCon1
        self._MaxCon2 = maxCon2
        self._MaxDil1 = maxDil1
        self._MaxDil2 = maxDil2
        self._CUT1 = cut1
        self._CUT2 = cut2
        self._TCV1 = tcv1
        self._TCV2 = tcv2
        
    def toString(self):
        return ("Mapping duration1: {duration1}, Mapping duration2: {duration2}, Average dilation1: {avgDil1}, Average dilation2: {avgDil2}, Maximum congestion1: {maxCon1}, Maximum congestion2: {maxCon2}, Maximum dilation1: {maxDil1}, Maximum dilation2: {maxDil2}, Cut1: {cut1}, Cut2: {cut2}, Total CV1: {tcv1}, Total CV2: {tcv2}.format(duration1=self._Duration1, duration2=self._Duration2, avgDil1=self._AvgDil1, avgDil2=self._AvgDil2, maxCon1=self._MaxCon1, maxCon2=self._MaxCon2, maxDil1=self._MaxDil1, maxDil2=self._MaxDil2, cut1=self._CUT1, cut2=self._CUT2, tcv1=self._TCV1, tcv2=self._TCV2)")

#graphDir = os.path.expanduser('~/c+/KarMa/data/graphs/extra_social/')
#partDir = os.path.expanduser('~/c+/KarMa/data/partitions/extra_social/')
#mapDir = os.path.expanduser('~/c+/KarMa/data/mappings/extra_social/')

graphDir = os.path.expanduser('~/c+/KarMa/data/graphs/walshaw/')
partDir = os.path.expanduser('~/c+/KarMa/data/partitions/walshaw/')
mapDir = os.path.expanduser('~/c+/KarMa/data/mappings/walshaw/')

exePath = os.path.expanduser('~/c+/KarMa/KarMa_sequential/bin/mapScenario')

Pset = [ProcessorGraph('256',  'grid',  '16x16', '1x1x1x1'  ),
#        ProcessorGraph('1024', 'grid',  '32x32', '1x1x1x1x1'),
        ProcessorGraph('512',  'grid',  '8x8x8', '1x1x1',   ),
        ProcessorGraph('256',  'torus', '16x16', '1x1x1x1'  ),
#        ProcessorGraph('1024', 'torus', '32x32', '1x1x1x1x1'),
        ProcessorGraph('512',  'torus', '8x8x8', '1x1x1'    ),
#        ProcessorGraph('256',  'tree',  '4x4x4x4'            , '64x16x4x1'          ),
#        ProcessorGraph('1024', 'tree',  '2x2x2x2x2x2x2x2x2x2', '1x1x1x1x1x1x1x1x1x1'),
#        ProcessorGraph('512',  'tree',  '8x8x8'              , '64x8x1'             ),
        ProcessorGraph('256', 'hypercube', '2x2x2x2x2x2x2x2', '1')
]

Aset = [#'initial',
        'random',
        #'drb',
        #'rcm',
        #'multisimple',
        #'multikahip',
        #'brandfass',
        #'brandfassL',
        #'hoefler',
        #'hoeflerL'
]

Cset = ['as-22july06',
        'as-skitter',
        'citationCiteseer',
        'coAuthorsCiteseer',
        'coAuthorsDBLP',
        'coPapersCiteseer',
        'coPapersDBLP',
        'email-EuAll',
        'loc-brightkite_edges',
        'loc-gowalla_edges',
        'p2p-Gnutella04',
        'PGPgiantcompo',
        'soc-Slashdot0902',
        'web-Google',
        'wiki-Talk'
]

Mset = ['fe_tooth',
        'fe_rotor',
        '598a',
        'fe_ocean',
        '144',
        'wave',
        'm14b',
        'auto'
]

statLines = []
stats = []

#for meshes
partDurations_256 = [
4.57077,
5.72353,
6.43707,
3.83098,
7.22077,
8.49545,
9.28305,
16.0473
]

#for meshes
partDurations_512 = [
6.99643,
8.93848,
11.3114,
5.51878,
12.4678,
13.0962,
15.4302,
23.268
]

#for extra_social
#partDurations_256 = [
#8.27015,
#2335.1,
#348.262,
#72.5771,
#233.959,
#356.104,
#865.136,
#4.77489,
#77.3281,
#483.792,
#5.5727,
#0.867453,
#266.771,
#109.612,
#1070.49
#]

#for extra_social
#partDurations_512 = [
#7.73432,
#4256.18,
#460.824,
#73.9928,
#298.871,
#344.216,
#1336.66,
#7.1657,
#48.6395,
#694.393,
#5.61128,
#1.18136,
#228.07,
#121.462,
#1284.41
#]

for procGraph in Pset:
    compare = None
    for mappingAlgo in Aset:
        graphCounter = 0;
        for commGraph in Mset:
            for seed in range(1,11):
                graphFile = graphDir + commGraph + '.graph'
                partFile = partDir + commGraph + '_partition_k=' + procGraph._K + '_seed=' + str(seed) + '_preconf=ecosocial'
                outFile = mapDir + commGraph + '_mapping_k=' + procGraph._K + '_type=' + procGraph._Typ + '_extensions=' + procGraph._Extensions + '_bandwidths=' + procGraph._Bandwidths + '_seed=' + str(seed) + '_preconf=ecosocial'
                process = subprocess.Popen([exePath, graphFile, '--input_partition=' + partFile, '--k=' + procGraph._K, '--proc_graph_type=' + procGraph._Typ, '--proc_graph_extensions=' + procGraph._Extensions, '--proc_graph_bandwidths=' + procGraph._Bandwidths, '--mapping_algo=' + mappingAlgo, '--filename_output=' + outFile, '--numberHierarchies=50'], stdout = subprocess.PIPE, stderr=subprocess.STDOUT)                    

                out, err = process.communicate()
                outString = out.decode("UTF-8")
                if outString.startswith( 'Error' ):
                    print (outString)
                    sys.exit(0)

                #find mapping duration
                keyword = outString.find("1:::Mapping duration:")
                start = keyword + 22
                end = outString.find("\n",keyword)
                if keyword == -1:
                    print (outString)
                    print (err)
                    sys.exit(0)
                duration1 = float(outString[start:end])
                keyword = outString.find("2:::Mapping duration:")
                start = keyword + 22
                end = outString.find("\n",keyword)
                if keyword == -1:
                    print (outString)
                    print (err)
                    sys.exit(0)
                duration2 = float(outString[start:end])

                #find maximum congestion
                keyword = outString.find("1:::Maximum congestion:")
                start = keyword + 24
                end = outString.find("\n",keyword)
                maxCongestion1 = int(outString[start:end])
                keyword = outString.find("2:::Maximum congestion:")
                start = keyword + 24
                end = outString.find("\n",keyword)
                maxCongestion2 = int(outString[start:end])

                #find maximum dilation
                keyword = outString.find("1:::Maximum dilation:")
                start = keyword + 22
                end = outString.find("\n",keyword)
                maxDilation1 = int(outString[start:end])
                keyword = outString.find("2:::Maximum dilation:")
                start = keyword + 22
                end = outString.find("\n",keyword)
                maxDilation2 = int(outString[start:end])

                #find average dilation
                keyword = outString.find("1:::Average dilation:")
                start = keyword + 22
                end = outString.find("\n",keyword)
                avgDilation1 = float(outString[start:end])
                keyword = outString.find("2:::Average dilation:")
                start = keyword + 22
                end = outString.find("\n",keyword)
                avgDilation2 = float(outString[start:end])

                #find Cut
                keyword = outString.find("1:::Cut is")
                start = keyword + 11
                end = outString.find("\n",keyword)
                Cut1 = float(outString[start:end])
                keyword = outString.find("2:::Cut is")
                start = keyword + 11
                end = outString.find("\n",keyword)
                Cut2 = float(outString[start:end])

                #find TCV
                keyword = outString.find("1:::TCV is")
                start = keyword + 11
                end = outString.find("\n",keyword)
                TCV1 = float(outString[start:end])
                keyword = outString.find("2:::TCV is")
                start = keyword + 11
                end = outString.find("\n",keyword)
                TCV2 = float(outString[start:end])

                if procGraph._K is '256':
                    duration1 = partDurations_256[graphCounter];
                if procGraph._K is '512':
                   duration1 = partDurations_512[graphCounter];
                #print duration1;

                line = StatLine(duration1, duration2, maxCongestion1, maxCongestion2, maxDilation1, maxDilation2, avgDilation1, avgDilation2, Cut1, Cut2, TCV1, TCV2)
                statLines.append(line)

            graphCounter = graphCounter + 1;
            #print(mappingAlgo + ' mapped ' + commGraph + ' on a ' + procGraph._Typ + ' with ' + procGraph._K + ' nodes.')
            stat = map_helper.analyze(statLines)
            statLines = []
            stats.append(stat)

        bigStat = map_helper.secondAnalyze(stats)
        stats = []
        bigStat.printString()

#        if compare is None:
#            compare = bigStat
#        else:
#            bigStat.compare(compare)
