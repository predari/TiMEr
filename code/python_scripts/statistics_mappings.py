#_author__ = 'alex'
#manipulator Roland Glantz

import sys
import subprocess
import map_helper
import os
from graphs import ProcessorGraph

class StatLine(object):
    def __init__(self, duration, maxCon, maxDil, avgDil):
        self._Duration = duration
        self._MaxCon = maxCon
        self._MaxDil = maxDil
        self._AvgDil = avgDil
        
    def toString(self):
        return ("Mapping duration: {duration}, Maximum congestion: {maxCon}, Maximum dilation: {maxDil}, Average dilation: {avgDil}".format(duration=self._Duration, maxCon=self._MaxCon, maxDil=self._MaxDil, avgDil=self._AvgDil))

graphDir = os.path.expanduser('~/c+/KarMa/data/graphs/extra_social/')
partDir = os.path.expanduser('~/c+/KarMa/data/partitions/extra_social/')
mapDir = os.path.expanduser('~/c+/KarMa/data/mappings/extra_social/')

exePath = os.path.expanduser('~/c+/KarMa/KarMa_sequential/bin/mapScenario')

Pset = [ProcessorGraph('256',  'grid',  '16x16', '1x1x1x1'  ),
        ProcessorGraph('1024', 'grid',  '32x32', '1x1x1x1x1'),
        ProcessorGraph('512',  'grid',  '8x8x8', '1x1x1',   ),
        ProcessorGraph('256',  'torus', '16x16', '1x1x1x1'  ),
        ProcessorGraph('1024', 'torus', '32x32', '1x1x1x1x1'),
        ProcessorGraph('512',  'torus', '8x8x8', '1x1x1'    ),
        ProcessorGraph('256',  'tree',  '4x4x4x4'            , '64x16x4x1'          ),
        ProcessorGraph('1024', 'tree',  '2x2x2x2x2x2x2x2x2x2', '1x1x1x1x1x1x1x1x1x1'),
        ProcessorGraph('512',  'tree',  '8x8x8'              , '64x8x1'             )
]

Mset = ['initial',
        'random',
        'drb',
        'rcm',
        'multisimple',
        'multikahip',
        'brandfass',
        'brandfassL',
        'hoefler',
        'hoeflerL'
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

statLines = []
stats = []

for procGraph in Pset:
    compare = None
    for mappingAlgo in Mset:
        for commGraph in Cset:
            for seed in range(1,21):
                graphFile = graphDir + commGraph + '.graph'
                partFile = partDir + commGraph + '_partition_k=' + procGraph._K + '_seed=' + str(seed) + '_preconf=ecosocial'
                outFile = mapDir + commGraph + '_mapping_k=' + procGraph._K + '_type=' + procGraph._Typ + '_extensions=' + procGraph._Extensions + '_bandwidths=' + procGraph._Bandwidths + '_seed=' + str(seed) + '_preconf=ecosocial'
                process = subprocess.Popen([exePath, graphFile, '--input_partition=' + partFile, '--k=' + procGraph._K, '--proc_graph_type=' + procGraph._Typ, '--proc_graph_extensions=' + procGraph._Extensions, '--proc_graph_bandwidths=' + procGraph._Bandwidths, '--mapping_algo=' + mappingAlgo, '--filename_output=' + outFile], stdout = subprocess.PIPE, stderr=subprocess.STDOUT)                    

                out, err = process.communicate()
                outString = out.decode("UTF-8")
                if outString.startswith( 'Error' ):
                    print (outString)
                    sys.exit(0)

                #find mapping duration
                keyword = outString.find("Mapping duration:")
                start = keyword + 18
                end = outString.find("\n",keyword)
                if keyword == -1:
                    print (outString)
                    print (err)
                    sys.exit(0)
                duration = float(outString[start:end])

                #find maximum congestion
                keyword = outString.find("Maximum congestion:")
                start = keyword + 20
                end = outString.find("\n",keyword)
                maxCongestion = int(outString[start:end])

                #find maximum dilation
                keyword = outString.find("Maximum dilation:")
                start = keyword + 18
                end = outString.find("\n",keyword)
                maxDilation = int(outString[start:end])

                #find average dilation
                keyword = outString.find("Average dilation:")
                start = keyword + 18
                end = outString.find("\n",keyword)
                avgDilation = float(outString[start:end])
                #print(avgDilation)

                line = StatLine(duration, maxCongestion, maxDilation, avgDilation)
                statLines.append(line)

            #print(mappingAlgo + ' mapped ' + commGraph + ' on a ' + procGraph._Typ + ' with ' + procGraph._K + ' nodes.')
            stat = map_helper.analyze(statLines)
            statLines = []
            stats.append(stat)

        bigStat = map_helper.secondAnalyze(stats)
        stats = []
        bigStat.printString()

        if compare is None:
            compare = bigStat
        else:
            bigStat.compare(compare)
