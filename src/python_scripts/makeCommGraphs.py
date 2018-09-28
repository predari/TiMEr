#_author__ = 'alex'
#manipulator Roland Glantz

import sys
import subprocess
import map_helper
import os
from graphs import ProcessorGraph

graphDir = os.path.expanduser('~/c+/KarMa/data/graphs/walshaw/')
partDir = os.path.expanduser('~/c+/KarMa/data/partitions/walshaw/')
commDir = os.path.expanduser('~/c+/KarMa/data/commGraphs/walshaw/')

exePath = os.path.expanduser('~/c+/KarMa/KarMa_sequential/bin/makeCommGraphs')

#Cset = ['as-22july06',
#        'as-skitter',
#        'citationCiteseer',
#        'coAuthorsCiteseer',
#        'coAuthorsDBLP',
#        'coPapersCiteseer',
#        'coPapersDBLP',
#        'email-EuAll',
#        'loc-brightkite_edges',
#        'loc-gowalla_edges',
#        'p2p-Gnutella04',
#        'PGPgiantcompo',
#        'soc-Slashdot0902',
#        'web-Google',
#        'wiki-Talk'
#]

Cset = [#'fe_tooth',
#        'fe_rotor',
#        '598a',
#        'fe_ocean',
#        '144',
#        'wave',
#        'm14b',
        'auto',
]

#Kset = [256, 512, 1024]
Kset = [1024]

for commGraph in Cset:
    for k in Kset:
        for seed in range(1,21):
            graphFile = graphDir + commGraph + '.graph'
            partFile = partDir + commGraph + '_partition_k=' + str(k) + '_seed=' + str(seed) + '_preconf=ecosocial'
            outFile = commDir + commGraph + '_k=' + str(k) + '_seed=' + str(seed) + '_preconf=ecosocial'
#            print (graphFile)
#            print (partFile)
#            print (outFile)
#            process = subprocess.Popen([exePath, graphFile, '--input_partition=' + partFile, '--k=' + str(k), '--filename_output=' + outFile], stdout = subprocess.PIPE, stderr=subprocess.STDOUT)                    
            process = subprocess.Popen([exePath, graphFile, '--k=' + str(k), '--filename_output=' + outFile], stdout = subprocess.PIPE, stderr=subprocess.STDOUT)                    
