__author__ = 'alex'
#manipulator Roland Glantz

class ProcessorGraph(object):
    def __init__(self, k, typ, extensions, bandwidths):
        self._K = k
        self._Typ = typ
        self._Extensions = extensions
        self._Bandwidths = bandwidths

class CommunicationGraph(object):
    def __init__(self, graph, suffix):
        self._Graph = graph
        self._Suffix = suffix
