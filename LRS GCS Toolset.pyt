import arcpy, IdentifySparseVertices
import importlib
importlib.reload(IdentifySparseVertices)

from IdentifySparseVertices import *

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "LRS GCS Toolset"
        self.alias = "LRSGCSToolset"

        # List of tool classes associated with this toolbox
        self.tools = [IdentifySparseVertices]

