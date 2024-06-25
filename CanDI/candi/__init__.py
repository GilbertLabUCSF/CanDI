from . import data

data = data.Data() #Global object data instantiated on import required for access by GeneQuery Objects

from .candi import (Gene, CellLine, Organelle, Cancer, CellLineCluster, GeneCluster)
