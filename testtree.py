from ete3 import TreeStyle, RectFace, faces, Tree, NodeStyle
from ete3 import NCBITaxa

from Bio import Phylo
ncbi = NCBITaxa()
tree = ncbi.get_topology([2, 33208], intermediate_nodes=True)
print(tree.get_ascii(attributes=["sci_name"]))