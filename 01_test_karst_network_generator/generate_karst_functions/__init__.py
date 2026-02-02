# mypackages/__init__.py
"""
Package pour la génération et la manipulation de réseaux karstiques.
"""

from .function_generate_karst import KarstNetworks
from .function_graf_simplification import  GRF_objt, get_position_center, get_sgs_values, simplify_graph, export_conduits, plot_conduits
from .function_write_inputs_cfp import write_kexch, write_geoheight, write_node_heads, write_pipe_info, write_network_info

#définit explicitement ce qui sera importé si quelqu’un fait: from mypackages import *
__all__ = ["KarstNetworks", "GRF_objt", "get_position_center", "get_sgs_values", "simplify_graph", "export_conduits", "plot_conduits","write_geoheight", "write_kexch",
           "write_node_heads", "write_pipe_info", "write_network_info"]
