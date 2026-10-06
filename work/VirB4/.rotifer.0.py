# IPython log file

# Initialization code
get_ipython().run_line_magic('matplotlib', 'inline')

# Python libraries
import os
import sys
import pickle
import networkx as nx
from io import StringIO
import matplotlib.pyplot as plt

# IPython
import ipywidgets as widgets
from qgridnext import show_grid
from IPython.display import display
output = widgets.Output()
plt.close()

# ROTIFER
import rotifer
from rotifer.db import ncbi
from rotifer.db import uniprot
from rotifer.genome import utils as rgu
from rotifer.core import functions as rcf
from rotifer.interval import utils as riu
from rotifer.genome import io as rgio
from rotifer.genome.data import NeighborhoodDF
from rotifer.devel.beta import hmmer as rdbh
from rotifer.devel.beta import sequence as rdbs
from rotifer.devel.alpha import gian_func as rdag

# Subroutines
def initialize(dummy):
    with output:
        output.clear_output()  # Clean previous click
        try:
            i = 0
            while os.path.exists(f".rotifer.{i}.py"):
                i = i + 1
            get_ipython().run_line_magic('logstart', f".rotifer.{i}.py")
            if os.path.exists(".rotifer"):
                rcf.restore_session(".rotifer")
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def save_session(dummy):
    with output:
        try:
            output.clear_output()
            rcf.save_session(".rotifer")
            print(f"Session saved!")
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

# IPython widgets: create
initButton = widgets.Button(description="Start", button_style="success")
saveButton = widgets.Button(description="Save", button_style="success")

# Connect callbacks
initButton.on_click(initialize)
saveButton.on_click(save_session)

# Bind and show
buttons = widgets.HBox([initButton, saveButton])
display(buttons)
display(output)
get_ipython().run_line_magic('who', '')
get_ipython().run_line_magic('logstate', '')
get_ipython().run_line_magic('connect_info', '')
get_ipython().run_line_magic('who', '')
from rotifer.devel.alpha import rodolfo as rdar
cagE = dict(data=rdbs.sequence("WP_000496009.1"))
cagE = dict(data=rdbs.sequence(["WP_000496009.1"]))
cagE
get_ipython().run_line_magic('who', '')
cagE["psi1"] = rdar.psiblast(cagE["data"], num_aln=5000, aln=False, db="/scratch/global/databases/fadb/uniprot/uniref50")
type(cagE["psi1"])
cagE["psi1"][0]
cagE["psi1"][0].columns
page(cagE["psi1"][1])
page(cagE["psi1"][1])
