# IPython log file

# Run this cell to activate or reset the panel buttons

# Python libraries
import os
import sys
import pickle
import networkx as nx
from io import StringIO
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
plt.close()

# IPython
import ipywidgets as widgets
from qgridnext import show_grid
from IPython.display import display

# ROTIFER
import rotifer
from rotifer.db import ncbi
from rotifer.db import uniprot
from rotifer.genome import utils as rgu
from rotifer.core import functions as rcf
from rotifer.interval import utils as riu
from rotifer.genome.data import NeighborhoodDF
from rotifer.devel.beta import hmmer as rdbh
from rotifer.devel.beta import sequence as rdbs
from rotifer.devel.alpha import gian_func as rdag

# Subroutines
def new_id():
    from datetime import datetime
    return datetime.now().strftime("%Y%m%d%H%M%S")

def new_session(widget):
    sesText.value = new_id()

def get_session_path():
    from pathlib import Path
    return Path(".rotifer", sesText.value)

def on_dropdown_change(change):
  if change['type'] == 'change' and change['name'] == 'value':
    sesText.value = change["new"]

def show_who(widget):
    with output:
        try:
            output.clear_output()
            get_ipython().run_line_magic('who','')
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def show_logfile(widget):
    with output:
        try:
            output.clear_output()
            print("Session ID     :",sesText.value)
            get_ipython().run_line_magic('logstate','')
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def show_connection(widget):
    with output:
        try:
            output.clear_output()
            print(get_jupyter_server_list_output())
            print("\nThis kernel's connection info:")
            get_ipython().run_line_magic('connect_info','')
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def clean_output(widget):
    with output:
        try:
            output.clear_output()
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def get_jupyter_server_list_output():
    server_loaders = []

    # Check for modern or legacy server utility packages
    try:
        from jupyter_server.serverapp import list_running_servers
        server_loaders.append(list_running_servers)
    except ImportError:
        pass
    try:
        from notebook.notebookapp import list_running_servers
        server_loaders.append(list_running_servers)
    except ImportError:
        pass

    active_servers = []
    for load_servers in server_loaders:
        for server in load_servers():
            base_url = server.get("url", "")
            token = server.get("token", "")
            root_dir = server.get("root_dir", server.get("notebook_dir", "Unknown Directory"))

            # Format URL with token authentication key if it exists
            if token:
                full_url = f"{base_url}?token={token}"
            else:
                full_url = base_url

            active_servers.append(f"{full_url} :: {root_dir}")

    # Standard terminal layout mirroring the CLI tool
    if active_servers:
        output_text = "Currently running servers:\n" + "\n".join(active_servers)
    else:
        output_text = (
            "Currently running servers:\n"
            "No active standalone web servers found.\n"
            "(If you are using VS Code Native layout, it executes via direct Python subprocesses.)"
        )
    return output_text

# Bind and show
def display_panel():
    buttons = widgets.HBox([saveButton, clsButton, whoButton, logButton, connButton, initButton])
    display(buttons)
    buttons2 = widgets.HBox([sesText,sesCombo,newButton])
    display(buttons2)
    display(output)

# Get data
def get_saved_sessions():
    from pathlib import Path
    sessions = []
    if os.path.exists(".rotifer"):
        if os.path.isdir(".rotifer"):
            sessions = [ item for item in Path('.rotifer').iterdir() if item.is_dir() and (item / "session.dat").exists() ]
        else:
            print(".rotifer exists but it is not a directory!!!")
    sessions = [ item.name for item in sorted(sessions, key=lambda x: (x / "session.dat").stat().st_mtime, reverse=True) ]
    return sessions

# Operate on data: load, save, etc
def initialize(widget):
    with output:
        output.clear_output()  # Clean previous click
        try:
            sesdir = get_session_path()
            if not sesdir.exists():
                print(f"Directory for session {sesText.value} created at {str(sesdir)}")
                sesdir.mkdir(parents=True, exist_ok=True)
            if get_ipython().logger.logfile == None:
                i = 0
                logfile = sesdir / f"log.{i}.py"
                while logfile.exists():
                    i = i + 1
                    logfile = sesdir / f"log.{i}.py"
                get_ipython().run_line_magic('logstart',str(logfile))
                del(logfile)
            rcf.restore_session(sesdir / "session")
            get_ipython().run_line_magic('who','')
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def save_session(widget):
    with output:
        try:
            output.clear_output()
            path = get_session_path()
            path.mkdir(parents=True, exist_ok=True)
            rcf.save_session(path / "session")
            sesCombo.options = get_saved_sessions()
            sesCombo.value = sesText.value
            print(f"Session saved!")
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

# Initialization before UI generation
_saved_sessions = get_saved_sessions()
if "_ROTIFER_SESSION_ID" not in globals():
    _ROTIFER_SESSION_ID = new_id()

# IPython widgets: create
output = widgets.Output()
newButton = widgets.Button(
    description="New session",
    button_style="success",
    icon="plus",
    tooltip="Create a new session"
)
initButton = widgets.Button(
    description="Reload",
    button_style="success",
    icon="refresh",
    tooltip="Replace data in memory with the last saved session"
)
saveButton = widgets.Button(
    description="Save",
    button_style="success",
    icon="floppy-o",
    tooltip="Save current session data"
)
whoButton  = widgets.Button(
    description="List",
    button_style="info",
    icon="list",
    tooltip="List variables, modules and functions"
)
logButton  = widgets.Button(
    description="Logstate",
    button_style="info",
    icon="history",
    tooltip="Show logging file name and info"
)
connButton = widgets.Button(
    description="Connection",
    button_style="info",
    icon="server",
    tooltip="Click to print the list of jupyter servers with auth tokens"
)
clsButton  = widgets.Button(
    description="Clear",
    button_style="info",
    icon="circle-thin",
    tooltip="Clear output"
)
sesCombo = widgets.Dropdown(
    options=_saved_sessions,
    description='Saved sessions:',
    disabled=False,
    style={'description_width': 'initial'},
    icon="check"
)
sesText = widgets.Text(
    placeholder='Rename current session',
    description='Session ID:',
    value=_saved_sessions[0] if len(_saved_sessions) > 0 else _ROTIFER_SESSION_ID,
    disabled=False
)

# Connect callbacks
newButton.on_click(new_session)
initButton.on_click(initialize)
saveButton.on_click(save_session)
whoButton.on_click(show_who)
logButton.on_click(show_logfile)
connButton.on_click(show_connection)
clsButton.on_click(clean_output)

# Attach the observer to the dropdown
sesCombo.observe(on_dropdown_change, names='value')

# Main
display_panel()
get_ipython().system('ls -lhrt')
get_ipython().system('ls -a')
get_ipython().system('less .rotifer.0.py')
cagE = dict(data=rdbs.sequence(["WP_000496009.1"]))
cagE["psi1"] = rdar.psiblast(cagE["data"], num_aln=5000, aln=False, db="/scratch/global/databases/fadb/uniprot/uniref50")
# Run this cell to activate or reset the panel buttons

# Python libraries
import os
import sys
import pickle
import networkx as nx
from io import StringIO
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
plt.close()

# IPython
import ipywidgets as widgets
from qgridnext import show_grid
from IPython.display import display

# ROTIFER
import rotifer
from rotifer.db import ncbi
from rotifer.db import uniprot
from rotifer.genome import utils as rgu
from rotifer.core import functions as rcf
from rotifer.interval import utils as riu
from rotifer.genome.data import NeighborhoodDF
from rotifer.devel.beta import hmmer as rdbh
from rotifer.devel.beta import sequence as rdbs
from rotifer.devel.alpha import gian_func as rdag
from rotifer.devel.alpha import epsoares as rdae
from rotifer.devel.alpha import rodolfo as rdar

# Subroutines
def new_id():
    from datetime import datetime
    return datetime.now().strftime("%Y%m%d%H%M%S")

def new_session(widget):
    sesText.value = new_id()

def get_session_path():
    from pathlib import Path
    return Path(".rotifer", sesText.value)

def on_dropdown_change(change):
  if change['type'] == 'change' and change['name'] == 'value':
    sesText.value = change["new"]

def show_who(widget):
    with output:
        try:
            output.clear_output()
            get_ipython().run_line_magic('who','')
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def show_logfile(widget):
    with output:
        try:
            output.clear_output()
            print("Session ID     :",sesText.value)
            get_ipython().run_line_magic('logstate','')
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def show_connection(widget):
    with output:
        try:
            output.clear_output()
            print(get_jupyter_server_list_output())
            print("\nThis kernel's connection info:")
            get_ipython().run_line_magic('connect_info','')
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def clean_output(widget):
    with output:
        try:
            output.clear_output()
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def get_jupyter_server_list_output():
    server_loaders = []

    # Check for modern or legacy server utility packages
    try:
        from jupyter_server.serverapp import list_running_servers
        server_loaders.append(list_running_servers)
    except ImportError:
        pass
    try:
        from notebook.notebookapp import list_running_servers
        server_loaders.append(list_running_servers)
    except ImportError:
        pass

    active_servers = []
    for load_servers in server_loaders:
        for server in load_servers():
            base_url = server.get("url", "")
            token = server.get("token", "")
            root_dir = server.get("root_dir", server.get("notebook_dir", "Unknown Directory"))

            # Format URL with token authentication key if it exists
            if token:
                full_url = f"{base_url}?token={token}"
            else:
                full_url = base_url

            active_servers.append(f"{full_url} :: {root_dir}")

    # Standard terminal layout mirroring the CLI tool
    if active_servers:
        output_text = "Currently running servers:\n" + "\n".join(active_servers)
    else:
        output_text = (
            "Currently running servers:\n"
            "No active standalone web servers found.\n"
            "(If you are using VS Code Native layout, it executes via direct Python subprocesses.)"
        )
    return output_text

# Bind and show
def display_panel():
    buttons = widgets.HBox([saveButton, clsButton, whoButton, logButton, connButton, initButton])
    display(buttons)
    buttons2 = widgets.HBox([sesText,sesCombo,newButton])
    display(buttons2)
    display(output)

# Get data
def get_saved_sessions():
    from pathlib import Path
    sessions = []
    if os.path.exists(".rotifer"):
        if os.path.isdir(".rotifer"):
            sessions = [ item for item in Path('.rotifer').iterdir() if item.is_dir() and (item / "session.dat").exists() ]
        else:
            print(".rotifer exists but it is not a directory!!!")
    sessions = [ item.name for item in sorted(sessions, key=lambda x: (x / "session.dat").stat().st_mtime, reverse=True) ]
    return sessions

# Operate on data: load, save, etc
def initialize(widget):
    with output:
        output.clear_output()  # Clean previous click
        try:
            sesdir = get_session_path()
            if not sesdir.exists():
                print(f"Directory for session {sesText.value} created at {str(sesdir)}")
                sesdir.mkdir(parents=True, exist_ok=True)
            if get_ipython().logger.logfile == None:
                i = 0
                logfile = sesdir / f"log.{i}.py"
                while logfile.exists():
                    i = i + 1
                    logfile = sesdir / f"log.{i}.py"
                get_ipython().run_line_magic('logstart',str(logfile))
                del(logfile)
            rcf.restore_session(sesdir / "session")
            get_ipython().run_line_magic('who','')
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def save_session(widget):
    with output:
        try:
            output.clear_output()
            path = get_session_path()
            path.mkdir(parents=True, exist_ok=True)
            rcf.save_session(path / "session")
            sesCombo.options = get_saved_sessions()
            sesCombo.value = sesText.value
            print(f"Session saved!")
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

# Initialization before UI generation
_saved_sessions = get_saved_sessions()
if "_ROTIFER_SESSION_ID" not in globals():
    _ROTIFER_SESSION_ID = new_id()

# IPython widgets: create
output = widgets.Output()
newButton = widgets.Button(
    description="New session",
    button_style="success",
    icon="plus",
    tooltip="Create a new session"
)
initButton = widgets.Button(
    description="Reload",
    button_style="success",
    icon="refresh",
    tooltip="Replace data in memory with the last saved session"
)
saveButton = widgets.Button(
    description="Save",
    button_style="success",
    icon="floppy-o",
    tooltip="Save current session data"
)
whoButton  = widgets.Button(
    description="List",
    button_style="info",
    icon="list",
    tooltip="List variables, modules and functions"
)
logButton  = widgets.Button(
    description="Logstate",
    button_style="info",
    icon="history",
    tooltip="Show logging file name and info"
)
connButton = widgets.Button(
    description="Connection",
    button_style="info",
    icon="server",
    tooltip="Click to print the list of jupyter servers with auth tokens"
)
clsButton  = widgets.Button(
    description="Clear",
    button_style="info",
    icon="circle-thin",
    tooltip="Clear output"
)
sesCombo = widgets.Dropdown(
    options=_saved_sessions,
    description='Saved sessions:',
    disabled=False,
    style={'description_width': 'initial'},
    icon="check"
)
sesText = widgets.Text(
    placeholder='Rename current session',
    description='Session ID:',
    value=_saved_sessions[0] if len(_saved_sessions) > 0 else _ROTIFER_SESSION_ID,
    disabled=False
)

# Connect callbacks
newButton.on_click(new_session)
initButton.on_click(initialize)
saveButton.on_click(save_session)
whoButton.on_click(show_who)
logButton.on_click(show_logfile)
connButton.on_click(show_connection)
clsButton.on_click(clean_output)

# Attach the observer to the dropdown
sesCombo.observe(on_dropdown_change, names='value')

# Main
display_panel()
cagE["psi1"] = rdar.psiblast(cagE["data"], num_aln=5000, aln=False, db="/scratch/global/databases/fadb/uniprot/uniref50")
get_ipython().run_line_magic('who', '')
cagE['psi1'][0]
# Run this cell to activate or reset the panel buttons

# Python libraries
import os
import sys
import pickle
import networkx as nx
from io import StringIO
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
plt.close()

# IPython
import ipywidgets as widgets
from qgridnext import show_grid
from IPython.display import display

# ROTIFER
import rotifer
from rotifer.db import ncbi
from rotifer.db import uniprot
from rotifer.genome import utils as rgu
from rotifer.core import functions as rcf
from rotifer.interval import utils as riu
from rotifer.genome.data import NeighborhoodDF
from rotifer.devel.beta import hmmer as rdbh
from rotifer.devel.beta import sequence as rdbs
from rotifer.devel.alpha import gian_func as rdag
from rotifer.devel.alpha import epsoares as rdae
from rotifer.devel.alpha import rodolfo as rdar

# Subroutines
def new_id():
    from datetime import datetime
    return datetime.now().strftime("%Y%m%d%H%M%S")

def new_session(widget):
    sesText.value = new_id()

def get_session_path():
    from pathlib import Path
    return Path(".rotifer", sesText.value)

def on_dropdown_change(change):
  if change['type'] == 'change' and change['name'] == 'value':
    sesText.value = change["new"]

def show_who(widget):
    with output:
        try:
            output.clear_output()
            get_ipython().run_line_magic('who','')
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def show_logfile(widget):
    with output:
        try:
            output.clear_output()
            print("Session ID     :",sesText.value)
            get_ipython().run_line_magic('logstate','')
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def show_connection(widget):
    with output:
        try:
            output.clear_output()
            print(get_jupyter_server_list_output())
            print("\nThis kernel's connection info:")
            get_ipython().run_line_magic('connect_info','')
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def clean_output(widget):
    with output:
        try:
            output.clear_output()
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def get_jupyter_server_list_output():
    server_loaders = []

    # Check for modern or legacy server utility packages
    try:
        from jupyter_server.serverapp import list_running_servers
        server_loaders.append(list_running_servers)
    except ImportError:
        pass
    try:
        from notebook.notebookapp import list_running_servers
        server_loaders.append(list_running_servers)
    except ImportError:
        pass

    active_servers = []
    for load_servers in server_loaders:
        for server in load_servers():
            base_url = server.get("url", "")
            token = server.get("token", "")
            root_dir = server.get("root_dir", server.get("notebook_dir", "Unknown Directory"))

            # Format URL with token authentication key if it exists
            if token:
                full_url = f"{base_url}?token={token}"
            else:
                full_url = base_url

            active_servers.append(f"{full_url} :: {root_dir}")

    # Standard terminal layout mirroring the CLI tool
    if active_servers:
        output_text = "Currently running servers:\n" + "\n".join(active_servers)
    else:
        output_text = (
            "Currently running servers:\n"
            "No active standalone web servers found.\n"
            "(If you are using VS Code Native layout, it executes via direct Python subprocesses.)"
        )
    return output_text

# Bind and show
def display_panel():
    buttons = widgets.HBox([saveButton, clsButton, whoButton, logButton, connButton, initButton])
    display(buttons)
    buttons2 = widgets.HBox([sesText,sesCombo,newButton])
    display(buttons2)
    display(output)

# Get data
def get_saved_sessions():
    from pathlib import Path
    sessions = []
    if os.path.exists(".rotifer"):
        if os.path.isdir(".rotifer"):
            sessions = [ item for item in Path('.rotifer').iterdir() if item.is_dir() and (item / "session.dat").exists() ]
        else:
            print(".rotifer exists but it is not a directory!!!")
    sessions = [ item.name for item in sorted(sessions, key=lambda x: (x / "session.dat").stat().st_mtime, reverse=True) ]
    return sessions

# Operate on data: load, save, etc
def initialize(widget):
    with output:
        output.clear_output()  # Clean previous click
        try:
            sesdir = get_session_path()
            if not sesdir.exists():
                print(f"Directory for session {sesText.value} created at {str(sesdir)}")
                sesdir.mkdir(parents=True, exist_ok=True)
            if get_ipython().logger.logfile == None:
                i = 0
                logfile = sesdir / f"log.{i}.py"
                while logfile.exists():
                    i = i + 1
                    logfile = sesdir / f"log.{i}.py"
                get_ipython().run_line_magic('logstart',str(logfile))
                del(logfile)
            rcf.restore_session(sesdir / "session")
            get_ipython().run_line_magic('who','')
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

def save_session(widget):
    with output:
        try:
            output.clear_output()
            path = get_session_path()
            path.mkdir(parents=True, exist_ok=True)
            rcf.save_session(path / "session")
            sesCombo.options = get_saved_sessions()
            sesCombo.value = sesText.value
            print(f"Session saved!")
        except Exception as e:
            print(f"ERROR: {type(e).__name__}: {e}")

# Initialization before UI generation
_saved_sessions = get_saved_sessions()
if "_ROTIFER_SESSION_ID" not in globals():
    _ROTIFER_SESSION_ID = new_id()

# IPython widgets: create
output = widgets.Output()
newButton = widgets.Button(
    description="New session",
    button_style="success",
    icon="plus",
    tooltip="Create a new session"
)
initButton = widgets.Button(
    description="Reload",
    button_style="success",
    icon="refresh",
    tooltip="Replace data in memory with the last saved session"
)
saveButton = widgets.Button(
    description="Save",
    button_style="success",
    icon="floppy-o",
    tooltip="Save current session data"
)
whoButton  = widgets.Button(
    description="List",
    button_style="info",
    icon="list",
    tooltip="List variables, modules and functions"
)
logButton  = widgets.Button(
    description="Logstate",
    button_style="info",
    icon="history",
    tooltip="Show logging file name and info"
)
connButton = widgets.Button(
    description="Connection",
    button_style="info",
    icon="server",
    tooltip="Click to print the list of jupyter servers with auth tokens"
)
clsButton  = widgets.Button(
    description="Clear",
    button_style="info",
    icon="circle-thin",
    tooltip="Clear output"
)
sesCombo = widgets.Dropdown(
    options=_saved_sessions,
    description='Saved sessions:',
    disabled=False,
    style={'description_width': 'initial'},
    icon="check"
)
sesText = widgets.Text(
    placeholder='Rename current session',
    description='Session ID:',
    value=_saved_sessions[0] if len(_saved_sessions) > 0 else _ROTIFER_SESSION_ID,
    disabled=False
)

# Connect callbacks
newButton.on_click(new_session)
initButton.on_click(initialize)
saveButton.on_click(save_session)
whoButton.on_click(show_who)
logButton.on_click(show_logfile)
connButton.on_click(show_connection)
clsButton.on_click(clean_output)

# Attach the observer to the dropdown
sesCombo.observe(on_dropdown_change, names='value')

# Main
display_panel()
rcf.print_everything()
rpf.print_everything()
from rotifer.pandas import functions as rfp
from rotifer.pandas import functions as rpf
del rfp
rpf.print_everything()
cagE['psi1'][0]
get_ipython().run_line_magic('page', '')
cagE["psi1"][0]
from IPython.core.page import page
page(cagE['psi1'][1])
cagE['psi1'][0].head(5)
cagE['psi1'][0]['evalue'].max()
page(cagE['psi1'][1])
cagE['psi1'][0].hit.nunique()
get_ipython().run_line_magic('pinfo2', 'rdbs.sequence')
cagE['psi1seqs'] = rdbs.sequence(cagE['psi1'][0].hit.drop_duplicates().to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['psi1seqs'].align()._view_sequence()
cagE['psi1seqs'].add_cluster(cascade=[(100,100),(80,70),(80,50),(80,30),(0,0)], inplace=True)
cagE['psi1seqs'].nunique()
cagE['psi1seqs'].df.nunique()
cagE['psi1seqs'].df.c0i0.value_counts()
cagE['psi1'][0].head()
cagE['psi1seqs'].df.query('id == "UniRef50_Q48252"').c0i0
# We should investigate evalue and bitscore
# The result consists into a tuple, 0 is the hits table, 1 is the out of psiblast

cagE["psi1"][0].head(10)
page(cagE["psi1"][1]).head(10)
page(cagE["psi1"][1])
#page(cagE["psi1"][1])
cagE["psi1"][0].evalue.hist()
# Let's see usefull statistics!

import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

cagE["psi1"][0].evalue.apply(
    lambda x: -np.log10(max(x, 1e-300))
).hist(ax=axes[0], bins=30)

cagE["psi1"][0].bitscore.hist(ax=axes[1], bins=30)
cagE["psi1"][0].length.hist(ax=axes[2], bins=30)

for ax, title in zip(axes, ["-log10(E-value)", "Bit score", "Alignment length"]):
    ax.set_title(title)
    ax.set_xlabel(title)
    ax.set_ylabel("Frequency")

plt.tight_layout()
plt.show()
import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

cagE["psi1"][0].evalue.apply(
    lambda x: -np.log10(max(x, 1e-300))
).hist(ax=axes[0], bins=30)

cagE["psi1"][0].bitscore.hist(ax=axes[1], bins=30)
cagE["psi1"][0].length.hist(ax=axes[2], bins=30)

for ax, title in zip(axes, ["-log10(E-value)", "Bit score", "Alignment length"]):
    ax.set_title(title)
    ax.set_xlabel(title)
    ax.set_ylabel("Frequency")

fig.text(
    0.99, 0.99,
    f"Unique hits: {cagE['psi1'][0].hit.nunique()}",
    ha="right", va="top", fontsize=11
)

plt.tight_layout(rect=[0, 0, 1, 0.92])
plt.show()
import matplotlib.pyplot as plt

cols = ["c100i100", "c80i70", "c80i50", "c80i30", "c0i0"]
counts = clusters[cols].nunique()

fig, ax = plt.subplots(figsize=(9, 5))

bars = ax.bar(counts.index, counts.values)

ax.bar_label(bars, padding=3)
ax.set_xlabel("Clustering threshold")
ax.set_ylabel("Number of unique clusters")
ax.set_title("Cluster distribution")
ax.set_ylim(0, counts.max() * 1.1)

plt.tight_layout()
plt.show()
import matplotlib.pyplot as plt

cols = ["c100i100", "c80i70", "c80i50", "c80i30", "c0i0"]
counts = cagE['psi1seqs']df[cols].nunique()

fig, ax = plt.subplots(figsize=(9, 5))

bars = ax.bar(counts.index, counts.values)

ax.bar_label(bars, padding=3)
ax.set_xlabel("Clustering threshold")
ax.set_ylabel("Number of unique clusters")
ax.set_title("Cluster distribution")
ax.set_ylim(0, counts.max() * 1.1)

plt.tight_layout()
plt.show()
import matplotlib.pyplot as plt

cols = ["c100i100", "c80i70", "c80i50", "c80i30", "c0i0"]
counts = cagE['psi1seqs'].df[cols].nunique()

fig, ax = plt.subplots(figsize=(9, 5))

bars = ax.bar(counts.index, counts.values)

ax.bar_label(bars, padding=3)
ax.set_xlabel("Clustering threshold")
ax.set_ylabel("Number of unique clusters")
ax.set_title("Cluster distribution")
ax.set_ylim(0, counts.max() * 1.1)

plt.tight_layout()
plt.show()
import matplotlib.pyplot as plt
import math

cluster_cols = ["c100i100", "c80i70", "c80i50", "c80i30", "c0i0"]
top_hit = cagE["psi1"][0].iloc[0]["hit"]

row = clusters.loc[clusters["id"] == top_hit]
if row.empty:
    raise ValueError(f"Top hit {top_hit!r} was not found in clusters['id'].")

n = len(cluster_cols)
ncols = 3
nrows = math.ceil(n / ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
axes = axes.flatten()

for ax, col in zip(axes, cluster_cols):
    cluster_sizes = clusters[col].value_counts()
    size_distribution = cluster_sizes.value_counts().sort_index()
    # index = cluster size, value = number of clusters with that size

    top_cluster = row.iloc[0][col]
    top_cluster_size = cluster_sizes[top_cluster]

    bars = ax.bar(size_distribution.index, size_distribution.values)

    # highlight the bar corresponding to the top-hit cluster size
    if top_cluster_size in size_distribution.index:
        idx = list(size_distribution.index).index(top_cluster_size)
        bars[idx].set_hatch("///")
        bars[idx].set_edgecolor("red")
        bars[idx].set_linewidth(1.5)

    ax.set_title(col)
    ax.set_xlabel("Cluster size")
    ax.set_ylabel("Number of clusters")

    ax.text(
        0.98, 0.98,
        f"top-hit cluster size: {top_cluster_size}",
        transform=ax.transAxes,
        ha="right", va="top"
    )

for ax in axes[len(cluster_cols):]:
    ax.set_visible(False)

plt.tight_layout()
plt.show()
import matplotlib.pyplot as plt

# Best PSI-BLAST hit
top_hit = cagE["psi1"][0].sort_values(
    ["evalue", "bitscore"],
    ascending=[True, False]
).iloc[0]["hit"]

seqs = cagE["psi1seqs"].df  # Reference, not a new DataFrame

fig, axes = plt.subplots(2, 3, figsize=(15, 9))

for ax, col in zip(axes.flat, cols):

    # Number of sequences per cluster
    sizes = seqs[col].value_counts()

    # Distribution of cluster sizes
    distribution = sizes.value_counts().sort_index()

    # Cluster containing the best PSI-BLAST hit
    cluster_id = seqs.loc[seqs["id"] == top_hit, col].iloc[0]
    top_size = sizes[cluster_id]

    # Plot
    bars = ax.bar(distribution.index, distribution.values)

    # Stripe the corresponding size
    idx = distribution.index.get_loc(top_size)
    bars[idx].set_hatch("///")
    bars[idx].set_edgecolor("red")
    bars[idx].set_linewidth(1.5)

    ax.set_title(f"{col} (n={len(sizes)} clusters)")
    ax.set_xlabel("Cluster size")
    ax.set_ylabel("Number of clusters")

    ax.text(
        0.98, 0.98,
        f"Top-hit cluster: {top_size} sequences",
        transform=ax.transAxes,
        ha="right", va="top"
    )

    ax.set_yscale("log")

axes.flat[-1].axis("off")

plt.tight_layout()
plt.show()
import matplotlib.pyplot as plt
cols = ["c80i70", "c80i50", "c80i30", "c0i0"]
# Best PSI-BLAST hit
top_hit = cagE["psi1"][0].sort_values(
    ["evalue", "bitscore"],
    ascending=[True, False]
).iloc[0]["hit"]

seqs = cagE["psi1seqs"].df  # Reference, not a new DataFrame

fig, axes = plt.subplots(2, 3, figsize=(15, 9))

for ax, col in zip(axes.flat, cols):

    # Number of sequences per cluster
    sizes = seqs[col].value_counts()

    # Distribution of cluster sizes
    distribution = sizes.value_counts().sort_index()

    # Cluster containing the best PSI-BLAST hit
    cluster_id = seqs.loc[seqs["id"] == top_hit, col].iloc[0]
    top_size = sizes[cluster_id]

    # Plot
    bars = ax.bar(distribution.index, distribution.values)

    # Stripe the corresponding size
    idx = distribution.index.get_loc(top_size)
    bars[idx].set_hatch("///")
    bars[idx].set_edgecolor("red")
    bars[idx].set_linewidth(1.5)

    ax.set_title(f"{col} (n={len(sizes)} clusters)")
    ax.set_xlabel("Cluster size")
    ax.set_ylabel("Number of clusters")

    ax.text(
        0.98, 0.98,
        f"Top-hit cluster: {top_size} sequences",
        transform=ax.transAxes,
        ha="right", va="top"
    )

    ax.set_yscale("log")

axes.flat[-1].axis("off")

plt.tight_layout()
plt.show()
import numpy as np
import matplotlib.pyplot as plt

top_hit = cagE["psi1"][0].sort_values(
    ["evalue", "bitscore"],
    ascending=[True, False]
).iloc[0]["hit"]

seqs = cagE["psi1seqs"].df

fig, axes = plt.subplots(2, 3, figsize=(15, 9))

for ax, col in zip(axes.flat, cols):

    # Cluster sizes, sorted largest to smallest
    sizes = seqs[col].value_counts()

    # Identify the cluster containing the top hit
    top_cluster = seqs.loc[seqs["id"] == top_hit, col].iloc[0]
    top_rank = sizes.index.get_loc(top_cluster)
    top_size = sizes.iloc[top_rank]

    # One bar per cluster
    bars = ax.bar(np.arange(1, len(sizes) + 1), sizes.values, width=1)

    # Highlight the specific cluster containing the top hit
    bars[top_rank].set_hatch("///")
    bars[top_rank].set_edgecolor("red")
    bars[top_rank].set_linewidth(1.5)

    # Mark its position (useful when there are many clusters)
    ax.axvline(top_rank + 1, color="red", linestyle="--", alpha=0.7)

    ax.set_title(f"{col} (n={len(sizes)} clusters)")
    ax.set_xlabel("Clusters (ranked by size)")
    ax.set_ylabel("Number of hits")

    ax.text(
        0.98, 0.98,
        f"Top-hit cluster: {top_size} hits\nRank: {top_rank + 1}",
        transform=ax.transAxes,
        ha="right", va="top"
    )

axes.flat[-1].axis("off")

plt.tight_layout()
plt.show()
import numpy as np
import matplotlib.pyplot as plt

top_hit = cagE["psi1"][0].sort_values(
    ["evalue", "bitscore"],
    ascending=[True, False]
).iloc[0]["hit"]

seqs = cagE["psi1seqs"].df

# Automatically adjust grid size
ncols = 2 if len(cols) <= 4 else 3
nrows = (len(cols) + ncols - 1) // ncols

fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4 * nrows), squeeze=False)

for ax, col in zip(axes.flat, cols):

    # Cluster sizes, sorted largest to smallest
    sizes = seqs[col].value_counts()

    # Identify the cluster containing the best PSI-BLAST hit
    top_cluster = seqs.loc[seqs["id"] == top_hit, col].iloc[0]
    top_rank = sizes.index.get_loc(top_cluster)
    top_size = sizes.iloc[top_rank]

    # One bar per cluster
    bars = ax.bar(np.arange(1, len(sizes) + 1), sizes.values, width=1)

    # Highlight the cluster
    bars[top_rank].set_hatch("///")
    bars[top_rank].set_edgecolor("red")
    bars[top_rank].set_linewidth(1.5)

    ax.axvline(top_rank + 1, color="red", linestyle="--", alpha=0.7)

    ax.set_title(f"{col} (n={len(sizes)} clusters)")
    ax.set_xlabel("Clusters (ranked by size)")
    ax.set_ylabel("Number of sequences")

    ax.text(
        0.98, 0.98,
        f"Highlighted cluster size: {top_size} sequences\nRank: {top_rank + 1}",
        transform=ax.transAxes,
        ha="right", va="top"
    )

# Remove any unused panels
for ax in list(axes.flat)[len(cols):]:
    fig.delaxes(ax)

plt.tight_layout()
plt.show()
import numpy as np
import matplotlib.pyplot as plt

top_hit = cagE["psi1"][0].sort_values(
    ["evalue", "bitscore"],
    ascending=[True, False]
).iloc[0]["hit"]

seqs = cagE["psi1seqs"].df

ncols = 2 if len(cols) <= 4 else 3
nrows = (len(cols) + ncols - 1) // ncols

fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4 * nrows), squeeze=False)

for ax, col in zip(axes.flat, cols):
    # cluster sizes sorted from largest to smallest
    sizes = seqs[col].value_counts()

    # cluster containing the top PSI-BLAST hit
    top_cluster = seqs.loc[seqs["id"] == top_hit, col].iloc[0]
    top_rank = sizes.index.get_loc(top_cluster) + 1   # 1-based rank
    top_size = sizes.loc[top_cluster]

    ranks = np.arange(1, len(sizes) + 1)

    # rank-size curve
    ax.plot(ranks, sizes.values, linewidth=2)

    # highlight the cluster containing the top hit
    ax.scatter(top_rank, top_size, s=60, zorder=3)
    ax.axvline(top_rank, color="red", linestyle="--", alpha=0.7)

    ax.set_title(f"{col} (n={len(sizes)} clusters)")
    ax.set_xlabel("Clusters (ranked by size)")
    ax.set_ylabel("Number of sequences")

    ax.text(
        0.98, 0.98,
        f"Top-hit cluster rank: {top_rank}\nCluster size: {top_size}",
        transform=ax.transAxes,
        ha="right", va="top"
    )

    # optional: helps when the distribution is very skewed
    ax.set_yscale("log")

# remove unused panels
for ax in list(axes.flat)[len(cols):]:
    fig.delaxes(ax)

plt.tight_layout()
plt.show()
cagE['psi1seqs'].df.query('id == "UniRef50_Q48252"')[['c80i70','c80i50','c80i30','c0i0']]
cagE['psi1seqs'].df.c80i70.value_counts()
get_ipython().run_line_magic('page', '')
cagE['psi1'][0].head()
cagE['psi1seqs'].to_hist()
import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

cagE["psi1"][0].evalue.apply(
    lambda x: -np.log10(max(x, 1e-300))
).hist(ax=axes[0], bins=30)

cagE["psi1"][0].bitscore.hist(ax=axes[1], bins=30)
cagE["psi1"][0].length.hist(ax=axes[2], bins=30)

for ax, title in zip(axes, ["-log10(E-value)", "Bit score", "Alignment length"]):
    ax.set_title(title)
    ax.set_xlabel(title)
    ax.set_ylabel("Frequency")

fig.text(
    0.99, 0.99,
    f"Unique hits: {cagE['psi1'][0].hit.nunique()}",
    ha="right", va="top", fontsize=11
)

plt.tight_layout()
plt.show()
cagE['psi1'][0].head()
cagE['psi1seqs'].df.query('id == "UniRef50_Q48252"')
cagE['psi1'][0].sort_values(['iteration','evalue','bitscore'], ascending=[True,True,False]).head()
cagE['psi1seqs'].df.eval(f'a = id.isin({["UniRef50_Q48252"]})').groupby(['c80i30','c0i0']).agg(q=('a','any'), **{ x:(x,'nunique') for x in (['hit','c100i100','c80i70','c80i50']) })
cagE['psi1seqs'].df.eval(f'a = id.isin({["UniRef50_Q48252"]})').groupby(['c80i30','c0i0']).agg(q=('a','any'), **{ x:(x,'nunique') for x in (['id','c100i100','c80i70','c80i50']) })
cagE['psi1seqs'].df.eval(f'a = id.isin({["UniRef50_Q48252"]})').groupby(['c80i30','c0i0']).agg(q=('a','any'), **{ x:(x,'nunique') for x in (['id','c100i100','c80i70','c80i50']) }).sort_values(reversed(['id','c100i100','c80i70','c80i50']))
cagE['psi1seqs'].df.eval(f'a = id.isin({["UniRef50_Q48252"]})').groupby(['c80i30','c0i0']).agg(q=('a','any'), **{ x:(x,'nunique') for x in (['id','c100i100','c80i70','c80i50']) }).sort_values(list(reversed(['id','c100i100','c80i70','c80i50'])))
get_ipython().run_line_magic('page', '')
cagE['psi1seqs'].filter('c0i0 == "UniRef50_A0A6P2RT67"').align()._view_sequence()
cagE['psi1seqs'].df.eval(f'a = id.isin({["UniRef50_Q48252"]})').groupby(['c80i30','c0i0']).agg(q=('a','any'), **{ x:(x,'nunique') for x in (['id','c100i100','c80i70','c80i50']) }).sort_values(list(reversed(['id','c100i100','c80i70','c80i50'])))
get_ipython().run_line_magic('page', '')
cagE['psi1seqs'].df.nunique()
cagE['psi1seqs'].filter(f'id == c30i0').align()._view_sequence()
cagE['psi1seqs'].filter(f'id == c80i30').align()._view_sequence()
cagE['psi1seqs'].filter(f'id == c80i30').to_hist()
cagE['psi1seqs'].filter(f'id == c80i30', maxlen=2000, minlen=600).to_hist()
cagE['psi1seqs'].filter(f'id == c80i30', maxlen=1500, minlen=600).to_hist()
cagE['psi1seqs'].df['evalue'] = cagE['psi1seqs'].df.id.map(cagE['psi1'][0].set_index('hit').evalue)
cagE['psi1seqs'].df['evalue'] = cagE['psi1seqs'].df.id.map(cagE['psi1'][0].sort_values(['evalue']).drop_duplicates().set_index('hit').evalue)
cagE['psi1seqs'].df['evalue'] = cagE['psi1seqs'].df.id.map(cagE['psi1'][0].sort_values(['evalue']).drop_duplicates('hit').set_index('hit').evalue)
cagE['psi1seqs'].sort(['evalue'], inplace=True)
cagE['psi1seqs']._view_sequence()
cagE['psi1'][0]
get_ipython().run_line_magic('page', '')
cagE['psi1pfam'] = riu.architecture(rdbh.hmmscan(cagE['psi1seqs'], models_path="/scratch/global/databases/pfam/Pfam-A.hmm").filter('evalue <= 0.1'), **riu.config)
cagE['psi1pfam'] = riu.architecture(rdbh.hmmscan(cagE['psi1seqs'], models_path="/scratch/global/databases/pfam/Pfam-A.hmm").query('c_evalue <= 0.1'), **riu.config)
cagE['psi1pfam'] = riu.architecture(rdbh.hmmscan(cagE['psi1seqs'], models_path="/scratch/global/databases/pfam/Pfam-A.hmm").query('c_evalue <= 0.1'), **riu.config['hmmer'])
cagE['psi1pfam']
get_ipython().run_line_magic('page', '')
rdag
cagE['psi1hs'] = rdbh.hmmscan(cagE['psi1seqs'], models_path="/scratch/global/databases/pfam/Pfam-A.hmm")
cagE['psi1hs'].head()
cagE['psi1hs'].model.value_counts()
get_ipython().run_line_magic('page', '')
cagE['psi1'].hit.nunique()
cagE['psi1'][0].hit.nunique()
cagE['psi1hs'].model.value_counts()
get_ipython().run_line_magic('page', '')
get_ipython().run_line_magic('page', '')
cagE['psi1hsno'] = riu.filter_nonoverlapping_regions(cagE['psi1hs'], **riu.config['hmmer'])
cagE['psi1hsno'].model.value_counts()
get_ipython().run_line_magic('page', '')
cagE['psi1'][0].hit.nunique()
get_ipython().run_line_magic('page', '')
cagE['psi1hsno'].model.value_counts()
get_ipython().run_line_magic('page', '')
cagE['psi1hsno'].query('model != "P-loop_TraG"')
get_ipython().run_line_magic('page', '')
cagE['psi1hsno'].query('evalue <= 1 and model != "P-loop_TraG"')
get_ipython().run_line_magic('page', '')
cagE['psi1hsno'][~cagE['psi1hsno'].sequence.isin(cagE['psi1hsno'].query('model == "P-loop_TraG"').sequence)]
get_ipython().run_line_magic('page', '')
cagE['psi1hsno'][~cagE['psi1hsno'].sequence.isin(cagE['psi1hsno'].query('model == "P-loop_TraG"').sequence) & cagE['psi1hsno'].sequence.isin(cagE['psi1'].query('evalue <= 0.1').hit)]
cagE['psi1hsno'][(~cagE['psi1hsno'].sequence.isin(cagE['psi1hsno'].query('model == "P-loop_TraG"').sequence)) & cagE['psi1hsno'].sequence.isin(cagE['psi1'][0].query('evalue <= 0.1').hit)]
get_ipython().run_line_magic('page', '')
get_ipython().run_line_magic('page', '')
cagE['psi1hsno'][(~cagE['psi1hsno'].sequence.isin(cagE['psi1hsno'].query('model == "P-loop_TraG"').sequence)) & cagE['psi1hsno'].sequence.isin(cagE['psi1'][0].query('evalue <= 0.1').hit)].sequence.drop_duplicates()
cagE['psi1hsno'][(~cagE['psi1hsno'].sequence.isin(cagE['psi1hsno'].query('model == "P-loop_TraG"').sequence)) & cagE['psi1hsno'].sequence.isin(cagE['psi1'][0].query('evalue <= 0.1').hit)].sequence.drop_duplicates().nunique()
cagE['psi1seqs'].df.assign(pfp=cagE['psi1seqs'].df.id.isin(cagE['psi1hsno'][(~cagE['psi1hsno'].sequence.isin(cagE['psi1hsno'].query('model == "P-loop_TraG"').sequence)) & cagE['psi1hsno'].sequence.isin(cagE['psi1'][0].query('evalue <= 0.1').hit)].sequence.drop_duplicates().tolist())).groupby('c0i0').agg(pfp=('pfp','any'), **{ x:(x,'nunique') for x in ['id','c100i100','c80i70','c80i50','c80i30'] })
cagE['psi1seqs'].df.assign(pfp=cagE['psi1seqs'].df.id.isin(cagE['psi1hsno'][(~cagE['psi1hsno'].sequence.isin(cagE['psi1hsno'].query('model == "P-loop_TraG"').sequence)) & cagE['psi1hsno'].sequence.isin(cagE['psi1'][0].query('evalue <= 0.1').hit)].sequence.drop_duplicates().tolist())).groupby('c0i0').agg(pfp=('pfp','any'), **{ x:(x,'nunique') for x in ['id','c100i100','c80i70','c80i50','c80i30'] }).sort_values(list(reversed(['id','c100i100','c80i70','c80i50','c80i30'])))
cagE['psi1'][0].head()
cagE['psi1seqs'].df.query('id == "UniRef50_Q48252"').c80i0
cagE['psi1seqs'].df.query('id == "UniRef50_Q48252"').c0i0
cagE['psi1seqs'].df.query('id == "UniRef50_Q48252"').c80i30
cagE['psi1seqs'].df.query('id == "UniRef50_J0J980"').c80i30
cagE['psi1seqs'].df.query('id == "UniRef50_J0J980"').c0i0
cagE['psi1seqs'].df.query('id == "UniRef50_J0J980"').c80i30
cagE['psi1seqs'].df.assign(pfp=cagE['psi1seqs'].df.id.isin(cagE['psi1hsno'][(~cagE['psi1hsno'].sequence.isin(cagE['psi1hsno'].query('model == "P-loop_TraG"').sequence)) & cagE['psi1hsno'].sequence.isin(cagE['psi1'][0].query('evalue <= 0.1').hit)].sequence.drop_duplicates().tolist())).groupby('c0i0').agg(pfp=('pfp','any'), **{ x:(x,'nunique') for x in ['id','c100i100','c80i70','c80i50','c80i30'] }).sort_values(list(reversed(['id','c100i100','c80i70','c80i50','c80i30'])))
cagE['psi1seqs'].df.assign(pfp=cagE['psi1seqs'].df.id.isin(cagE['psi1hsno'][(~cagE['psi1hsno'].sequence.isin(cagE['psi1hsno'].query('model == "P-loop_TraG"').sequence)) & cagE['psi1hsno'].sequence.isin(cagE['psi1'][0].query('evalue <= 0.1').hit)].sequence.drop_duplicates().tolist())).groupby('c0i0').agg(pfp=('pfp','sum'), **{ x:(x,'nunique') for x in ['id','c100i100','c80i70','c80i50','c80i30'] }).sort_values(list(reversed(['id','c100i100','c80i70','c80i50','c80i30'])))
cagE['psi1seqs'].filter('c0i0 == "UniRef50_A0A2A2K1R1" and id == c80i50').align()._view_sequence()
cagE['psi1seqs'].df.assign(pfp=cagE['psi1seqs'].df.id.isin(cagE['psi1hsno'][(~cagE['psi1hsno'].sequence.isin(cagE['psi1hsno'].query('model == "P-loop_TraG"').sequence)) & cagE['psi1hsno'].sequence.isin(cagE['psi1'][0].query('evalue <= 0.1').hit)].sequence.drop_duplicates().tolist())).groupby('c0i0').agg(pfp=('pfp','sum'), **{ x:(x,'nunique') for x in ['id','c100i100','c80i70','c80i50','c80i30'] }).sort_values(list(reversed(['id','c100i100','c80i70','c80i50','c80i30'])))
cagE['psi1seqs'].filter('c0i0 == "UniRef50_A0A8T3YM87" and id == c80i50').align()._view_sequence()
cagE['psi1seqs'].df.assign(pfp=cagE['psi1seqs'].df.id.isin(cagE['psi1hsno'][(~cagE['psi1hsno'].sequence.isin(cagE['psi1hsno'].query('model == "P-loop_TraG"').sequence)) & cagE['psi1hsno'].sequence.isin(cagE['psi1'][0].query('evalue <= 0.1').hit)].sequence.drop_duplicates().tolist())).groupby('c0i0').agg(pfp=('pfp','sum'), **{ x:(x,'nunique') for x in ['id','c100i100','c80i70','c80i50','c80i30'] }).sort_values(list(reversed(['id','c100i100','c80i70','c80i50','c80i30'])))
gnc = ncbi.GeneNeighborhoodCursor(mirror="/scratch/global/databases/genomes")
mc
mc = uniprot.MappingCursor()
mc["UniRef50_A0A0R3Q3Q1"]
mc = uniprot.MappingCursor(host="davinci")
mc["UniRef50_A0A0R3Q3Q1"]
cagE['psi1seqs'].df.assign(pfp=cagE['psi1seqs'].df.id.isin(cagE['psi1hsno'][(~cagE['psi1hsno'].sequence.isin(cagE['psi1hsno'].query('model == "P-loop_TraG"').sequence)) & cagE['psi1hsno'].sequence.isin(cagE['psi1'][0].query('evalue <= 0.1').hit)].sequence.drop_duplicates().tolist())).groupby('c0i0').agg(pfp=('pfp','sum'), **{ x:(x,'nunique') for x in ['id','c100i100','c80i70','c80i50','c80i30'] }).sort_values(list(reversed(['id','c100i100','c80i70','c80i50','c80i30'])))
mc["T1BDM8"]
cagE['psi1'][0].head()
mc["Q48252"]
get_ipython().run_line_magic('pinfo2', 'mc.fetchall')
cagE['psi1u2n'] = mc.fetchall(cagE['psi1'][0].hit.drop_duplicates().tolist(), source="UniRef50", target=["RefSeq","EMBL-CDS"])
cagE['psi1u2n'].head()
cagE['psi1u2n']
mc.missing
get_ipython().run_line_magic('page', '')
mc["UPI00313314C7"]
mc.cursors['webapi["UPI00313314C7"]
mc.cursors['webapi']["UPI00313314C7"]
mc.cursors['webapi'].fetchall(["UPI00313314C7"], source="UniParc")
mc.cursors['webapi'].fetchall(["UPI00313314C7"], source="Proteomes")
mc.cursors['webapi'].fetchall(["UPI00313314C7"], source="Uniprot-KB")
mc.cursors['webapi'].fetchall(["UPI00313314C7"], source="UniParc")
mc.cursors['webapi'].fetchall(["UPI00313314C7"], source="UniParc", target="RefSeq")
fc
fc = uniprot.FastaCursor(local_fasta_path="/scratch/global/databases/fadb/uniprot/uniref50")
fc["UniRef50_UPI00313314C7"]
mc.cursors['clickhouse'].client
help(mc.cursors['clickhouse'].client)
mc.cursors['clickhouse'].client.query_df("SELECT * FROM idmapping LIMIT 3")
mc.cursors['clickhouse'].client.query_df("SELECT * FROM uniprot.idmapping LIMIT 3")
mc.cursors['clickhouse'].client.query_df("SELECT * FROM rotifer.idmapping LIMIT 3")
mc.cursors['clickhouse'].client.query_df("SELECT * FROM rotifer.idmapping WHERE id == 'UPI00313314C7'")
mc.cursors['clickhouse'].client.query_df("SELECT * FROM rotifer.idmapping WHERE id == 'UniRef50_UPI00313314C7'")
mc.cursors['clickhouse'].client.query_df("SELECT * FROM rotifer.idmapping WHERE id = 'UniRef50_UPI00313314C7'")
mc.cursors['clickhouse'].client.query_df("SELECT * FROM rotifer.idmapping WHERE id like '%UPI00313314C7%'")
mc.cursors['clickhouse'].client.query_df("SELECT * FROM rotifer.idmapping WHERE id like '%UPI0%'")
mc.cursors['clickhouse'].client.query_df("SELECT * FROM rotifer.idmapping LIMIT 3")
mc.cursors['clickhouse'].client.query_df("SELECT * FROM rotifer.idmapping WHERE id_type == 'UniParc'")
mc.cursors['clickhouse'].client.query_df("SELECT * FROM rotifer.idmapping WHERE id_type = 'UniParc'")
rdbs.sequence(["WP_000496009.1"], local_database_path=os.path.join(rotifer.config["data"],"fadb","nr","nr")))
rdbs.sequence(["WP_000496009.1"], local_database_path=os.path.join(rotifer.config["data"],"fadb","nr","nr"))
