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
get_ipython().run_line_magic('who', '')
cagE['ATPaseenv']._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[61] == "-"')._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[61] == "-"').slice((60,1490))._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[61] == "-"').slice((60,1490)).align(method='lisni')._view_sequence()
get_ipython().run_line_magic('pwd', '')
cagE['ATPaseenv'].filter('sequence.str[61] == "-"').slice((60,1490)).align(method='linsi')._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[61] == "-"').slice((60,1490)).align(method='linsi').to_hist()
cagE['ATPaseenv'].filter('sequence.str[61] != "-"').slice((55,1490)).align(method='linsi').to_hist()
cagE['ATPaseenv'].to_hist()
cagE['ATPaseenv'].df.columns
cagE['ATPaseenv'].df.id.value_counts().head()
cagE['ATPaseenv'].df.pid.value_counts().max()
cagE['ATPaseenv'].df.nunique()
cagE['ATPaseenv'] = cagE['ATPaseenv'].sort(['evalue']) 
cagE['ATPaseenv']._view_sequence()
cagE['ATPaseenv'].df.groupby('c0i0').c80i50.nunique()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698"').align()._view_sequence()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').align()._view_sequence()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').to_hist()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50')._view_sequence()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').columns
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.columns
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.id
mc
mc = uniprot.MappingCursor()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.drop_duplicates().tolist()
get_ipython().run_line_magic('who', '')
cagE['ATPasemc'] = mc.fetchall(cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.drop_duplicates().tolist(), source="UniRef50", target="TaxId")
mc = uniprot.MappingCursor(host="davinci")
cagE['ATPasemc'] = mc.fetchall(cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.drop_duplicates().tolist(), source="Uniref50", target="TaxId")
cagE['ATPasemc'] = mc.fetchall(cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.drop_duplicates().tolist(), source="UniRef50", target="TaxId")
cagE['ATPasemc'] = mc.fetchall(cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.drop_duplicates().tolist(), source="UniRef50", target="TaxID")
cagE['ATPasemc'] = mc.fetchall(cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.drop_duplicates().tolist(), source="UniRef50", target="Taxid")
mc.fetchall(cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.drop_duplicates().sample(5).tolist(), source="UniRef50")
from rotifer.pandas import functions as rpf
rpf.print_everything()
get_ipython().run_line_magic('page', '_')
mc.fetchall(cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.drop_duplicates().sample(5).tolist(), source="UniRef50", target="NCBI_TaxID")
tc = ncbi.TaxonomyCursor()
tc.fetchall(mc.fetchall(cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.drop_duplicates().sample(5).tolist(), source="UniRef50", target="NCBI_TaxID").target.drop_duplicates().tolist()))
tc.fetchall(mc.fetchall(cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.drop_duplicates().sample(5).tolist(), source="UniRef50", target="NCBI_TaxID").target.drop_duplicates().tolist())
cagE["ATPaseU2T"] = mc.fetchall(cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.drop_duplicates().sample(5).tolist(), source="UniRef50", target="NCBI_TaxID")
cagE["ATPaseU2T"] = mc.fetchall(cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.drop_duplicates().tolist(), source="UniRef50", target="NCBI_TaxID")
cagE["ATPaseU2T"]
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A3B0S2C7/444-698" and id == c80i50').df.pid.isin(cagE["ATPaseU2T"].source).value_counts()
cagE['ATPaseenv'].df.description.head()
cagE["ATPaseU2T"].target.nunique()
cagE["ATPaseenv"].df.description.head()
cagE["ATPaseenv"].df.description.str.extract("TaxID=([0-9]+)").head()
cagE["ATPaseenv"].df.description.str.extract("TaxID=([0-9]+)")
cagE["ATPaseenv"].df['taxid'] = cagE["ATPaseenv"].df.description.str.extract("TaxID=([0-9]+)")
cagE["ATPaseenv"].df['taxid'].value_counts()
cagE["ATPaseenv"].df['taxid'].value_counts(dropna=False)
get_ipython().run_line_magic('page', '')
cagE["ATPaseenv"].df.taxid.isna().value_counts()
(cagE["ATPaseenv"].df.taxid == "").value_counts()
cagE["ATPaseenv"].df.taxid.nunique()
del cagE["ATPaseU2T"]
cagE["ATPaseTax"] = tc.fetchall(cagE["ATPaseenv"].df.taxid.drop_duplicates().tolist())
cagE["ATPaseTax"]
get_ipython().run_line_magic('page', '')
cagE["ATPaseTax"].lineage.value_counts()
get_ipython().run_line_magic('page', '')
cagE["ATPaseTax"].classification.head()
cagE["ATPaseTax"].classification.extract("(\S+ota)\;")
cagE["ATPaseTax"].classification.str.extract("(\S+ota)\;")
cagE["ATPaseTax"].classification.str.extract("(\S+ati)\;")
cagE["ATPaseTax"].classification.head()
cagE["ATPaseTax"].classification.str.split("; ", expand=True)[0].value_counts()
cagE["ATPaseTax"]["superkingdom"] = cagE["ATPaseTax"].classification.str.split("; ", expand=True)[0]
cagE["ATPaseTax"]
get_ipython().run_line_magic('page', '')
cagE["ATPaseTax"] = tc.fetchall(cagE["ATPaseenv"].df.taxid.drop_duplicates().tolist())
cagE["ATPaseTax"]
get_ipython().run_line_magic('page', '')
cagE["ATPaseTax"].superkingdom.value_counts()
cagE["ATPaseTax"]["kingdom"] = cagE["ATPaseTax"].classification.str.extract("(\S+ati)\;")
cagE["ATPaseTax"]["kingdom"] = cagE["ATPaseTax"].classification.str.extract(r"(\S+ati)\;")
cagE["ATPaseTax"]["kingdom"].value_counts()
cagE["ATPaseTax"]["phylum"] = cagE["ATPaseTax"].classification.str.extract(r"(\S+ota)\;")
cagE["ATPaseTax"]["phylum"].value_counts()
cagE["ATPaseTax"]
get_ipython().run_line_magic('page', '')
cagE["ATPaseTax"].columns
cagE["ATPaseenv"].columns
cagE["ATPaseenv"].df.columns
cagE["ATPaseTax"].taxid.value_counts().max()
cagE["ATPaseenv"] = cagE["ATPaseenv"].df.merge(cagE["ATPaseTax"], on="taxid", how="left")
cagE["ATPaseenv"].df.isna()
z = rdbs.sequence()
z.df
z.df = cagE["ATPaseenv"].df
z.df = cagE["ATPaseenv"]
z.sample(3)._view_sequnece()
z.sample(3)._view_sequence()
cagE["ATPaseenv"] = z.copy()
del z
cagE["ATPaseenv"]._view_sequence()
cagE["ATPaseenv"].df.groupby(['c0i0','kingdom']).pid.nunique()
cagE["ATPaseenv"]._view_sequence()
cagE["ATPaseenv"].to_hist()
cagE["ATPaseenv"]._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[61] != "-"').slice((55,1490)).align(method='linsi')._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[61] != "-"').slice((40,1490)).sort(['c0i0','c80i30']).align()._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[61] != "-"').slice((40,1490)).sort(['c80i50','c80i70']).align()._view_groups(groupby='c80i30')
# We can inspect the query by either running:
#cagE['data'].view()
# or cagE['data']._view_sequence() on a terminal window
len(cagE["Atpasepsi1"][0])
cagE["Atpasepsi1"][0].hit.count()# We can check how many 
1350 - 1028
cagE['ATPaseseqs'].to_hist().filter(maxlen = 1500)
cagE['ATPaseseqs']._view_sequence()
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseseqs'].to_hist() # A few number of sequences are making our histogram hard to read, let's remove them
cagE['ATPaseseqs'].to_hist().filter(maxlen = 1500)
cagE['ATPaseseqs'].filter(maxlen = 1500).to_hist()
cagE['Atpasepsi1'][0].querystart.min()
cagE['Atpasepsi1'][0].queryend.min()
cagE['Atpasepsi1'][0].querstart.max()
cagE['Atpasepsi1'][0].querystart.max()
cagE['Atpasepsi1'][0].querystart.hist()
cagE['Atpasepsi1'][0].querystart.max()
cagE['ATPaseenv'].to_hist()
cagE['ATPaseenv'].df.c0i0.value_counts()
get_ipython().run_line_magic('who', '')
get_ipython().run_line_magic('pinfo2', 'rdag.extend_aln')
