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
get_ipython().run_line_magic('pinfo2', 'rdagf.add_cordinates_to_aln')
rdag.add_cordinates_to_aln
get_ipython().run_line_magic('pinfo2', 'rdag.add_cordinates_to_aln')
cagE['ATPaseenv'].df.columns
cagE['ATPaseenv'].df[['id','start','end']].head(15)
cagE['ATPaseenv'].df.assign(
    current_length=lambda d: d.sequence.str.replace('-', '', regex=False).str.len(),
    coordinate_length=lambda d: d.end - d.start + 1
)[['id', 'length', 'current_length', 'coordinate_length']].head(15)
get_ipython().run_line_magic('pinfo', 'rdag.add_cordinates_to_aln')
cagE['ATPaseenv'].df['bckp'] = cagE['ATPaseenv'].df.id
cagE['ATPaseenv'].df.id = cagE['ATPaseenv'].df.pid
cagE['ATPaseenv'].df.colunns
cagE['ATPaseenv'].df.columns
meta = cagE['ATPaseenv'].df[['description', 'query', 'hitstart','hitend', 'evalue', 'querycoverage', 'querystart', 'queryend','iteration', 'bitscore', 'og_length', 'pid', 'c100i100','c80i70', 'c80i50', 'c80i30', 'c0i0', 'taxid', 'organism','superkingdom', 'lineage', 'classification', 'alternative_taxids','kingdom', 'phylum', 'bckp']]
cagE['ATPaseenv'] = rdag.extend_aln(cagE['ATPaseenv'],n_terminal=50,c_terminal=50,full_length_seqobj=cagE['ATPaseseqs'])
cagE['ATPaseenv']._view_sequence()
cagE['ATPaseenv'].to_hist()
cagE['ATPaseenv'].sort(['length'])._view_sequence()
cagE['ATPaseenv'].sort(['length'],ascending=[False])._view_sequence()
cagE['ATPaseenv'].align().to_hist()
cagE['ATPaseenv'].align(method='linsi').to_hist()
cagE['ATPaseenv'].slice((1,-1)).align().to_hist()
get_ipython().run_line_magic('pinfo2', 'rdbs.sequence.slice')
cagE['ATPaseenv']._reset()
cagE['ATPaseenv'].to_hist()
# Let's check if something drastic changed
cagE['ATPaseenv'].to_hist()
cagE['ATPaseenv']._view_sequence()
cagE['ATPaseenv'].sort([241])._view_sequence()
cagE['ATPaseenv'].sort([1548])._view_sequence()
cagE['ATPaseenv'].sort([1538])._view_sequence()
cagE['ATPaseenv'].df.query('sequence.str[240] == "-"').count() # Notice the count is 0 based
cagE['ATPaseenv'].df.query('sequence.str[240] == "-"')['id'].count() # Notice the count is 0 based
cagE['ATPaseenv'].df.query('sequence.str[1537] == "-"')['id'].count()
cagE['ATPaseenv'].df.query('sequence.str[1537] == "-"')['sequence'].head(1)
cagE['ATPaseenv'].df.query('sequence.str[1537] != "-"')['sequence'].head(1)
cagE['ATPaseenv'].df.query('sequence.str[1537] != "-"')['sequence'].str[1537].head(1)
cagE['ATPaseenv'].df.query('sequence.str[240] == "-"')['sequence'].str[240].head(1)#['id'].count() # Notice the count is 0 based
cagE['ATPaseenv'].df.query('sequence.str[240] != "-"')['sequence'].str[240].head(1)#['id'].count() # Notice the count is 0 based
cagE['ATPaseenv'].df.query('sequence.str[1537] == "-"')['id'].count()
cagE['ATPaseenv'].df.query('sequence.str[240] == "-"')['id'].count() # Notice the count is 0 based
cagE['ATPaseenv'].df.query('sequence.str[240] == "-" and sequence.str[1537] == "-"')['id'].count()
cagE['ATPaseenv']._view_sequence()
cagE['ATPaseenv'].sort(['list'],id_list=['UniRef50_Q48252'])_view_sequence()
cagE['ATPaseenv'].sort(['list'],id_list=['UniRef50_Q48252'])._view_sequence()
cagE['ATPaseenv'].sort(['list'],id_list=['UniRef50_Q48252']).slice((210,1697))._view_sequence()
cagE['ATPaseenv'].sort(['list'],id_list=['UniRef50_Q48252']).slice((220,1687))._view_sequence()
cagE['ATPaseenv'].sort(['list'],id_list=['UniRef50_Q48252']).slice((260,1685))._view_sequence()
cagE['ATPaseenv'].sort(['list'],id_list=['UniRef50_Q48252']).slice((226,1685))._view_sequence()
cagE['ATPaseenv'].sort(['list'],id_list=['UniRef50_Q48252']).slice((226,1685)).align()._view_sequence()
cagE['ATPaseenv'] = cagE['ATPaseenv'].sort(['list'],id_list=['UniRef50_Q48252']).slice((226,1685)).align() # trimm down to the atpase domain, we need to try to get as close as possible
cagE['ATPaseenv'].sort(['list'],id_list=['UniRef50_Q48252')._view_sequence()
cagE['ATPaseenv'].sort(['list'],id_list=['UniRef50_Q48252'])._view_sequence()
cagE['ATPaseenv'].sort(['list',20],id_list=['UniRef50_Q48252'])._view_sequence()
cagE['ATPaseenv'].sort(['list',20],id_list=['UniRef50_Q48252'])..df.query('sequence.str[20]')['id'].count()
cagE['ATPaseenv'].sort(['list',20],id_list=['UniRef50_Q48252']).df.query('sequence.str[20]')['id'].count()
cagE['ATPaseenv'].df.query('sequence.str[20]')['id'].count()
cagE['ATPaseenv'].df.query('sequence.str[20]')['id'].count()
cagE['ATPaseenv'].df.query('sequence.str[20]== "-"')['id'].count()
cagE['ATPaseenv'].df.query('sequence.str[20]== "-"')['id'].count()
cagE['ATPaseenv'].sort(['list',20],id_list=['UniRef50_Q48252'])._view_sequence()
cagE['ATPaseenv'].df.query('sequence.str[16]== "-"')['id'].count()
cagE['ATPaseenv'].df.query('sequence.str[16]!= "G"')['id'].count()
cagE['ATPaseenv'].df.query('sequence.str[16]!= "G"')['sequence'].str[16].head(12)
cagE['ATPaseenv'].df.query('sequence.str[16] != "-"')['sequence'].str[16].head(12)
cagE['ATPaseenv'].df.query('sequence.str[16] != "-"')['sequence'].str[16:20].head(12)
cagE['ATPaseenv'].df.query('sequence.str[16] != "-"')['sequence'].str[13:26].head(12)
cagE['ATPaseenv'].df.query('sequence.str[16] != "-"')['sequence'].str[15:26].head(12)
cagE['ATPaseenv'].df.query('sequence.str[15]!= "G"')['id'].count()
cagE['ATPaseenv'].df.query('sequence.str[15]== "-"')['id'].count()
cagE['ATPaseenv'].df.query('sequence.str[15]== "-"')['id'].count()
cagE['ATPaseenv'].add_cluster(cascade=[(100,100),(80,70),(80,50),(80,30),(0,0)], inplace=True)
cagE['ATPaseenv'].df.nunique()
cagE['ATPaseenv'].to_hist()
cagE['ATPaseenv'].df.c0i0.value_counts()
#this is a temporary variable we can delete after 
ungap = cagE['ATPaseenv'].copy()
ungap.df['sequence'] = seqs.df.sequence.str.replace('-', '', regex=False)

# Note to self - add ungapping step on hmmscan functions
cagE['ATPaseenv']
#this is a temporary variable we can delete after 
ungap = cagE['ATPaseenv'].copy()
ungap.df['sequence'] = ungap.df.sequence.str.replace('-', '', regex=False)

# Note to self - add ungapping step on hmmscan functions
cageE
cagE['Atpasehs'] = rdbh.hmmscan(ungap, models_path="/scratch/global/databases/pfam/Pfam-A.hmm")
cagE['Atpasepfam'] = riu.architecture(cagE['Atpasehs'].query('c_evalue <= 0.1'), **riu.config['hmmer'])
cagE['Atpasepfam'].value_counts().head(10)
cagE['Atpasepfam'].architecture.value_counts().head(10)
# First, let's extract the taxonomic IDs from the description

cagE["ATPaseenv"].df['taxid'] = cagE["ATPaseenv"].df.description.str.extract("TaxID=([0-9]+)")
cagE["ATPaseTax"].head(5)
# Let's add the taxonomy to our sequence object
cagE["ATPaseenv"].df = cagE["ATPaseenv"].df.merge(cagE["ATPaseTax"], on="taxid", how="left")
# Exporting
cagE["ATPaseenv"].to_file('Atpasepsi1aln.fasta')
cagE['ATPasepfam'].head(2)
cagE['Atpasepfam'].head(2)
cagE['ATPaseenv'].df['pfam'] = cagE['ATPaseenv'].df.id.map(cagE['Atpasepfam'].set_index('sequence').architecture)
cagE['ATPaseenv']._view_sequence()
cagE['ATPaseenv']._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str.[1321] != "-"').slice((8,1456)).sort(['c80i50','c80i70']).align()._view_groups(groupby='c80i30')
cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').slice((8,1456)).sort(['c80i50','c80i70']).align()._view_groups(groupby='c80i30')
cagE["ATPaseenv"].df.superkingdom.value_counts()
cagE["ATPaseenv"].df.phylum.value_counts()
cagE["ATPaseenv"].df.kingdom.value_counts()
cagE["ATPaseenv"].df.c80i30.nunique()
cagE["ATPaseenv"].df.c80i50.nunique()
cagE["ATPaseenv"].df.groupby('c80i50')['phylum'].value_counts()
from rotifer.pandas import functions as rpf
rpf.print_everything()
cagE["ATPaseenv"].df.groupby('c80i50')['phylum'].value_counts()
get_ipython().run_line_magic('page', '')
cagE["ATPaseenv"].df.groupby(['c0i0','kingdom']).id.nunique()
cagE['ATPaseenv'].df.groupby('c0i0').c80i50.nunique()
cagE["ATPaseenv"].df.groupby('phylun')['c80i50'].value_counts()
cagE["ATPaseenv"].df.groupby('phylum')['c80i50'].value_counts()
cagE["ATPaseenv"].df.groupby('c80i50')['phylum'].value_counts().head(15)
cagE["ATPaseenv"].df.groupby(['c0i0','kingdom']).id.nunique()
cagE["ATPaseenv"].df.groupby(['c80i50','kingdom']).id.nunique()
cagE["ATPaseenv"].df.groupby(['c80i50','phylum']).id.nunique()
cagE["ATPaseenv"].df.groupby(['c0i0','phylum']).id.nunique()
cagE['ATPaseenv'].df[['id','superkingdom','kindgdom','phylum','c80i70','c80i50','c80i30','c0i0']].to_csv('meta)fortree.tsv',sep="\t",index=False)
cagE['ATPaseenv'].df[['id','superkingdom','kingdom','phylum','c80i70','c80i50','c80i30','c0i0']].to_csv('meta)fortree.tsv',sep="\t",index=False)
cagE['ATPaseenv'].df[['id','superkingdom','kingdom','phylum','c80i70','c80i50','c80i30','c0i0']].to_csv('meta_fortree.tsv',sep="\t",index=False)
cagE['ATPaseenv'].df.c80i50.value_counts().head(7)
cagE['ATPaseenv'].df.c80i70.value_counts().head(7)
cagE['ATPaseenv'].df.c80i30.value_counts().head(7)
cagE['ATPaseenv'].df.c80i30.value_counts().head(7).sum()
cagE['ATPaseenv'].df['motif_ok'] = (cagE['ATPaseenv'].df.sequence.str[20] != '-') & (cagE['ATPaseenv'].df.sequence.str[1321] != '-')
cagE['ATPaseenv'].df['motif_ok'].value_counts()
cagE['ATPaseenv'].df.groupby('motif_ok').c0i0.value_counts()
cagE['ATPaseenv'].df.groupby('motif_ok').c0i0.value_counts(normalize=True)
cagE['ATPaseenv'].df.groupby('c0i0').motif_ok.valuue_counts(normalize=True)
cagE['ATPaseenv'].df.groupby('c0i0').motif_ok.value_counts(normalize=True)
cagE['ATPaseenv'].df.groupby('c80i30').motif_ok.value_counts(normalize=True)
cagE['ATPaseenv'].df.groupby('c80i30').motif_ok.value_counts()
cagE['ATPaseenv'].df['evalue'] = cagE['ATPaseenv'].df.id.map(cagE['Atpasepsi1'].drop_duplicates('hit').set_index('hit').evalue)
cagE['ATPaseenv'].df['evalue'] = cagE['ATPaseenv'].df.id.map(cagE['Atpasepsi1'][0].drop_duplicates('hit').set_index('hit').evalue)
cagE['ATPaseenv'].df[['id','evalue','motif_ok','superkingdom','kingdom','phylum','c80i70','c80i50','c80i30','c0i0']].to_csv('meta_fortree.tsv',sep="\t",index=False)
cagE['ATPaseenv'].df.groupby('phylum').motif_ok.value_counts(normalize=True)
cagE['ATPaseenv'].df.groupby('phylum').motif_ok.value_counts()
for col in ['c80i70', 'c80i50', 'c80i30']:
    print(f"\n{col}")
    print(atp.df[col].value_counts().describe())
for col in ['c80i70', 'c80i50', 'c80i30']:
    print(f"\n{col}")
    print(cagE['ATPaseenv'].df[col].value_counts().describe())
for col in ['c0i0','c80i70', 'c80i50', 'c80i30']:
    print(f"\n{col}")
    print(cagE['ATPaseenv'].df[col].value_counts().describe())
cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').slice((8,1456)).filter('id == "c80i30"').sort(['c80i50','c80i70']).align()._view_sequences()
cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').slice((8,1456)).filter('id == "c80i30"').sort(['c80i50','c80i70']).align()._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').slice((8,1456)).filter('id == "c80i30"').align()._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').slice((8,1456)).filter('id == "c80i30"').align()._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"')._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').filter('id == "c80i30"').slice((8,1456)).align()._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').filter('id == c80i30').slice((8,1456)).align()._view_sequence()
cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').filter('id == c80i50').slice((8,1456)).align()._view_sequence()
# We are going to check if taking one of each of c80i50 after removing the sequences that are kills our diversity
cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').filter('id == c80i50').slice((8,1456)).align().df[[col]].value_counts().describe())
# We are going to check if taking one of each of c80i50 after removing the sequences that are kills our diversity
cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').filter('id == c80i50').slice((8,1456)).align().df[[col]].value_counts().describe()
cols = ['c80i70', 'c80i50', 'c80i30', 'phylum']

cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').filter('id == c80i50').slice((8,1456)).align().df[cols].apply(lambda s: s.value_counts().describe())
cols = ['c80i70', 'c80i50', 'c80i30','c0i0', 'phylum']

cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').filter('id == c80i50').slice((8,1456)).align().df[cols].apply(lambda s: s.value_counts().describe())
pd.DataFrame({
        'before': cagE['ATPaseenv'].df[cols].nunique(),
        'after': (
            cagE['ATPaseenv']
            .filter('sequence.str[20] != "-" and sequence.str[1321] != "-"')
            .filter('id == c80i50')
            .slice((8,1456))
            .align()
            .df[cols]
            .nunique()
        )
    })
    .assign(
        lost=lambda d: d['before'] - d['after'],
        retained_pct=lambda d: 100 * d['after'] / d['before']
    )
    .round(2)
)
(
    pd.DataFrame({
        'before': cagE['ATPaseenv'].df[cols].nunique(),
        'after': (
            cagE['ATPaseenv']
            .filter('sequence.str[20] != "-" and sequence.str[1321] != "-"')
            .filter('id == c80i50')
            .slice((8,1456))
            .align()
            .df[cols]
            .nunique()
        )
    })
    .assign(
        lost=lambda d: d['before'] - d['after'],
        retained_pct=lambda d: 100 * d['after'] / d['before']
    )
    .round(2)
)
import pandas as pd
(
    pd.DataFrame({
        'before': cagE['ATPaseenv'].df[cols].nunique(),
        'after': (
            cagE['ATPaseenv']
            .filter('sequence.str[20] != "-" and sequence.str[1321] != "-"')
            .filter('id == c80i50')
            .slice((8,1456))
            .align()
            .df[cols]
            .nunique()
        )
    })
    .assign(
        lost=lambda d: d['before'] - d['after'],
        retained_pct=lambda d: 100 * d['after'] / d['before']
    )
    .round(2)
)
import pandas as pd
(
    pd.DataFrame({
        'before': cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').df[cols].nunique(),
        'after': (
            cagE['ATPaseenv']
            .filter('sequence.str[20] != "-" and sequence.str[1321] != "-"')
            .filter('id == c80i50')
            .slice((8,1456))
            .align()
            .df[cols]
            .nunique()
        )
    })
    .assign(
        lost=lambda d: d['before'] - d['after'],
        retained_pct=lambda d: 100 * d['after'] / d['before']
    )
    .round(2)
)
import pandas as pd
(
    pd.DataFrame({
        'before': cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').df[cols].nunique(),
        'after': (
            cagE['ATPaseenv']
            .filter('sequence.str[20] != "-" and sequence.str[1321] != "-"')
            .filter('id == c80i50')
            .slice((8,1456))
            .df[cols]
            .nunique()
        )
    })
    .assign(
        lost=lambda d: d['before'] - d['after'],
        retained_pct=lambda d: 100 * d['after'] / d['before']
    )
    .round(2)
)
cagE['it1atpsel'] = cagE['ATPaseenv'].filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').filter('id == c80i50').slice((8,1456)).align()
cagE['ATPaseenv']
cagE['ATPaseenv'].filter('id != cagE["it1atpsel"].id.tolist()')._view_sequence()
cagE['ATPaseenv'].filter( f'id not in {cagE["it1atpsel"].id.tolist()}')._view_sequence()
cagE['ATPaseenv'].filter( f'id not in {cagE["it1atpsel"].df.id.tolist()}')._view_sequence()
cagE['ATPaseenv'].filter( f'id not in {cagE["it1atpsel"].df.id.tolist()}')sort([20])._view_sequence()
cagE['ATPaseenv'].filter( f'id not in {cagE["it1atpsel"].df.id.tolist()}').sort([20])._view_sequence()
cagE['ATPaseenv'].df.sequence.str[15:20].head(3)
cagE['ATPaseenv'].df.sequence.str[15:25].head(3)
cagE['ATPaseenv'].df.sequence.str[10:30].head(3)
cagE['ATPaseenv'].df.sequence.str[14:30].head(3)
cagE['ATPaseenv'].df.sequence.str[20:30].head(3)
cagE['ATPaseenv'].filter( f'id not in {cagE["it1atpsel"].df.id.tolist()}').sort([20])._view_sequence()
cagE['ATPaseenv'].filter( f'id not in {cagE["it1atpsel"].df.id.tolist()}').sort([21])._view_sequence()
cagE['ATPaseenv'].filter(f'c80i50 not in {cagE["it1atpsel"].df.c80i50.dropna().unique().tolist()}').filter('sequence.str[20] != "-" and sequence.str[1321] != "-"')._view_sequence()
# Let's recover what's lost and check it
temp = cagE['ATPaseenv'].filter(f'c80i50 not in {cagE["it1atpsel"].df.c80i50.dropna().unique().tolist()}').filter('sequence.str[20] != "-" and sequence.str[1321] != "-"').sort(['evalue'])

temp = temp.filter(f'id in {temp.df.groupby("c80i50", sort=False).head(1).id.tolist()}')
temp._view_sequence()
cagE['ATPaseenv'].df.query('c0i50 == "UniRef50_A0A485BXW6"')
cagE['ATPaseenv'].df.query('c80i50 == "UniRef50_A0A485BXW6"')
cagE['ATPaseenv'].df.query('c80i50 == "UniRef50_A0A485BXW6"').id.nunique()
cagE['ATPaseenv'].df.query('c80i50 == "UniRef50_A0A485BXW6"')._view_sequence()
cagE['ATPaseenv'].filter('c80i50 == "UniRef50_A0A485BXW6"')._view_sequence()
temp._view_sequence()
cagE['ATPaseenv'].filter('c80i50 == "UniRef50_UPI00068EA29C"')._view_sequence()
temp._view_sequence()
cagE['ATPaseenv'].filter('c80i50 == "UniRef50_UPI00068EA29C"')._view_sequence()
cagE['it1atpsel_rescued'] = rdbs.concat([cagE['it1atpsel'],temp]).align()
cagE['it1atpsel_rescued'].view_sequence()
cagE['it1atpsel_rescued']._view_sequence()
cagE['it1atpsel_rescued'].sort([550])._view_sequence()
cagE['it1atpsel_rescued'].sort([12])._view_sequence()
cagE['it1atpsel_rescued'].sort([12,15])._view_sequence()
cagE['it1atpsel_rescued'].sort([12,18])._view_sequence()
cagE['it1atpsel_rescued'].sort([12,17])._view_sequence()
cagE['it1atpsel_rescued'].sort([17])._view_sequence()
cagE['it1atpsel_rescued'].sort([17,'length'])._view_sequence()
cagE['it1atpsel_rescued'].filter(minlen = 230).sort([17,'length']).align()._view_sequence()
cagE['it1atpsel_rescued'].filter(minlen = 250).sort([17,'length']).align()._view_sequence()
