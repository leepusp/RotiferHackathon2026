# IPython log file

cagE = dict(data=rdbs.sequence(["WP_000496009.1"], local_database_path=os.path.join(rotifer.config["data"],"fadb","nr","nr")))
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
get_ipython().run_line_magic('who', '')
1+1
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
mc
get_ipython().run_line_magic('who', '')
mc = uniprot.MappingCursor()
mc = uniprot.MappingCursor(host="davinci")
mc["UPI0000135940"]
rpf.print_everything()
from rotifer.pandas import functions as rpf
rpf.print_everything()
fc
fc.cursors['easel'].path
fc["UPI00313314C7"]
fc["UPI00313314C7"][0].seq
from Bio.SeqUtils
from Bio.SeqUtils import CheckSum
help(CheckSum)
get_ipython().run_line_magic('page', '(help(CheckSum))')
get_ipython().run_line_magic('pinfo2', 'CheckSum')
CheckSum.crc64(fc["UPI00313314C7"][0])
mc["UPI0000135940"]
mc["UPI0000135940"].head()
CheckSum.crc64(fc["UPI00313314C7"][0])
CheckSum.crc64(fc["UPI00313314C7"][0]).replace("CRC-","")
mc[CheckSum.crc64(fc["UPI00313314C7"][0]).replace("CRC-","")]
gnc
gnc.cursors['mirror'].path
gnc["BAA09657.1"]
ic = ncbi.IPGCursor()
ic["BAA09657.1"]
cagE['psi1u2n']
cagE['psi1u2n'].sample()
cagE['psi1u2n'].sample().target
cagE['psi1u2n'].sample().target.iloc[0]
gnc[cagE['psi1u2n'].sample().target.iloc[0]]
get_ipython().run_line_magic('page', '')
cagE['psi1u2n'].head(20)
cagE['psi1u2n'].target_type.value_counts()
cagE['psi1u2n'].target_type.drop_duplicates().sort_values()
cagE['psi1u2n'].target_type.drop_duplicates().sort_values(ascending=False)
cagE['psi1u2n'].sort_values(['target_type','accession','target'], ascending=[False,True,True]).drop_duplicates('source')
get_ipython().run_line_magic('page', '')
cagE['psi1ndf'] = gnc.fetchall(cagE['psi1u2n'].sort_values(['target_type','accession','target'], ascending=[False,True,True]).drop_duplicates('source').target.drop_duplicates().tolist())
gnc.missing
get_ipython().run_line_magic('page', '')
cagE['psi1ndf'].shape
cagE['psi1ndf']
get_ipython().run_line_magic('page', '')
cagE['psi1ndf'].pid.nunique()
get_ipython().run_line_magic('pinfo2', 'rdar.full_annotate')
get_ipython().run_line_magic('pinfo2', 'rdar.add_arch_to_df')
rdar.config
cagE['psi1ndf'] = rdar.add_arch_to_df(cagE['psi1ndf'])
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
cagE['psi1ndf'].head()
cagE['psi1ndf'].head(15)
get_ipython().run_line_magic('page', '')
# First, let us load the sequence, for simplicity, we are creating a dictionary
# We chose as an starting point the protein Q9ZLT4, a cagE entry
cagE = dict(data=rdbs.sequence(["WP_000496009.1"], local_database_path=os.path.join(rotifer.config["data"],"fadb","nr","nr")))
cagE._view_sequence()
cagE['data']._view_sequence()
rotifer.config["data"]
cagE = dict(data=rdbs.sequence(["Q9ZLT4"], local_database_path=os.path.join(rotifer.config["data"],"fadb","uniprot","uniref50")))
cagE['data']._view_sequence()
os.path.join(rotifer.config["data"],"fadb","uniprot","uniref50")
cagE = dict(data=rdbs.sequence(["Q9ZLT4"], local_database_path="/scratch/global/databases/fadb/uniprot/uniref50")))
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
cagE = dict(data=rdbs.sequence(["Q9ZLT4"], local_database_path="/scratch/global/databases/fadb/uniprot/uniref50"))
cagE['data']._view_sequence()
cagE = dict(data=rdbs.sequence(["UniRef50_Q9ZLT4"], local_database_path="/scratch/global/databases/fadb/uniprot/uniref50"))
cagE['data']._view_sequence()
cagE = dict(data=rdbs.sequence(["Q9ZLT4"], local_database_path="/scratch/global/databases/fadb/uniprot/uniprotDB"))
cagE['data']._view_sequence()
cagE = dict(data=rdbs.sequence(["Q9ZLT4"], local_database_path=os.path.join(rotifer.config["data"],"fadb","uniprot","uniprotDB")))
cagE['data']._view_sequence()
#
#cagE['data']._view_sequence()
# We can inspect the alignment by either running:
cagE['data'].view()
# or cagE['data']._view_sequence() on a terminal window
cagE["data"].slice((592,891)).view_sequence()
cagE["data"].slice((592,891))._view_sequence()
# psiblast search on the ATPase and N-terminal domains
cagE['qNterm'] = cagE["data"].slice((1,548))
cagE['qATP'] = cagE["data"].slice((592,891))
cagE["Ntermpsi1"] = rdar.psiblast(cagE["qNterm"], num_aln=5000, aln=False, db="/scratch/global/databases/fadb/uniprot/uniref50")
cagE["Atpasepsi1"] = rdar.psiblast(cagE["qATP"], num_aln=5000, aln=False, db="/scratch/global/databases/fadb/uniprot/uniref50")
cagE.keys()
get_ipython().system('ls -a')
get_ipython().system('less .rotifer/20261103123045/log.1.py')
get_ipython().run_line_magic('pwd', '')
cagE["Atpasepsi1"][0].head(10)
cagE["Ntermpsi1"][0].head(10)
cagE["Atpasepsi1"][0].hit.nunique()
len(cagE["Atpasepsi1"][0])
cagE["Atpasepsi1"][0].hit.value_counts()
cagE["Atpasepsi1"][0].query('hit == "UniRef50_UPI003C7329C5"')
get_ipython().run_line_magic('page', '(cagE["Atpasepsi1"][1])')
from IPython.core.page import page
page(cagE["Atpasepsi1"][1])
cagE["Atpasepsi1"][0].hit.nunique()
cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.nunique()
#Let's check if there are any truncated results here
cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.size()
#Let's check if there are any truncated results here
cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.count()
get_ipython().run_line_magic('pinfo2', 'rdae.extract_envelope')
get_ipython().run_line_magic('pinfo2', 'rdae.extract_envelope')
get_ipython().run_line_magic('pinfo', 'rdae.extract_envelope')
get_ipython().run_line_magic('pinfo2', 'rdae.extract_envelope')
cagE["Atpasepsi1"][0].query('evalue <= 0.1').head(10)
cagE["Atpasepsi1"][0].query('evalue <= 0.1').tail(10)
#Let's check if there are any truncated results here
cagE["Atpasepsi1"][0].query('evalue <= 0.01').hit.count()
#Let's check if there are any truncated results here
cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.count()
# Let's see distributions
cage['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1'), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
# Let's see distributions
cage['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseseqs'].to_hist()
cagE['ATPaseseqs'].to_hist()
# Let's see distributions
cage['ATPaseenv'] = rdae.extract_envelope(cagE['ATPaseseqs'][0].query('evalue <= 0.1'), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
# Let's see distributions
cage['ATPaseenv'] = rdae.extract_envelope(cagE['ATPaseseqs'], local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
get_ipython().run_line_magic('pinfo2', 'rdae.extract_envelope')
# Let's see distributions
cage['ATPaseenv'] = rdae.extract_envelope(cagE['ATPaseseqs'].df, local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseseqs'].df.id.value_counts().head(4)
get_ipython().run_line_magic('pinfo', 'rdae.extract_envelope')
# Let's see distributions
cage['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1'), seqobj = cagE['ATPaseseqs'], local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE["Atpasepsi1"][0].query('evalue <= 0.1').head(5)
# Let's see distributions
cage['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1'), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
# Let's see distributions
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1'), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
get_ipython().run_line_magic('pinfo2', 'rdae.extract_envelope')
import importlib
import linecache

linecache.clearcache()
importlib.reload(rdae)
get_ipython().run_line_magic('pinfo2', 'rdae.extract_envelope')
# Let's see distributions
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1'), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseenv']df.head(2)
cagE['ATPaseenv'].df.head(2)
cagE['ATPaseenv']._view_sequence()
cagE['Atpasepsi1'].df.head(2)
cagE['ATPaseseqs'].df.head(2)
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
# Let's see distributions
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1'), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseseqs'].df.head(2)
cagE['ATPaseseqs'].to_hist()
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseenv'].to_hist()
cagE['ATPaseenv'].df.columns
cagE['ATPaseenv'].df.pid.value_counts()
cagE['ATPaseenv'].df.pid.value_counts().head(6)
cagE["Atpasepsi1"][0].query('hit == "UniRef50_UPI000834519E"')
cagE['ATPaseenv'].df.id.value_counts().head(4)
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseseqs'].to_hist()
# Let's see distributions
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1'), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseenv'].to_hist()
cagE['ATPaseenv'].df.id.value_counts().head(4)
cagE["Atpasepsi1"][0].query('hit == "UniRef50_A0A6B8KKW5"')
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseseqs'].to_hist()
# Let's see distributions
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1').rename(columns={'hit': 'sequence','lenght':'og_length'}), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseenv'].to_hist()
cagE['ATPaseenv'].df.id.value_counts().head(4)
cagE['ATPaseenv'].align()._view_sequence()
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseseqs'].align()._view_sequence()
cagE['ATPaseenv'].align().sort(['length'])._view_sequence()
cagE['ATPaseenv'].align().sort([42])._view_sequence()
cagE['ATPaseseqs'].align()._view_sequence()
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
get_ipython().run_line_magic('pinfo2', 'rdae.extract_envelope')
# Let's see distributions
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1').rename(columns={'hit': 'sequence','lenght':'og_length'}), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], expand = 20, local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseenv'].to_hist()
cagE['ATPaseenv'].align()._view_sequence()
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
python --version
get_ipython().system('python --version')
cagE['ATPaseseqs'].align()._view_sequence()
cagE['ATPaseseqs'].align().sort([3930])._view_sequence()
cagE['ATPaseenv'].align()._view_sequence()
cagE['ATPaseenv'].align().sort([71])._view_sequence()
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
# Let's extract the envelope
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1').rename(columns={'hit': 'sequence','lenght':'og_length'}), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], expand = 50, local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseenv'].align()._view_sequence()
cagE['ATPaseenv'].align().sort([46])._view_sequence()
cagE['ATPaseenv'].align()._view_sequence()
cagE['ATPaseenv'].align().sort([200])._view_sequence()
cagE['ATPaseseqs'].filter(keep = cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list())align()._view_sequence()
cagE['ATPaseseqs'].filter(keep = cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list()).align()._view_sequence()
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseseqs'].filter(keep = cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list()).align()._view_sequence()
cagE['ATPaseseqs'].filter(keep = cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list()).align().sort([3930])._view_sequence()
# Let's extract the envelope
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1').rename(columns={'hit': 'sequence','lenght':'og_length'}), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], expand = 150, local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseenv'].align().sort([200])._view_sequence()
cagE['ATPaseenv'].align().sort([4700])._view_sequence()
cagE['ATPaseenv'].align().sort([470])._view_sequence()
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE["Atpasepsi1"][0].query('evalue <= 0.1').hist()
cagE["Atpasepsi1"][0].query('evalue <= 0.1')['querystart'].hist()
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseseqs'].filter(keep = cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list()).align().sort([3930])._view_sequence()
cagE['ATPaseseqs'].filter(keep = cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(),remove = ['UniRef50_F9ZMR9']).align().sort([3930])._view_sequence()
get_ipython().run_line_magic('pinfo', 'rdbs.sequence.filter')
cagE['ATPaseseqs'].filter(keep = cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list()).filter(remove = ['UniRef50_F9ZMR9']).align().sort([3930])._view_sequence()
# Let's extract the envelope
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1').rename(columns={'hit': 'sequence','lenght':'og_length'}), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], expand = 250, local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseenv'].align().sort([470])._view_sequence()
cagE['ATPaseenv'].align().sort([830])._view_sequence()
cagE['ATPaseenv'].align().sort([828])._view_sequence()
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseseqs'].filter(keep = cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list()).align().sort([3930])._view_sequence()
cagE["Atpasepsi1"][0].query('evalue <= 0.1')['querystart'].max()
# Let's extract the envelope
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1').rename(columns={'hit': 'sequence','lenght':'og_length'}), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], expand = 270, local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseenv'].align()._view_sequence()
cagE['ATPaseenv'].align().sort([1112])._view_sequence()
cagE['ATPaseenv'].align().sort([1113])._view_sequence()
cagE['ATPaseenv'].align().sort([1012])._view_sequence()
cagE['Atpasepsi1'][0].query('hit == "UniRef50_UPI00402ADB7E"')
cagE['Atpasepsi1'][0].query('hit == "UniRef50_F9ZMR9"')
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseseqs'].filter(keep = cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list()).align().sort([3930])._view_sequence()
cagE['ATPaseenv'].df.columns
cagE['ATPaseenv'].df[['id','hitstart','hitend','querystart']]
get_ipython().run_line_magic('page', '')
cagE["Atpasepsi1"][0].query('evalue <= 0.1')['querystart'].hist()
cagE["Atpasepsi1"][0].query("evalue <= 0.1").assign(late_start=hits.querystart >= 150,
    aligned_length=hits.queryend - hits.querystart + 1
).groupby("late_start")[
    ["querystart", "queryend", "aligned_length", "hitstart", "hitend"]
].describe()
cagE["Atpasepsi1"][0].query("evalue <= 0.1").assign(late_start=hits.querystart >= 150,
    aligned_length=cagE["Atpasepsi1"][0].query("evalue <= 0.1").queryend - cagE["Atpasepsi1"][0].query("evalue <= 0.1").querystart + 1
).groupby("late_start")[
    ["querystart", "queryend", "aligned_length", "hitstart", "hitend"]
].describe()
cagE["Atpasepsi1"][0].query("evalue <= 0.1").assign(late_start=cagE["Atpasepsi1"][0].query("evalue <= 0.1").querystart >= 150,
    aligned_length=cagE["Atpasepsi1"][0].query("evalue <= 0.1").queryend - cagE["Atpasepsi1"][0].query("evalue <= 0.1").querystart + 1
).groupby("late_start")[
    ["querystart", "queryend", "aligned_length", "hitstart", "hitend"]
].describe()
cagE['ATPaseenv'].df.id.head(4)
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseseqs'].filter(keep = cagE["Atpasepsi1"][0].query('evalue <= 0.1 and querystart >=150').hit.to_list()).align()._view_sequence()
cagE['ATPaseseqs'].filter(keep = cagE["Atpasepsi1"][0].query('evalue <= 0.1 and querystart >=150').hit.to_list()).align().sort([1370])._view_sequence()
cagE['ATPaseseqs'].filter(keep = cagE["Atpasepsi1"][0].query('evalue <= 0.1 and querystart >=150').hit.to_list()).align().sort([1370])._view_sequence()
cagE['ATPaseenv'].align().sort([1012])._view_sequence()
# Let's extract the envelope
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1').rename(columns={'hit': 'sequence','length':'og_length'}), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], expand = 270, local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseenv'].filter(f"pid in {cagE['Atpasepsi1'][0].query('evalue <= 0.1 and querystart >= 150').hit.to_list()!r}").align().sort([1370])._view_sequence()
cagE['ATPaseenv'].filter(f"pid in {cagE['Atpasepsi1'][0].query('evalue <= 0.1 and querystart >= 150').hit.to_list()!r}").align().sort([261])._view_sequence()
cagE['ATPaseenv'].align()._view_sequence()
get_ipython().run_line_magic('pinfo2', 'sort')
get_ipython().run_line_magic('pinfo2', 'sort')
get_ipython().run_line_magic('pinfo', 'rdbs.sequence.sort')
cagE['ATPaseenv'].align().sort(['list'],id_list=['UniRef50_Q48252/1-300'])._view_sequence()
cagE['ATPaseenv'].align().sort(['list'],id_list=['UniRef50_Q48252/592-891'])._view_sequence()
cagE['ATPaseenv'].align().sort(['list'],id_list=['UniRef50_Q48252/592-891']).slice((1000,2445))._view_sequence()
cagE['ATPaseenv'].align().sort(['list'],id_list=['UniRef50_Q48252/592-891']).slice((1000,2445)).align()._view_sequence()
cagE['ATPaseenv'] = cagE['ATPaseenv'].align().sort(['list'],id_list=['UniRef50_Q48252/592-891']).slice((1000,2445)).align()
cagE['ATPaseenv'].add_cluster(cascade=[(100,100),(80,70),(80,50),(80,30),(0,0)], inplace=True)
cagE['ATPaseenv'].to_hist()
cagE['ATPaseenv']df.describe()
cagE['ATPaseenv'].df.describe()
cagE['ATPaseenv'].df.nunique()
cagE['ATPaseenv'].df.c0i0.value_counts()
cagE['ATPaseenv'].df.c80i30.value_counts()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A8T3YM87" and id == c80i50').align()._view_sequence()
cagE['ATPaseenv'].df.c0i0.value_counts()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align()._view_sequence()
get_ipython().run_line_magic('pinfo2', 'rdar.rdar.add_arch_to_df')
get_ipython().run_line_magic('pinfo2', 'rdar.add_arch_to_df')
get_ipython().run_line_magic('pinfo2', 'rdae.add_arch_to_df')
get_ipython().run_line_magic('pinfo2', 'rdae.add_arch_to_df')
cagE['Atpasehs'] = rdbh.hmmscan(cagE['ATPaseenv'], models_path="/scratch/global/databases/pfam/Pfam-A.hmm")
cagE['Atpasepfam'] = riu.architecture(rdbh.hmmscan(cagE['Atpasehs'], models_path="/scratch/global/databases/pfam/Pfam-A.hmm").query('c_evalue <= 0.1'), **riu.config['hmmer'])
cagE['Atpasehs'].head()
cagE['Atpasepfam'] = riu.architecture(cagE['Atpasehs'].query('c_evalue <= 0.1'), **riu.config['hmmer'])
cagE['Atpasepfam'].head()
cagE['Atpasepfam'].architecture.value_counts()
cagE['Atpasehs'].head()
cagE['Atpasehs'].sequence.nunique()
cagE['ATPaseenv'].df.nunique()
cagE['Atpasehs'].sequence.nunique()
display(
    cagE['Atpasehs'].assign(
        actual_length=lambda d: d.sequence.map(
            cagE['ATPaseenv'].df.set_index('id').sequence.str.len()
        )
    )[['sequence', 'actual_length', 'aln_target_length', 'estart', 'eend']].head(10)
)
get_ipython().run_line_magic('pinfo2', 'rdbh.hmmscan')
get_ipython().run_line_magic('pinfo2', 'rdbh.hmmscan')
cagE['Atpasehs'].assign(
    actual_length=lambda d: d.sequence.map(
        cagE['ATPaseenv'].df.set_index('id').sequence.str.len()
    )
)[
    ['sequence', 'actual_length', 'aln_target_length', 'estart', 'eend']
].head(10)
cagE['ATPaseenv'].to_hist()
cagE['ATPaseenv'].df.assign(
    actual_length=lambda d: d.sequence.str.len(),
    gaps=lambda d: d.sequence.str.count('-')
)[['id', 'length', 'actual_length', 'gaps']].head(15)
get_ipython().run_line_magic('pinfo2', 'rdbh._load_sequences')
get_ipython().run_line_magic('pinfo2', 'rdbh.hmmscan')
get_ipython().run_line_magic('pinfo2', 'rdbh._load_sequences')
import importlib
import linecache

linecache.clearcache()
importlib.reload(rdbh)
get_ipython().run_line_magic('pinfo2', 'rdbh._load_sequences')
cagE['Atpasehs'] = rdbh.hmmscan(cagE['ATPaseenv'], models_path="/scratch/global/databases/pfam/Pfam-A.hmm")
cagE['Atpasehs'].head(40)
get_ipython().run_line_magic('page', '')
cagE['Atpaseungap'] = cagE['ATPaseenv'].copy()
cagE['Atpaseungap'].df['sequence'] = cagE['Atpaseungap'].df.sequence.str.replace('-', '', regex=False) 
cagE['Atpasehs'] = rdbh.hmmscan(cagE['ATPaseungap'], models_path="/scratch/global/databases/pfam/Pfam-A.hmm")
cagE['Atpasehs'] = rdbh.hmmscan(cagE['Atpaseungap'], models_path="/scratch/global/databases/pfam/Pfam-A.hmm")
cagE['Atpasehs'].head(40)
get_ipython().run_line_magic('page', '')
gaphscan = rdbh.hmmscan(cagE['ATPaseenv'], models_path="/scratch/global/databases/pfam/Pfam-A.hmm")
gaphscan.head(10)
cagE['Atpasehs'].head(10)
cagE['Atpasehs'].sequence.nunique()
cagE['Atpasepfam'] = riu.architecture(cagE['Atpasehs'].query('c_evalue <= 0.1'), **riu.config['hmmer'])
cagE['Atpasepfam'].architecture.value_counts()
cagE['Atpasepfam'].architecture.value_counts().head(5)
cagE['ATPaseenv'].df['evalue'] = cagE['psi1seqs'].df.id.map(cagE['psi1'][0].sort_values(['evalue']).drop_duplicates('hit').set_index('hit').evalue)
cagE['ATPaseenv'].df['pfam'] = cagE['ATPaseenv'].df.id.map(cagE['Atpasepfam'].set_index('sequence').architecture)
cagE['ATPaseenv']._view_sequence()
cagE['ATPaseenv'].sort(['evalue','pfam'])._view_sequence()
cagE['ATPaseenv'].sort(['evalue','pfam'])._view_groups(groupby='c0i0',sample=5)
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align()._view_sequence()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().filter('cagE["ATPaseenv"].df.sequence[10] =="-"')._view_sequence()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().filter(f'cagE["ATPaseenv"].df.sequence[10] =="-"')._view_sequence()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().filter('sequence[10] =="-"')._view_sequence()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().filter('sequence.str[10]=="-"')._view_sequence()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().filter('sequence.str[10]=="-"').sort(['evalue'])._view_sequence()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().filter(f'cagE["ATPaseenv"].df.sequence[10] =="-"').c80i50
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().filter(f'cagE["ATPaseenv"].df.sequence.str[10] =="-"').df.c80i50
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().filter('cagE["ATPaseenv"].df.sequence.str[10] =="-"').df.c80i50
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().filter('sequence.str[10] =="-"').df.c80i50
cagE['ATPaseenv'].filter(' c80i50 == "UniRef50_A0A378QLD4/1-167"').align()._view_sequence()
cagE['ATPaseenv'].filter('c80i50 == "UniRef50_A0A090MWC5/8-253"').align()._view_sequence()
cagE['ATPaseenv'].filter('c80i50 == "UniRef50_A0A090MWC5/8-253"').align().sort(['evalue'])._view_sequence()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').drop_duplicates('c80i50').id.tolist()}").align()._view_sequence()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').drop_duplicates('c80i50').id.tolist()}").align().sort([10])._view_sequence()
cagE['ATPaseenv'].filter('c80i50 == "  UniRef50_A0A3B0S2C7/444-698"').align()._view_sequence()
cagE['ATPaseenv'].filter('c80i50 == "UniRef50_A0A3B0S2C7/444-698"').align()._view_sequence()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().sort([10])._view_sequence()
cagE['ATPaseenv'].c80i50.value_counts(15)
cagE['ATPaseenv'].df.c80i50.value_counts(15)
cagE['ATPaseenv'].df.c80i50.value_counts().head(15)
cagE['ATPaseenv'].df.c80i50.nunique()
cagE['ATPaseenv'].df.c0i0.value_counts().head(15)
cagE['ATPaseenv'].df.c0i0.value_counts()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A519BBK4/467-757" and id == c80i50').align().sort()._view_sequence()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0A519BBK4/467-757" and id == c80i50').align().sort([10])._view_sequence()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0A519BBK4/467-757"').sort_values('evalue').drop_duplicates('c80i50').id.tolist()}").align().sort([10])._view_sequence()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().sort([10]).df.evalue.max()
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().sort([10,'evalue'])
cagE['ATPaseenv'].filter('c0i0 == "UniRef50_A0AAW9DAR3/433-671" and id == c80i50').align().sort([10,'evalue'])._view_sequence()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0A519BBK4/467-757"').sort_values('evalue').drop_duplicates('c80i50').id.tolist()}").align().sort([10])._view_sequence()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').drop_duplicates('c80i50').id.tolist()}").align().sort([10])._view_sequence()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').drop_duplicates('c80i50').id.tolist()}").align().sort(['evalue'])._view_sequence()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue')..groupby('c80i50').head(2).id.tolist()}").align().sort(['evalue'])._view_sequence()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(2).id.tolist()}").align().sort(['evalue'])._view_sequence()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(2).id.tolist()}").align().sort(['evalue']).filter('sequence.str[9] != "-"').sort(['evalue'])._view_sequence()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(2).id.tolist()}").align().sort(['evalue']).filter('sequence.str[9] != "-"').sort(['evalue']).c80i50.nunique()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(2).id.tolist()}").align().sort(['evalue']).filter('sequence.str[9] != "-"').sort(['evalue']).df.c80i50.nunique()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(2).id.tolist()}").align().sort(['evalue']).df.c80i50.nunique()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(3).id.tolist()}").align().sort(['evalue']).filter('sequence.str[9] != "-"').sort(['evalue']).df.c80i50.nunique()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(5).id.tolist()}").align().sort(['evalue']).filter('sequence.str[9] != "-"').sort(['evalue']).df.c80i50.nunique()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(15).id.tolist()}").align().sort(['evalue']).filter('sequence.str[9] != "-"').sort(['evalue']).df.c80i50.nunique()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(2).id.tolist()}").align().sort(['evalue']).filter('sequence.str[9] != "-"').sort(['evalue']).df.c80i50.nunique()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(2).id.tolist()}").align().sort(['evalue']).filter('sequence.str[9] != "-"').align().sort(['evalue'])._view_sequence()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(2).id.tolist()}").align().sort(['evalue']).filter('sequence.str[9] != "-"').align().sort(['evalue']).filter('~c80i50.duplicated()')._view_sequence()
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(2).id.tolist()}").align().sort(['evalue']).filter('sequence.str[9] != "-"').align().sort(['evalue']).filter('~c80i50.duplicated()').df.c80i50.value_counts().head(4)
cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(2).id.tolist()}").align().sort(['evalue']).filter('sequence.str[9] != "-"').align().sort(['evalue']).filter('~c80i50.duplicated()').align()._view_sequence()
cagE['ATPasec0i01'] = cagE['ATPaseenv'].filter(f"id in {cagE['ATPaseenv'].df.query('c0i0 == "UniRef50_A0AAW9DAR3/433-671"').sort_values('evalue').groupby('c80i50').head(2).id.tolist()}").align().sort(['evalue']).filter('sequence.str[9] != "-"').align().sort(['evalue']).filter('~c80i50.duplicated()').align()
cagE['ATPasec0i01']._view_sequence()
cagE['ATPasec0i01'].sort([452])._view_sequence()
cagE['ATPasec0i01'].sort([453])._view_sequence()
cagE['ATPasec0i01'].filter('sequence.str[454] != "-"').sort([453])._view_sequence()
cagE['ATPasec0i01'].filter('sequence.str[454] != "-"').filter(remove = ['UniRef50_UPI0029C012DE/63-279']).align().sort(['evalue'])._view_sequence()
cagE['ATPasec0i01'].filter('sequence.str[454] != "-"').filter(remove = ['UniRef50_UPI0029C012DE/63-279','UniRef50_A0A1R3L3S7/46-244']).align().sort(['evalue'])._view_sequence()
cagE['ATPasec0i01'].filter('sequence.str[454] != "-"').filter(remove = ['UniRef50_UPI0029C012DE/63-279','UniRef50_A0A1R3L3S7/46-244']).align().slice((2,-2)).sort(['evalue'])._view_sequence()
cagE['ATPasec0i01'] = cagE['ATPasec0i01'].filter('sequence.str[454] != "-"').filter(remove = ['UniRef50_UPI0029C012DE/63-279','UniRef50_A0A1R3L3S7/46-244']).align().slice((2,-2)).sort(['evalue'])
cagE['Ntermseqs'] = rdbs.sequence(cagE["Ntermpsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['Ntermseqs'].align()._view_sequence()
cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.nunique()
cagE["Ntermpsi1"][0].query('evalue <= 0.1').hit.nunique()
cagE["Ntermpsi1"][0].query('evalue <= 0.1').hit.count()
cagE["Ntermpsi1"][0].query('evalue <= 0.5').hit.nunique()
cagE["Ntermpsi1"][0].hit.nunique()
cagE["Ntermpsi1"][0].query('evalue <= 0.1').hit.nunique()
cagE["Ntermpsi1"][0].query('evalue <= 0.1').hit.count()
cagE['ATPasec0i01'].df.pfam.value_counts()
cagE["Ntermpsi1"][0].query('evalue <= 0.1').queryend.max()
cagE["Ntermpsi1"][0].query('evalue <= 0.1').querystart.max()
cagE["Ntermpsi1"][0].query('evalue <= 0.1').querystart.min()
cagE["Ntermpsi1"][0].query('evalue <= 0.1').querystart.value_counts()
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Ntermpsi1"][0].query('evalue <= 0.1').rename(columns={'hit': 'sequence','length':'og_length'}), start ='hitstart', end = 'hitend', seqobj = cagE['Ntermseqs'], expand = 10, local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
# Let's see distributions
cagE['ATPaseseqs'] = rdbs.sequence(cagE["Atpasepsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
# Let's extract the envelope, we are being generous because the insertion
cagE['ATPaseenv'] = rdae.extract_envelope(cagE["Atpasepsi1"][0].query('evalue <= 0.1').rename(columns={'hit': 'sequence','length':'og_length'}), start ='hitstart', end = 'hitend', seqobj = cagE['ATPaseseqs'], expand = 270, local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['ATPaseenv'] = cagE['ATPaseenv'].align().sort(['list'],id_list=['UniRef50_Q48252/592-891']).slice((1000,2445)).align() #trimming down to the atpase domain
cagE['ATPaseenv'].to_hist()
cagE['ATPaseenv'].add_cluster(cascade=[(100,100),(80,70),(80,50),(80,30),(0,0)], inplace=True)
cagE['ATPaseenv'].df.c0i0.value_counts()
cagE['ATPaseenv']._view_sequence()
cagE['ATPaseenv'].df.nunique()
cagE['ATPaseenv'].to_hist()
cagE['Atpaseungap'].align()._view_sequence()
cagE['ATPaseenv']._view_sequence()
cagE['NTermenv'] = rdae.extract_envelope(cagE["Ntermpsi1"][0].query('evalue <= 0.1').rename(columns={'hit': 'sequence','length':'og_length'}), start ='hitstart', end = 'hitend', seqobj = cagE['Ntermseqs'], expand = 10, local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['Ntermseqs'] = rdbs.sequence(cagE["Ntermpsi1"][0].query('evalue <= 0.1').hit.to_list(), local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
cagE['NTermenv'] = rdae.extract_envelope(cagE["Ntermpsi1"][0].query('evalue <= 0.1').rename(columns={'hit': 'sequence','length':'og_length'}), start ='hitstart', end = 'hitend', seqobj = cagE['Ntermseqs'], expand = 10, local_database_path='/scratch/global/databases/fadb/uniprot/uniref50')
# Fiz caca, deixa eu tentar ver o que rolou

old = cagE['Atpaseungap'].df.set_index('id')
new = cagE['ATPaseenv'].df.set_index('id')

print("Same IDs:", set(old.index) == set(new.index))

common = old.index.intersection(new.index)

print("Different amino acid sequences:", (
    old.loc[common, 'sequence'].str.replace('-', '', regex=False)
    != new.loc[common, 'sequence'].str.replace('-', '', regex=False)
).sum())
pd.crosstab(
    old.loc[common, 'c0i0'],
    new.loc[common, 'c0i0']
)
import pandas as pd
pd.crosstab(old.loc[common, 'c0i0'], new.loc[common, 'c0i0'])
from sklearn.metrics import adjusted_rand_score

old = cagE['Atpaseungap'].df.set_index('id')
new = cagE['ATPaseenv'].df.set_index('id')

print(
    "Same row order:",
    cagE['Atpaseungap'].df.id.reset_index(drop=True).equals(
        cagE['ATPaseenv'].df.id.reset_index(drop=True)
    )
)

for col in ['c100i100', 'c80i70', 'c80i50', 'c80i30', 'c0i0']:
    print(
        col,
        "Old:", old[col].nunique(),
        "New:", new[col].nunique(),
        "ARI:", round(adjusted_rand_score(
            old[col].fillna('MISSING'),
            new.loc[old.index, col].fillna('MISSING')
        ), 3)
    )
(
    old['c100i100'] !=
    new.loc[old.index, 'c100i100']
).sum()
import copy

a = copy.deepcopy(cagE['Atpaseungap'])
b = copy.deepcopy(cagE['ATPaseenv'])

a.df = a.df.sort_values('id').reset_index(drop=True)
b.df = b.df.sort_values('id').reset_index(drop=True)

cascade = [(100,100),(80,70),(80,50),(80,30),(0,0)]

a.add_cluster(cascade=cascade, inplace=True)
b.add_cluster(cascade=cascade, inplace=True)

for col in ['c80i70', 'c80i50', 'c80i30', 'c0i0']:
    print(
        col,
        adjusted_rand_score(
            a.df.set_index('id')[col],
            b.df.set_index('id').loc[a.df.id, col]
        )
    )
exit()
