
# Adding module location to sys.path
import sys
sys.path.append('/home/oth/anaconda3/lib/python3.8/site-packages')
sys.path.append('/home/oth/BHS/PYT/1.project/SBI_PYT/PPI_main') # Path to I_O_args script
sys.path.append('/home/oth/BHS/PYT/1.project/dash_test.py')
sys.path.append('/home/oth/BHS/PYT/1.project/test_project.py')
import os
import numpy
import scipy
import parser as args
from Bio.PDB import *
from Bio.Seq import Seq
from Bio import pairwise2 as pw2

# DASH IMPORTS
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from urllib.parse import quote as urlquote
from flask import Flask, send_from_directory
import plotly.express as px
import webbrowser
import base64
import io
import datetime

#############
#   DASH    #
#############
# -*- coding: utf-8 -*-

# UPLOADING, DOWNLOADING FILES
#######################################################################

# Creating the upload directory, if it does not exist
UPLOAD_DIRECTORY = "./project/app_uploaded_files"

if not os.path.exists(UPLOAD_DIRECTORY):
    os.makedirs(UPLOAD_DIRECTORY)

# Flask server creation
server = Flask(__name__)
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css'] # External CSS stylesheet
app = dash.Dash(__name__, external_stylesheets=external_stylesheets, server=server)

@server.route("/download/<path:path>")
def download(path):
    return send_from_directory(UPLOAD_DIRECTORY, path, as_attachment=True)


# LAYOUT
#######################################################################

app.layout = html.Div([
    dcc.Markdown(children="""
    #     **pytprot**
    ## A web-app for multicomplex building
    """, style={'textAlign' : 'left', 'marginLeft': 50, 'marginTop': 80, 'marginRight': 10,
                'color':'black', 'width': '80%'}),

    # MODEL INPUT

    html.Div([
        dcc.Markdown(children="""
        ##### Introduce the **pdb files** to build the model with:
        """, style={'width': '80%', 'marginLeft': 50, 'textAlign': 'left', 'display':'inline-block'}),

        # UPLOAD
        html.Div([
                     dcc.Upload(
                id='upload-data',
                children=html.Div(['Drag and Drop or ', html.A('Select PDB Files')]),
                style={'width': '70%',
                       'height': '60px',
                       'lineHeight': '60px',
                       'borderWidth': '1px',
                       'borderStyle': 'dashed',
                       'borderRadius': '5px',
                       'textAlign': 'center',
                       'margin': '10px',
                       'marginLeft': 50,
                       'background-color': 'white'},
                multiple=True,
                ),
        dcc.Markdown(children="""
        #### **File list**
        """,
                     style={'marginLeft': 50}),
        html.Ul(id="file-list", style={'marginLeft': 50}),
        html.Button('Run complex building', id='submit-run', n_clicks=0,
                    style={'marginTop': 20, 'marginLeft': 50})
        ],
            ),

    # RUN COMPLEX BUILDING Button
    html.Div(id="start-button", children='Run model',
                      style={'marginTop': 10, 'marginLeft': 50, 'color':'red'})

        ]),

    html.Br(),

    # PARAMETERS
    #
    dcc.Markdown(children="""Specify the different **parameters**:""",
                 style={'textAlign': 'left', 'marginLeft': 50, 'marginTop': 30}),

    dcc.Checklist(
        options=[
            {'label': 'Stoichiometry', 'value': 'sto'},
            {'label': 'Forcing', 'value': 'for'},
            {'label': 'verbose', 'value': 'ver'}
        ],
        labelStyle={'float': 'center', 'display':'inline-block', 'marginLeft': 50, 'marginBottom': 300})
    ],
    style={
        'columnCount': 1,
        'background-color':'beige'})



# FUNCTIONS
#######################################################################
port = 8000

def open_browser():
    """
    Everytime the code is run, the browser opens
    :return:
    """
    webbrowser.open_new("http://localhost:{}".format(port))

def save_file(name, content):
    """Decodes, stores uploaded files"""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(os.path.join(UPLOAD_DIRECTORY, name), "wb") as fp:
        fp.write(base64.decodebytes(data))

def uploaded_files():
    """List files in the upload directory."""
    files = []
    for filename in os.listdir(UPLOAD_DIRECTORY):
        path = os.path.join(UPLOAD_DIRECTORY, filename)
        if os.path.isfile(path):
            files.append(filename)
    return files

def file_download_link(filename):
    """A Plotly Dash that downloads a file"""
    location = "/download/{}".format(urlquote(filename))
    return html.A(filename, href=location)

# CALLBACK
#######################################################################

@app.callback(
    dash.dependencies.Output('start-button', 'children'),
    [dash.dependencies.Input('submit-run', 'n_clicks')]
)
def main_code_exec(n_clicks):
    if n_clicks > 0:
        with open(os.getcwd()) as infile:
            exec(infile.read())
    if n_clicks > 2:
        return f"You are already running a macrocomplex. To run another one, you must wait until the current one is finished."


input_files = []

@app.callback(
    Output("file-list", "children"),
    [Input("upload-data", "filename"), Input("upload-data", "contents")],
)
def update_output(uploaded_filenames, uploaded_file_contents):
    """Save uploaded files and regenerate the file list."""
    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file(name, data)

    global input_files
    input_files = [file for file in uploaded_files()]

    if len(input_files) == 0:
        return "No files found!"
    else:
        return input_files

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

def input_ids():
    print("These are the input files:\n")
    print(input_files)

print(f"#######################\nYou are now running the pytprot module..........\n#######################")
input_ids()
#pdblist = update_output()
#print(pdblist)

# ARGPARSER
#options = args.output_argparser() # Object that contains all argument parsing functions
#a = args.check_infile_names(options.infile) # Checking correct file-naming



# DUMMY INPUT (This works, as a backup)
inputfiles_basic = [file for file in os.listdir('/home/oth/BHS/PYT/1.project/SBI_PYT/tests/examples/1gzx')]
path_basic = '/home/oth/BHS/PYT/1.project/SBI_PYT/tests/examples/1gzx/'


# List of files
pdblist = []
#for pdbfile in os.listdir(infiles_dash):
#    pdblist.append(str(options.infile+pdbfile)) # List of PDBs with their path

for file in inputfiles_basic:
    pdblist.append(str(path_basic+file)) # List of PDBs with their path

print(pdblist)



# PDBLIST IS THE OBJECT USED AS AN INPUT


# Dictionary with files and their PDB models
i = 0
dictmodels = {}
for pairs in pdblist:
    #final_name = pairs.split("/")[-1] # Better to not split it, actually
    final_name = pairs
    dictmodels[final_name]= []
    parser = PDBParser(PERMISSIVE=True, QUIET=False)
    structure = parser.get_structure(f"TEST_{i}", pairs)
    model = structure[0]
    model.id = i
    dictmodels[final_name] = model # Path+file as KEY, model object as VALUE
    i += 1

#print(dictmodels.items())


# Stechiometry reading

"""stechio_dict = {}
stech = options.stechiometry
total_subunits = 0

try:
    for line in open(stech):
        letter = line.split(":")[0].rstrip()
        num = int(line.split(":")[1].rstrip())
        stechio_dict[letter] = num
        total_subunits = sum(stechio_dict.values())
except TypeError:
    pass"""

####################
#   FUNCTIONS      #
####################


def similar_chains(input_pdb_dict):
    """
    Given an dictionary of PDB files as Keys, and their Structure objects as values, this function
    makes pairwise alignments between the chains, keeping only those with a 95% or higher
    similarity.
    """
    similar_chains = {}
    ppb = PPBuilder()
    model_list = [model for model in input_pdb_dict.values()]
    total_chains = []
    for model in model_list:
        for chain in model:
            total_chains.append(chain)
    for c, chain in enumerate(total_chains):
        for pp in ppb.build_peptides(chain):
            seq1 = pp.get_sequence()
            #print(chain)
            for chain2 in total_chains[c+1:]:
                if chain2.id != chain.id:
                    for pp in ppb.build_peptides(chain2):
                        seq2 = pp.get_sequence()
                        #print(f"\t\t{chain2}")
                        alignments = pw2.align.globalxx(seq1, seq2, score_only = True)
                        #print(pw2.format_alignment(*alignments[0]))
                        max_length = min(len(seq1), len(seq2))
                        identity_perc = round(alignments / max_length, 2)
                        #print(identity_perc)
                        if identity_perc > 0.95:
                            similar_chains.setdefault(chain, chain2)
    return similar_chains


# Alignment funciton test:
#for x, y in similar_chains(dictmodels).items():
#    print(x, y)

# IMPROVEMENTS:
# Crear un archivo con los alineamientos. O a lo mejor ocupa mucha memoria. Ponerlo como optional argument.
# INPUT AS A LIST, this way it's easire to avoid repeating unnecessary comparisons


def check_type(input_pdb_file):
    """
    Given a PDB file, it checks if it pertains to a DNA or protein structure.
    :param input_pdb_file:
    :return: structure object
    """
    i = 0
    final_name = input_pdb_file.split("/")[-1] # Only keep name, not full path
    structure = PDBParser(PERMISSIVE=True, QUIET=True).get_structure(f"model_{i}", input_pdb_file)
    atoms = structure[0].get_atoms()
    for x in atoms:
        if x.get_id() == "CA":
            print(f"{final_name} is a PROTEIN file") # Print optional
            return f"Protein"
        if x.get_id() == "P":
            print(f"{final_name} is a DNA file") # Print optional
            return f"DNA"

def change_chain_id(input_pdb_files):
    """
    Given a list of .pdb files, the Chain IDs are converted from a letter format to a numeric format
    :param input_pdb_file: 
    :return: A list of Structures
    """
    struct_list = []
    i = 1
    j = 0
    while i < len(input_pdb_files)*2: # As we will always have 2 chains per file
        for file in input_pdb_files:
            structure = PDBParser(PERMISSIVE=True, QUIET=True).get_structure(f"model_{j}", file)
            chains = structure[0].get_chains()
            j += 1
            for ch in chains:
                ch.id = i # replace chain ID with number
                i += 1 # i will keep increasing depending on the number of models
            struct_list.append(structure)
            break

    # IMPROVEMENT: Cambiar la letra por su posición en el abecedario? De esa manera
    # si tenemos varias "A" todas se cambiarían por 1, pero no sé si esto es
    # realmente necesario o no

    return struct_list

# Chain ID change test
#
"""for x in change_chain_id(pdblist):
    print(x.get_id())
    for ch in x.get_chains():
        print(ch.get_id())"""

"""for x,y in dictmodels.items():
    print(x)
    for ch in y.get_chains():
        print(ch)"""


models = [model for id, model in dictmodels.items()]

def common_chain_res(chain1, chain2):
    """
    Given a pair of chains, obtains the common residues and their respective CA atoms.
    Returns a tuple with two lists: The first one contains the list of atoms corresponding
    to the first chain. The second list contains the list of atoms of the second chain.
    :param chain1:
    :param chain2:
    :return:
    """
    res1 = [res for res in chain1 if res["CA"]]
    res2 = [res for res in chain2 if res["CA"]]
    common_res1 = [res1 for res1, res2 in zip(res1, res2) if res1.get_resname() == res2.get_resname()]
    common_res2 = [res2 for res1, res2 in zip(res1, res2) if res1.get_resname() == res2.get_resname()]

    chain1_atoms_list = [res["CA"] for res in common_res1]
    chain2_atoms_list = [res["CA"] for res in common_res2]
    common_atoms = (chain1_atoms_list, chain2_atoms_list)

    return common_atoms

# Common atoms test

chain1 = []
for chains in models[0].get_chains():
    chain1.append(chains)

chain2 = []
for chains in models[1].get_chains():
    chain2.append(chains)

#common_atoms = common_chain_res(chain2[1], chain1[1]) # Between chains C and B


def superimpose(chain1, chain2):
    """
    Given a pair of chains, these are superimposed. The first input is the
    fixed chain, and the second one the moving chain. The rotran matrix is computed, applied to the
    second chain, and then the RMSD is computed.
    :param :
    :return: RMSD score between structures
    """

    # Set fixed, moving models
    fixed_chain = chain1
    moving_chain = chain2

    # Create the LIST OF ATOMS to be aligned
    fixed_atoms = [res["CA"] for res in chain1]
    moving_atoms = [res["CA"] for res in chain2]

    # When superimposing chains are not equally sized
    if len(fixed_atoms) != len(moving_atoms):
        common_atoms = common_chain_res(fixed_chain, moving_chain)

        imposer = Superimposer()
        imposer.set_atoms(common_atoms[0], common_atoms[1])
        imposer.apply(fixed_chain.get_atoms())
        print(imposer.rms)

    # If they are the same size
    else:
        imposer = Superimposer()
        imposer.set_atoms(fixed_atoms, moving_atoms)
        imposer.apply(fixed_chain.get_atoms())
        print(imposer.rms)

# Superimpose test
#superimpose(chain2[1], chain1[0])


# FOR THE DASH APP TEST #
#if __name__ == '__main__':
    #dt.Timer(10, dt.open_browser).start()
    #dt.app.run_server(debug=True, port=dt.port) # Hot reloading
#    app.run_server(debug=True)
#    print("Running the server")


from modeller import *






