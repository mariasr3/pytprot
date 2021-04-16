#!/usr/bin/python3

# Adding module location to sys.path
import sys, os, re, numpy, scipy, io, datetime, webbrowser, base64, time, shutil
#sys.path.append('/home/oth/anaconda3/lib/python3.8/site-packages/')
#sys.path.append('/home/oth/BHS/PYT/1.project/pytprot/pytprot') # Path to the main script folder
#sys.path.append('/home/oth/BHS/PYT/1.project/dash_test.py')
#sys.path.append('/home/oth/BHS/PYT/1.project/test_project.py')
from Bio.PDB import *
from Bio.Seq import Seq
from Bio import pairwise2 as pw2
import inputfunctions, modelfunctions, chainfunctions

import dash, dash_component_unload
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from urllib.parse import quote as urlquote
from flask import Flask, send_from_directory
#import plotly.express as px

#############
#   DASH    #
#############
# -*- coding: utf-8 -*-

# UPLOADING, DOWNLOADING FILES
#######################################################################

# Creating the upload directory, if it does not exist
UPLOAD_DIRECTORY = "./app_uploaded_files/"
FINAL_MODEL = UPLOAD_DIRECTORY + "built_models"

if os.path.isdir(UPLOAD_DIRECTORY):
    shutil.rmtree(UPLOAD_DIRECTORY)
elif not os.path.exists(UPLOAD_DIRECTORY):
    os.makedirs(UPLOAD_DIRECTORY)

if not os.path.exists(FINAL_MODEL):
    os.makedirs(FINAL_MODEL)



# Flask server creation
server = Flask(__name__)
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css'] # External CSS stylesheet
app = dash.Dash(__name__, external_stylesheets=external_stylesheets, server=server)

@app.server.route("/downloads/<path:path>")
def download(path):
    root_dir = os.getcwd()
    return send_from_directory(os.path.join(root_dir, 'downloads'), path, as_attachment=True)


# LAYOUT
#######################################################################

app.layout = html.Div([
    dcc.Markdown(children="""
    #     **pytprot**
    ## A web-app for multicomplex building
    """, style={'textAlign' : 'left', 'marginLeft': 50, 'marginTop': 80, 'marginRight': 10,
                'color':'black', 'width': '80%'}),

    html.Hr(style={"width": "80%", 'borderColor': 'black', 'align':'left'}),

    # MODEL INPUT

    html.Div([
        dcc.Markdown(children="""
        ##### Introduce the **pdb files** to build the model:
        **(Optional)** You can also include the  **Stoichiometry .txt** file here
        """, style={'width': '80%', 'marginLeft': 50, 'textAlign': 'left', 'display':'inline-block'}),

        # UPLOAD
        html.Div([
                     dcc.Upload(
                id='upload-data',
                children=html.Div(['Drag and Drop or ', html.A('Select PDB (and .txt) files')]),
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

            html.Br(),

            html.Div([
                dcc.Markdown(children=""" 
                ###### Multicomplex assembly parameters:
                **NOTE:** If specific values are wanted, these must be indicated **before** loading the PDB files. 
                """),

                html.Div([

                    html.Div([
                        dcc.Markdown(children=""" **Contact distance: ** """),
                        dcc.Input(id="dist-contact", type="number"),
                    ]),

                    html.Div([
                        dcc.Markdown(children="""**Number of contacts: **"""),
                        dcc.Input(id="num-contact", type="number"),
                    ])

                ], style={'columnCount': 2}),

                html.Div(id="parameters-out"),

                dcc.Markdown(children="""*By default, the contact distance is set to 12 (Å) 
                                        and the number of contacts is set to 8 contacts*"""),

            ], style={'marginLeft': 50, 'columnCount': 1, 'background-color':'beige', 'width':'70%'}),


        html.Br(),


        dcc.Markdown(children=""" **Upload information** """, style={'marginLeft': 50}),
        html.Ul(id="file-list",
                style={'marginLeft': 50,
                       'background-color': 'white',
                       'width':'30%',
                       'borderWidth': '1px',
                       'borderStyle':'dashed',
                       'padding':'10px'
                       })],
            ),

        html.Br(),

        dcc.Markdown(children=""" 
        ##### **Models built** 
        *In order to save any model, you can right-click on the name > Save link as...*
        """, style={'marginLeft': 50}),

        html.Ul(children="PDB files processed", id="pdb-process", style={'marginLeft': 50}),

          ]),

    html.Div(id="model-build"),

    html.Br(),

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

def uploaded_pdb_files():
    """List files in the upload directory."""
    files = []
    print("Searching for files....")
    for filename in os.listdir(UPLOAD_DIRECTORY):
        #print(filename)
        path = os.path.join(UPLOAD_DIRECTORY, filename)
        if os.path.isfile(path):
            #print(path, "THIS IS A FILE")
            files.append(filename)

    print(files)

    return files

def file_download_link(filename):
    """A Plotly Dash that downloads a file"""
    location = "/app_uploaded_files/built_models/{}".format(urlquote(filename))
    return html.A(filename, href=location)



# MODEL PARSING FUNCTIONS
def check_name(uploaded_files_list):
    """ Given the uploaded data files, checks the ones that follow the naming convention"""

    pattern = r"\w+(\.\w+\.\w+)?_(\w)_\w+\.pdb+(\.gz)?"
    if len([file for file in uploaded_files_list if file.endswith(".pdb")]) == 1:
        return [file for file in uploaded_files_list if file.endswith(".pdb")]
    matches = [f for f in uploaded_files_list if re.match(pattern, f) is not None]
    if len(matches) != 0:
        return matches


def get_pdb_structures(list_paths):
    """ Given a list of pdb paths, it converts them o PDB.Bio.Structures and fetches the stoichiometry file"""

    print("Fetching PDB files...")

    list_pdb = [pdb for pdb in list_paths if pdb.endswith(".pdb")]


    if len(list_pdb) != 1:
        str_dict = {}
        i = 0
        for pdb in list_paths:
            if pdb.endswith(".pdb"):
                structure_name = (pdb.split("/"))[-1]
                parser = PDBParser(PERMISSIVE=True, QUIET=True)
                structure = parser.get_structure(structure_name, UPLOAD_DIRECTORY+pdb)
                str_dict[UPLOAD_DIRECTORY+pdb] = structure  # Path+file as KEY, structure object as VALUE
                i += 1

        return str_dict

    elif len(list_pdb) == 1:
        str_dict = inputfunctions.macrocomplex_parser(UPLOAD_DIRECTORY, interact_dist=dist_contacts, interact_contacts=num_contacts,
                                                     redundant_dist=1.9)
        return str_dict


def stoich_parser():

    stoich_dict = {}

    print("Parsing stoichiometry")
    for filename in os.listdir(UPLOAD_DIRECTORY):
        path = os.path.join(UPLOAD_DIRECTORY, filename)
        if filename.endswith(".txt"):
            stoich_dict = inputfunctions.stoichiometry_parser(path)

    if len(stoich_dict) == 0:
        return None
    else:
        return stoich_dict


# CALLBACKS
#######################################################################


@app.callback(
    Output("file-list", "children"),
    [Input("upload-data", "filename"), Input("upload-data", "contents")],
)
def update_output(uploaded_filenames, uploaded_file_contents):
    """Save uploaded files and regenerate the file list."""
    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file(name, data)
    else:
        return "No files have been uploaded"

    stoich_file = [file for file in uploaded_pdb_files() if file.endswith("txt")]
    input_files = check_name(uploaded_pdb_files())

    print("THESE ARE THE INPUT FILES:", input_files)

    try:
        complex_name = input_files[0].split("_")[0]
    except TypeError:
        complex_name = input_files[0].split(".")[0] # for a marcocomplex input

    ### CHANGE ALL THIS TO ANOTHER DASH OBJECT WITH A DIFFERENT ID, but
    ### CHAIN IT TO THE OUTPUT OF THIS CALLBACK

    if len(stoich_file) == 1:
        if len(input_files) == 1:
            return (html.P([f"The structure name is: {complex_name}", html.Br(),
                        f"It is a multicomplex input.", html.Br(),
                        f"Stoichiometry file found.", html.Br(),
                        "Begin complex building..."]))
        else:
            return (html.P([f"The structure name is: {complex_name}", html.Br(),
                        f"Contains {len(input_files)} pairs of interacting chains.", html.Br(),
                        f"Stoichiometry file found.", html.Br(),
                        "Begin complex building..."]))

    elif len(stoich_file) == 0:
        if len(input_files) > 1:
            (html.P([f"The structure name is: {complex_name}", html.Br(),
                     f"It is a multicomplex input.", html.Br(),
                     f"No stoichiometry file found.", html.Br(),
                     "Begin complex building..."]))
        else:
            return (html.P([f"The structure name is: {complex_name}", html.Br(),
                        f"Contains {len(input_files)} pairs of interacting chains.", html.Br(),
                        f"No stoichiometry file found.", html.Br(),
                        "Begin complex building..."]))


@app.callback(
    Output("parameters-out", "children"),
    [Input("dist-contact", "value"), Input("num-contact", "value")]
)
def parameter_parse(dist_contact, number_contacts):
    """ Fetch the relevant parameteres for multicomplex building """
    global dist_contacts
    global num_contacts

    dist_contacts = dist_contact
    num_contacts = number_contacts

    if dist_contacts is None and num_contacts is None:
        dist_contacts = 12
        num_contacts = 8
        return None
    elif dist_contacts is not None and num_contacts is not None:
        return f"The contact distance is set to {dist_contact} (Å). The number of contacts is set to {num_contacts}."



@app.callback(
    Output("pdb-process", "children"),
    Input("file-list", "children")
)
def pdb_processing(input_files):
    """ Given a list of PDBs, these are processed and their names returned"""


    pdb_files = get_pdb_structures(uploaded_pdb_files())
    stoich_dict = stoich_parser()
    equiv_chains, pdb_dict = inputfunctions.chain_processing(pdb_files)


    # different chains

    diff_chains = set()
    for struct in pdb_dict.values():
        for chain in struct.get_chains():
            diff_chains.add(chain)

    # similar chains

    print("Looking for similar chains...")

    ##### SIMILAR SEQUENCE CHAINS
    ### Obtain high-sequence (>95%) similarity chains from the interacting pairs
    try:
        similar_chains_prot = chainfunctions.similar_chains(pdb_dict, "Protein")  # Protein
        similar_chains_dna = chainfunctions.similar_chains(pdb_dict, "DNA")  # Nucleic acids

        # Merge both dictionaries
        similar_chains = similar_chains_prot
        similar_chains.update(similar_chains_dna)
    except IndexError:
        return None


    # unicommons
    unicommon = chainfunctions.unique_common_chains(similar_chains)
    if stoich_dict is not None:
        unicommon = chainfunctions.unicommon_completer(unicommon, stoich_dict, diff_chains)

    # stoich check:
    if stoich_dict is not None:
        stch_check = chainfunctions.stoichiometry_check_simple(unicommon, stoich_dict)
        if stch_check == True:
            pass
        else:
            print("Stoichiometry provided is incorrect")

        final_model = modelfunctions.model_construction(unicommon, pdb_dict, equiv_chains,
                                                        stoichiometry_input=stoich_dict,
                                                        forcing=True)
    else:
        final_model = modelfunctions.model_construction(unicommon, pdb_dict, equiv_chains,
                                                        forcing=False)

    if len(check_name(uploaded_pdb_files())) == 1:
        modelfunctions.save_model(final_model, uploaded_pdb_files(), outdir=FINAL_MODEL)
    else:
        modelfunctions.save_model(final_model, uploaded_pdb_files()[1:], outdir=FINAL_MODEL)


    new_model_chains = [ch for ch in final_model.get_chains()]

    print(f"NEW MODEL WITH:\n{len(new_model_chains)} CHAINS\nWHICH ARE: {new_model_chains}")
    print("Saving the structure...")

    return [html.Li(file_download_link(filename)) for filename in os.listdir(FINAL_MODEL)]

# FOR THE DASH APP TEST #
if __name__ == '__main__':
    #dt.Timer(10, dt.open_browser).start()
    #dt.app.run_server(debug=True, port=dt.port) # Hot reloading
    app.run_server(debug=True)
    print("Running the server")
