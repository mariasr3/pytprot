#############
#    DASH   #
#############
# -*- coding: utf-8 -*-
import sys
sys.path.append('/home/oth/anaconda3/lib/python3.8/site-packages')
sys.path.append('/home/oth/BHS/PYT/1.project/SBI_PYT/PPI_main') # Path to I_O_args script
sys.path.append('/home/oth/BHS/PYT/1.project/dash_test.py')
sys.path.append('/home/oth/BHS/PYT/1.project/test_project.py')

import I_O_args as args
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from urllib.parse import quote as urlquote
from flask import Flask, send_from_directory
import os
import plotly.express as px
import webbrowser
from threading import Timer
import base64
import io
import datetime




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
        dcc.Markdown("""
        ##### Introduce here the **full path** where the **input models** are stored:
        """),

        # Input box + Submit button
        html.Div([
            dcc.Input(id="input-on-submit", value='(e.g.) /home/usr/.../input_folder', type='text',
                           style={'width': '80%', 'textAlign': 'center'}),
            html.Button('Submit', id='submit-val', n_clicks=0)],
                 style={'width': '80%', 'textAlign': 'left', 'display':'inline-block'}),

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
                       'marginLeft': 0,
                       'background-color': 'white'},
                multiple=True,
                ),
        dcc.Markdown("""
        #### **File list**
        """),
        html.Ul(id="file-list"),

        # Output area
        html.Div(id="container-button-basic",
                 children='Enter your full input path',
                 style={'lineHeight': '65px', 'borderWidth': '1px',
                        'background-color':'white', 'width':'76.5%', 'textAlign': 'center'}),
    ], style={'marginTop': 60, 'marginLeft': 50}),

    # RUN COMPLEX BUILDING Button
    html.Div(
             html.Button('Run complex building', id='submit-run', n_clicks=0,
                         style={'marginTop': 20, 'marginLeft': 50})
             ),
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

def PDB_input_parser(_uploaded_files):
    pdblist = []
    pdblist.append(_uploaded_files)
    return pdblist


# CALLBACK
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

    files = uploaded_files()
    if len(files) == 0:
        return [html.Li("No files yet!")]
    else:
        return f"This will be the final model: {file_download_link(files)}"


@app.callback(
    dash.dependencies.Output('container-button-basic', 'children'),
    [dash.dependencies.Input('submit-val', 'n_clicks')],
    [dash.dependencies.State('input-on-submit', 'value')]
)
def input_path(n_clicks, value):
    """
    Checks, within a given directory, if there are any .PDB files with the correct
    naming (Pairs of interacting chains). Then, it returns them as a list.
    :param n_clicks:
    :param value:
    :return:
    """
    finalpaths = []
    if n_clicks == 0:
        return ""
    if n_clicks != 0:
        try:
            filenames = [os.listdir(value)]
            files_str = ""
            filenames = args.check_infile_names(value)
            path = str(value)
            for file in filenames:
                filename = file.split("/")[-1]
                files_str += str(filename) + ", "
                final_path = path + "/" + str(file)
                finalpaths.append(final_path)
            if files_str == "":
                return f"This directory does not contain any pairs of interacting chains"
            else:
                return f"The models employed will be:\n{files_str}\n"

        except NotADirectoryError:
            return f"{value} is not a directory"
        except FileNotFoundError:
            return f"The directory {value} does not exist"

    return finalpaths

@app.callback(
    dash.dependencies.Output('start-button', 'children'),
    [dash.dependencies.Input('submit-run', 'n_clicks')]
)
def main_code_exec(n_clicks):
    if n_clicks > 0:
        test_path = '/home/oth/BHS/PYT/1.project/test_project.py'
        with open(test_path) as infile:
            exec(infile.read())
    if n_clicks > 1:
        return f"You are already running a macrocomplex. To run another one, you must wait until the current one is finished."



