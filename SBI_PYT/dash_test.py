#############
#    DASH   #
#############
# -*- coding: utf-8 -*-

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd
import webbrowser
from threading import Timer

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css'] # External CSS stylesheet


port = 8050

def open_browser():
    webbrowser.open_new("http://localhost:{}".format(port))

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    dcc.Markdown(children="""
    #     **pytprot**
    ## A web-app for multicomplex building
    """, style={'textAlign' : 'left', 'marginLeft': 50, 'marginTop': 80, 'marginRight': 10,
                'color':'black', 'width': '80%'}),

    dcc.Markdown(children="""Introduce here your **input** models:""",
                 style={'marginTop': 90, 'marginLeft': 50}),
    dcc.Upload(
        id='upload-data',
        children=html.Div(['Drag and Drop or ', html.A('Select PDB Files')]),
        style={'width': '50%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px',
            'marginLeft': 50,
            'background-color':'white'},
        multiple=True),
    html.Div(id='output-data-upload'),

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

