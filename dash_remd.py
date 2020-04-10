# -*- coding: utf-8 -*-
import dash
#import sys
import numpy as np
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
#sys.path.append('.')
#sys.path.append('..')
#import util
#from server import server
import random
import pandas as pd
import math
import copy
from flask import Flask
import plotly.graph_objects as go
server = Flask('my app')
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
external_scripts = ['https://yyrcd-1256568788.cos.na-siliconvalley.myqcloud.com/yyrcd/2020-03-21-iframeResizer.contentWindow.min.js']
app = dash.Dash(__name__,
                external_stylesheets=external_stylesheets,server=server)
               # external_scripts=external_scripts)
     #           (server=server,
      #          routes_pathname_prefix='/test1/')
all_frame=np.zeros((500,10,3))
all_f=pd.read_csv('std_all_frame_40.csv')

for i in range(1,501):
    for d in range(10):
        all_frame[i-1][d][0]=all_f.iloc[i][d*3]
        all_frame[i-1][d][1]=all_f.iloc[i][d*3+1]
        all_frame[i-1][d][2]=all_f.iloc[i][d*3+2]


app.layout = html.Div([#dcc.Input(id='my-id',value='initial value',type='text'),
   # dcc.Markdown(Markdown_text, dangerously_allow_html=True),
    html.Div([dcc.Graph(id='graph-with-slider')],style={'width': '80%', 'display': 'inline-block', 'padding': '0 1'}),
    #html.Div([dcc.Graph(id='population')],style={'width':'49%','display': 'inline-block'}),
    html.Div(dcc.Slider(
        id='year-slider',
        min=1,
        max=500,
        value=2,
        marks={str(m): str(m) for m in range(1,501,10)},
        step=1
    ),style={'width': '80%'})
 #   , 'padding': '0px 20px 20px 20px'})
  #  dcc.Interval(
  #      id='interval-component',
  #      interval=2 * 1000,  # in milliseconds
  #      n_intervals=0
  #  )
])

@app.callback(
    dash.dependencies.Output('graph-with-slider', 'figure'),
    [dash.dependencies.Input('year-slider', 'value')])
def update_figure(selected_number):
    #print(r_data[selected_number])
    x = all_frame[selected_number][:,0]
    y = all_frame[selected_number][:,1]
    z = all_frame[selected_number][:,2]
    data = [go.Scatter3d(x=x,y=y,z=z,text='diameter:1.5',mode='markers',marker={'size':8})]
    data.append(go.Scatter3d(x=[0,15,None,15,15,None,15,0,None,0,0,None,0,0,None,15,15,None,15,15,None,0,0,None,0,15,None,15,15,None,15,0,None,0,0],
                                y=[0,0,None,0,15,None,15,15,None,15,0,None,0,0,None,0,0,None,15,15,None,15,15,None,0,0,None,0,15,None,15,15,None,15,0],
                                z=[0,0,None,0,0,None,0,0,None,0,0,None,0,15,None,0,15,None,0,15,None,0,15,None,15,15,None,15,15,None,15,15,None,15,15],mode='lines'
                                ,line=dict(color='firebrick',width=5)))
#    N = 200
#    x = np.linspace(0, 12, N)
#    k = 1
#   # print(sec)
#    w = 1 - 0.3 * sec
#    b = 0
#    y = a * np.sin(k * x + w) + b

    return {
        'data': data,
        'layout': dict(
            xaxis={'range': [0,15]},
            yaxis={'range': [0,15]},
            zaxis={'range': [0,15]},
            margin={'l': 20, 'b': 80, 't': 10, 'r': 10},
            #height=450
        )
    }



if __name__ == '__main__':
    app.server.run(port=8000, host='127.0.0.1')
