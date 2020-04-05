# -*- coding: utf-8 -*-
import dash
#import sys
import numpy as np
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
#sys.path.append('.')
#sys.path.append('..')
import util
#from server import server
import random
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


Markdown_text = r"""
## Describing the Free Particle
We start by describing a free particle: a particle that is not under the influence of a potential.   
 As any other particle, the state of a ***Free Particle*** is described with a ket $\left|\psi(x)\right>$. In order to learn about this particle (measure its properties) we must construct *Hermitian Operators*. For example, what is the **momentum operator $\hat P$**?
 
The momentum operator most obey the following property (eigenfunction/eigenvalue equation):
$$\hat P \left| \psi_k(x) \right> =p\left | \psi_k(x)\right>  \tag{1}$$ 
where *p* is an eigenvalue of a *Hermitian operator* and therefore it is a real number.
In the $x$ representation, using the momentum operator as $\hat P =-i\hbar \frac{\partial }{\partial x}$, we can solve equation 1 by proposing a function to represent $\left| \psi_k(x) \right>$ as $\psi_k(x) = c\ e^{ikx}$, where $k$ is a real number.
Let's see if it works:  
$$\hat P \psi_k(x) =p \psi_k(x)$$ 
$$-i\hbar \frac{\partial {c\ e^{ikx}}}{\partial x} =-i\hbar\ c\ ik\ e^{ikx} $$ 
$$\hbar k\ c\ e^{ikx} = \hbar k\ \psi_k(x) \tag{2}$$
with $p=\hbar k$
-------------------
"""

sigma   = 1
epsilon = 2
L=3*sigma
delta=0.1
kb=1.38e-23*6.022e23/(1000*4.184)
#start = time.time()
def cutoff(r):
    if r > L:
        r -= L
    elif r < 0:
        r += L
    return r

def energy(r,epsilon):
    energy = 4*epsilon*(((sigma/r)**12)-((sigma/r)**6))
    return energy

def move(diff_e,e,r,temp):
    if diff_e < 0:
        
        e += diff_e
    else:
        rand = random.random()
        if math.exp(-diff_e/(kb*temp)) > rand:
            
            e += diff_e
        else:
            r = origi
            e = pre_e
    return r,e
#initialize system

temp1=[50,100,200,300,500,1000]
#initialize system
steps=20000
r = random.uniform(0,L)
r_data=np.zeros((6,steps))
e_data=np.zeros((6,steps))
#run mc with different temperature
acc_ratio=[]
for i in range(6):
    acc_step = 0
    e = energy(r,epsilon)
    for step in range(0,steps):
        origi = copy.deepcopy(r)
        pre_e = e
        r += random.uniform(-1,1)*delta
        r = cutoff(r)
        new_e = energy(r,epsilon)
        diff_e = new_e - pre_e
        r,e = move(diff_e,e,r,temp1[i])
        if r != origi:
            acc_step += 1
        r_data[i][step]=r
        e_data[i][step]=e
    acc_ratio.append(round(acc_step/steps,2))

#acc_ratio = acc_step/steps
#print(r_data)
Markdown_text = util.convert_latex(Markdown_text)
app.layout = html.Div([#dcc.Input(id='my-id',value='initial value',type='text'),
   # dcc.Markdown(Markdown_text, dangerously_allow_html=True),
    html.Div([dcc.Dropdown(id='temp',options=[{'label':i,'value':i} for i in temp1],value=temp1[3])],style={'width': '49%', 'display': 'inline-block'}),
    html.Div([dcc.Graph(id='graph-with-slider')],style={'width': '49%', 'display': 'inline-block', 'padding': '0 1'}),
    html.Div([dcc.Graph(id='population')],style={'width':'49%','display': 'inline-block'}),
    html.Div(dcc.Slider(
        id='year-slider',
        min=0,
        max=step+1,
        value=0,
        marks={str(m): str(m) for m in range(0,step+1,2000)},
        step=100
    ),style={'width': '45%', 'padding': '0px 20px 20px 20px'})
  #  dcc.Interval(
  #      id='interval-component',
  #      interval=2 * 1000,  # in milliseconds
  #      n_intervals=0
  #  )
])


@app.callback(
    dash.dependencies.Output('graph-with-slider', 'figure'),
    [dash.dependencies.Input('temp', 'value'),
     dash.dependencies.Input('year-slider', 'value')])
def update_figure(selected_temp,selected_step):
    #print(r_data[selected_number])
    x = [r_data[temp1.index(selected_temp)][selected_step]]
    y = [e_data[temp1.index(selected_temp)][selected_step]]
    data = [dict(x=x,y=y,text='nmd',mode='markers',marker={'size':20})]
    #print(data)
    x_1 = np.arange(0,10,0.001)
    y_1 = [energy(i,epsilon) for i in x_1]
    #print(y_1)
    data.append(dict(x=x_1,y=y_1,mode='lines',line={'color':'firebrick','width':4}))
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
            xaxis={'range': [0, 5]},
            yaxis={'range': [-3,3]},
            margin={'l': 40, 'b': 30, 't': 10, 'r': 0},
            height=450
        )
    }
@app.callback(
    dash.dependencies.Output('population','figure'),
    [dash.dependencies.Input('temp', 'value'),
     dash.dependencies.Input('year-slider','value')]
)
def update_p_figure(selected_temp,selected_step):
    x = r_data[temp1.index(selected_temp)][0:selected_step]
   
    data = [go.Histogram(x=x,nbinsx=1000, histnorm='percent')]
    #for i in range(len(temp1)-1):
    #   if i != temp1.index(selected_temp):
    #        data.append(go.Histogram(x=r_data[i],nbinsx=1000,histnorm='probability density',opacity=0.15*(i+1)))
    x_2 = np.arange(0,10,0.001)
    y_2 = [energy(i,epsilon) for i in x_2]
    data.append(dict(x=x_2,y=y_2,mode='lines',line={'color':'firebrick','width':4}))
    return {
        'data': data,
        'layout': dict(
            xaxis={'range': [0, 5]},
            yaxis={'range': [-2,4]},
            height= 450,
            margin={'l': 20, 'b': 30, 'r': 10, 't': 10}
        )
    }


if __name__ == '__main__':
    app.run_server(debug=True)
