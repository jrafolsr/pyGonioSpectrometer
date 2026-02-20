#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 09:47:49 2019

@author: JOANRR
"""
#%%
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
from dash.dependencies import Input, Output, State
from time import sleep
from os.path import join as pjoin, isdir
import plotly.graph_objs as go
from instrumentation import list_spectrometers, list_ports
import numpy as np
from flask import request
from pyGonioSpectrometer import find_symmetry, GonioLogger
from gui_init import INITIAL_INTEGRATION_TIME, INITIAL_NSPECTRA,INITIAL_STEP, INITIAL_MAX_ANGLE, INITIAL_PATH, INITIAL_FILENAME, PORT
import seaborn as sns
from pathlib import Path
global gonio
import traceback

gonio = None

#LEN_WAVELENGTHS = 2028 # Just necessary if Mattias' saving scheme

TRACE_SPECTRA = [go.Scatter(x=[], y=[], name = 'counts', mode = 'lines')]
TRACE_SRI = [go.Scatter(x=[], y=[], name = 'sri', mode = 'markers',\
             xaxis = 'x2', yaxis = 'y2')]
RESET_TRACES = [go.Scatter(x=[], y=[], name = 'sri', mode = 'markers',\
             xaxis = 'x2', yaxis = 'y2'),\
                go.Scatter(x=[], y=[], name = 'counts', mode = 'lines')
                 ]


LSPECTROMETERS = list_spectrometers()
SRI = []

# Some useful funcions
def calculate_sri(wavelengths, intensity):
    wl_min = 380
    wl_max = 780
    ff = (wavelengths >= wl_min) & ( wavelengths <= wl_max)
    sri = np.trapz(intensity[ff], wavelengths[ff])
    return sri/1E6


# A function to shutdown the server
def shutdown():
    func = request.environ.get('werkzeug.server.shutdown')
    if func is None:
        raise RuntimeError('Not running with the Werkzeug Server')
    func()
    
#%%    

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

plot_data = RESET_TRACES

plot_layout = dict(margin =  {'l': 60, 'r': 60, 'b': 60, 't': 20},\
                   legend =  {'x': 0, 'y': 1, 'xanchor': 'left'},\
                   xaxis = dict(title =  "Wavelength (nm)",\
                            range = [350, 850],
                            domain=[0, 0.6]),\
                   yaxis = dict(range  =  [0, None],\
                            title =  "Counts"),
                   xaxis2 = dict(title =  "Angle (°)",\
                            # range = [-90, 90],
                            domain=[0.7, 1]),
                   yaxis2 = dict(title =  "Radiant Intensity (a.u.)",\
                            # range = [0, None],
                            anchor="x2"))

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.title = 'Goniospectrometer'
app.layout = html.Div(children =  [
        html.Div(id='live-update-text', className = 'row', children  = [
        html.Div(className = 'column left', children = [
            daq.Indicator(id='my-daq-indicator',
              value=True,
              color="#FF6633",
              size = 25, style = {'width': '50px', 'display': 'inline-block', 'vertical-align':'middle'}
              ),
        
            html.H4('Goniospectrometer', style = {'width': '40%', 'display': 'inline-block','vertical-align':'middle'}),
            html.Span([], id = 'trigger-update-graph', hidden  = True),
            dcc.Graph(id='live-update-graph', 
                figure= { "data": plot_data,
                          "layout": plot_layout
                          }
                ),
            
          html.Span('Spectrometer ID'),
            dcc.Dropdown(id  = 'dropdown-spectrometers',
                options = [{'label' : str(name), 'value': name.serial_number} for name in  LSPECTROMETERS],
                value = None if LSPECTROMETERS == [] else  LSPECTROMETERS[0].serial_number,
                placeholder = 'No detected spectrometers',
                style = {'width' : '200'},
                searchable = False
            ),
#        html.Span('Arduino COM port:'),
#        dcc.Dropdown(id  = 'dropdown-arduino',
#            options = [{'label' : name, 'value': name} for name in LPORTS],
#            value = ARDUINO_PORT if ARDUINO_PORT in LPORTS else None,
#            placeholder = 'No detected ports',
#            style = {'width' : '200'},
#            searchable = False
#            )            
            ],    
        ),
        
        html.Div(className =  'column middle', children = [
            daq.PowerButton(
               id='power-button',
               color =  "#FF5E5E",
               size = 60,
               on = False,
             ), 
            daq.StopButton(id='button-adquire',
               disabled = True,
#               title = 'Adquire a single spectra',
               buttonText = 'acquire',
               n_clicks = 0,
               ),
            daq.StopButton(id='button-set-bkg',
               disabled = True,
#               title = 'Sets the backgorund for the single spectra',
               buttonText = 'set bkg',
               n_clicks = 0,
               ),
            daq.StopButton(id='button-start',
               disabled = True,
#               title = 'Starts the goniomeasurement',
               buttonText = 'start',
               n_clicks = 0,
               ),
            daq.StopButton(id='button-update',
               disabled = True,
#               title = 'Updates the plot',
               buttonText = 'update',
               n_clicks = 0,
               ),
            daq.StopButton(id='button-clear',
               disabled = True,
#               title = 'Clears the plot',
               buttonText = 'clear',
               n_clicks = 0,
               ),            
            daq.StopButton(id='button-move-left',
               disabled = True,
#               title = 'Moves the gonio 0.1125° left',
               buttonText = 'Move CW',
               n_clicks = 0,
               ),
         daq.StopButton(id='button-move-right',
               disabled = True,
#               title = 'Moves the gonio 0.1125° right',
               buttonText = 'move CCW',
               n_clicks = 0,
               ),
         daq.StopButton(id='button-move-shutter',
               disabled = True,
#               title = 'Opens/closes the shutter',
               buttonText = 'shutter',
               n_clicks = 0,
               ),
         daq.StopButton(id='button-refresh-ports',
               disabled = False,
#               title = 'Scans for new available resources',
               buttonText = 'refresh',
               n_clicks = 0,
               ),
         daq.StopButton(id='button-autozero',
               disabled = True,
#               title = 'Automatically correctes the offset using the current  visible data',
               buttonText = 'auto-zero',
               n_clicks = 0,
               ),
         html.Span(id = 'motor-movement', children = [], hidden = True)
        ]),  
        html.Div(id  = 'div-inputs', className = 'column right', children = [
          html.Div([
            html.Div('Integration time (ms): ', id = 'label-it', style = {'display': 'inline-block'}),
            daq.PrecisionInput(
              id='integration-time',
              labelPosition = 'top',
              precision = 4,
              min = 1,
              max = 10000,
              value = INITIAL_INTEGRATION_TIME,
              style = {'display': 'inline-block'}
              )
            ]),
          html.Div(['Number of spectra: ',
              dcc.Input(
                  id = "input-n-spectra",
                  type = 'number',
                  value = INITIAL_NSPECTRA,
                  size = '5',
                  debounce = True)
              ]),
          html.Div(['Step angle (°): ',
              dcc.Input(
                  id = "input-step-angle",
                  type = 'number',
                  value = INITIAL_STEP,
                  size = '5',
                  debounce = True)
              ]),
          html.Div(['Max. angle (°): ',
              dcc.Input(
                  id = "input-max-angle",
                  type = 'number',
                  value = INITIAL_MAX_ANGLE,
                  size = '5',
                  debounce=True)
              ]),
          html.Div(['Folder: ',          
              dcc.Input(id="folder-input",
                        type="text",
                        placeholder="Folder",
                        value = INITIAL_PATH,
                        size = '40'),
              html.Span(id = 'folder-exist', children = '')
              ]),
          html.Div(['Filename: ',
              dcc.Input(id="filename-input",
                        type= "text",
                        placeholder= "Filename",
                        size = '40',
                        value = INITIAL_FILENAME)
              ]),
            ])
       ])
    ])
    
# Enable

@app.callback([Output('button-adquire', 'disabled'),
               Output('button-start', 'disabled'),
               Output('integration-time', 'disabled'),
               Output('button-move-shutter', 'disabled'),
               Output('button-clear', 'disabled'),
               Output('button-update', 'disabled'),
               Output('button-set-bkg', 'disabled'),
               Output('button-autozero', 'disabled'),
               Output('button-move-right', 'disabled'),\
               Output('button-move-left', 'disabled')],
              [Input('power-button', 'on')],
              [State('dropdown-spectrometers', 'value'),
               State('folder-input', 'value'),
               State('filename-input', 'value'),
               State('integration-time', 'value'),
               State('input-n-spectra','value'),
               ],
              prevent_initial_call = True)
def enable_buttons(on, resource_spectrometer, folder, filename, integration_time, n_spectra):
    global gonio, WAVELENGTHS
    n_buttons = 10
    
    try:
        if on:
    
            gonio = GonioLogger(filename, folder=Path(folder), integration_time = integration_time, n_spectra = n_spectra)
            
            print('INFO: Instrument is configured and ready')
            WAVELENGTHS = gonio.wavelengths
            sleep(0.250)
            
            buttons_state = False

        else:
            
            print('INFO: Instrument is off')
            gonio.shutdown()
            sleep(1)
            buttons_state = True

    except Exception as e:
        print(e)
        traceback.print_exc()
        print('ERROR: An error occured in starting the instrument')
        buttons_state = True
        
    return n_buttons * [buttons_state]  
      
# Multiple components can update everytime interval gets fired.
@app.callback(Output('live-update-graph', 'figure'),
              [Input('button-adquire', 'n_clicks'),
               Input('button-update', 'n_clicks'),
               Input('button-clear', 'n_clicks')],
              [State('live-update-graph', 'figure')],
              prevent_initial_call = True)
def update_graph(n_adq, n_upd, n_clr, figure):
    # Collect some data
    global gonio, SRI

     # Determine which button has been clicked
    ctx = dash.callback_context

    if not ctx.triggered:
        button_id = 'No clicks yet'
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]


    if button_id == 'button-adquire':
        figure['data'] = []
        
        temp = gonio.flame.get_averaged_intensities()  

        SRI.append(calculate_sri(gonio.wavelengths, temp - gonio.background))
        figure['data'].append(go.Scatter(x=list(range(len(SRI))), y=SRI, name = 'counts', mode = 'markers', xaxis = 'x2', yaxis = 'y2'))
        
        figure['data'].append(go.Scatter(x=gonio.wavelengths, y=temp, name = 'counts', mode = 'lines'))


    elif button_id == 'button-update':
        figure['data'] = []
        
        data = gonio.current_angular_scan
        angles = np.array([el[0] for el in  data])
        integrated_sri = [calculate_sri(gonio.wavelengths, el[1]) for el in data]
        
        angles_unique = np.unique(np.round(np.abs(angles), 2))
        colors = sns.color_palette('rainbow', n_colors=len(angles_unique))
        rgb_to_hex = lambda rgb: '#{:02x}{:02x}{:02x}'.format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))
        
        color_dict = dict((f'{abs(key):.1f}', rgb_to_hex(value)) for key, value in zip(angles_unique, colors))
        colors = [color_dict.get(f'{abs(angle):.1f}', 'black') for angle in angles]
        
        figure['data'].append(go.Scatter(x = angles, y = integrated_sri, name = 'counts', mode = 'markers', xaxis = 'x2', yaxis = 'y2',\
              marker=dict(color=colors, size=10 )))
        
        if len(data):
            for angle, sri in gonio.current_angular_scan:
                key = f'{abs(angle):.1f}'
                figure['data'].append(go.Scatter(x = gonio.wavelengths, y = sri, name = f'{angle:.1f}', mode = 'lines', line=dict(color=color_dict.get(key, 'black'))))
        else:
            figure['data'].append(go.Scatter(x=[], y=[], name = 'counts', mode = 'lines'))

        
    elif button_id == 'button-clear':
        print('INFO: Clearing the plot')
        SRI = []
        gonio.current_angular_scan = []
        figure['data'] = [go.Scatter(x=[], y=[], name = 'counts', mode = 'lines'),\
                          go.Scatter(x = [], y = [], name = 'intensity', mode = 'markers', xaxis = 'x2', yaxis = 'y2')]
    else:
        pass
        
    return figure

@app.callback(Output('trigger-update-graph', 'children'),
              [Input('button-start', 'n_clicks')],
              [State('folder-input', 'value'),
               State('filename-input', 'value'),
               State('input-max-angle','value'),
               State('input-step-angle','value'),
               State('integration-time','value'),
               State('input-n-spectra','value')],
              prevent_initial_call = True)
def run_measurement(n, folder, filename, angle_max, angle_step, integration_time, n_spectra):
    global gonio
    gonio.filename = filename
    gonio.folder = Path(folder)
    gonio.angle_max, gonio.angle_step, gonio.integration_time, gonio.n_spectra = angle_max, angle_step, integration_time, n_spectra
    gonio.take_dark_spectra()
    gonio.take_gonio_measurement(suffix = '', plot=False)
    
    print('INFO: Measurement DONE!')
    
    return ' '

@app.callback(Output('folder-exist', 'children'),
              [Input('folder-input', 'value')])
def check_folder(value):

    if not isdir(pjoin(value)):
        msg = 'ERROR: Folder does not exist'
        print(msg)
    else:
        msg = ''

    return msg

@app.callback([Output('dropdown-spectrometers', 'options'),
                Output('dropdown-spectrometers', 'value')],
              [Input('button-refresh-ports', 'n_clicks')],
              prevent_initial_call = True)
def refresh_ports(n_ports):
    global LSPECTROMETERS,LPORTS
    LSPECTROMETERS = list_spectrometers()
    LPORTS = list_ports()
    
#    options_arduino = [{'label' : name, 'value': name} for name in LPORTS]
    options_spec = [{'label' : str(name), 'value': name.serial_number} for name in  LSPECTROMETERS]
    
    value_spec = None if LSPECTROMETERS == [] else  LSPECTROMETERS[0].serial_number
#    value_arduino = ARDUINO_PORT if ARDUINO_PORT in LPORTS else None

    return options_spec, value_spec


@app.callback(Output('label-it', 'children'),
              [Input('integration-time', 'value'), 
               Input('input-n-spectra', 'value')],
              prevent_initial_call = True)
def set_integration_time(integration_time, n_spectra):
    
    #gonio.integration_time, gonio.n_spectra = integration_time, n_spectra
    gonio.flame.config(integration_time, n_spectra)   
    
    print(f'INFO: Integration time set to {integration_time:.4g} ms x N = {n_spectra: 3d}')

    return f'Integration time is {integration_time:.4g} ms x N = {n_spectra: 3d}'

# Obs! Potnential infinite loops, but seems o work for some reason....
@app.callback(Output('input-step-angle','value'),
              Input('input-step-angle','value'), # Use State to read the current value
              prevent_initial_call=True)
def set_step_angle(step):
    if step is None:
        return step
    
    base = 2.7
    new_step = round(round(step/base)*base,2)
    print(f'INFO: Setting step to {new_step:.2f} deg ')
    return new_step


@app.callback(Output('input-max-angle','value'),
               Input('input-max-angle','value'),
               State('input-step-angle','value'),
              prevent_initial_call = True)
def set_max_angle(max_angle, base):
    if max_angle is None:
        return max_angle
    
    n = round(max_angle/base)
    new_max_angle = round(n*base,2)
    while new_max_angle > 90:
        n -=1
        new_max_angle = round(n*base,2)
        if new_max_angle <= base:
            break
    print(f'INFO: Setting max angle to {new_max_angle:.2f} deg ')
    return new_max_angle


@app.callback(Output('motor-movement', 'children'),
              [Input('button-move-shutter', 'n_clicks'),
              Input('button-set-bkg', 'n_clicks'),
              Input('button-autozero', 'n_clicks'),
              Input('button-move-right', 'n_clicks'),\
              Input('button-move-left', 'n_clicks')],
               [State('input-step-angle','value')])
def gonio_and_spectra_functions(nshutter, nbkg, nautozero, move_left, move_right, step_angle):
    global gonio
    # Determine which button has been clicked
    ctx = dash.callback_context

    if not ctx.triggered:
        button_id = 'No clicks yet'
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if button_id == 'button-move-shutter':
        if gonio.gonio.shutter_is_closed:
            gonio.gonio.open_shutter()
        else:
            gonio.gonio.close_shutter()
            
    elif button_id == 'button-set-bkg':
        print('INFO: Background spectra set')
        gonio.background = gonio.flame.get_averaged_intensities()   
    elif button_id == 'button-autozero':
        data = gonio.current_angular_scan

        if len(data)> 3:
            angles = np.array([el[0] for el in  data])
            integrated_sri = [calculate_sri(gonio.wavelengths, el[1]) for el in data]
            x0 = find_symmetry(angles, integrated_sri)
            offset_angle = -x0
            
            out_angle = gonio.gonio.move_angle(np.round(offset_angle, 4), correct_drift = False)
            gonio.gonio.steps_counter = 0
            print(f'INFO: Zero offset is {offset_angle:.4f}° and moved by {out_angle:.4f}')
        else: 
            print(f'ERROR: No or too few data to fit')
    elif button_id == 'button-move-right':
        gonio.gonio.move_angle(np.round(step_angle, 4), correct_drift = False)
    elif button_id == 'button-move-left':
        gonio.gonio.move_angle(-np.round(step_angle, 4), correct_drift = False)    
    else:
        pass
    return
        

if __name__ == '__main__':
    try:
        app.run_server(debug = True, port = PORT)
    except KeyboardInterrupt as e:
        print(e)
    finally:  
        if gonio is not None:
            gonio.shutdown()