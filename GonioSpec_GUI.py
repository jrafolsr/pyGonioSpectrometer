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
from time import sleep, time
from os.path import join as pjoin, isdir
import plotly.graph_objs as go
from instrumentation import SpectraMeasurement, list_spectrometers, list_ports, ArduinoMotorController
import numpy as np
from winsound import Beep
from datetime import datetime
from scipy.optimize import curve_fit
import getopt
import sys

from GonioSpec_init import SATURATION_COUNTS, INITIAL_INTEGRATION_TIME, INITIAL_NSPECTRA, ARDUINO_PORT,INITIAL_STEP, INITIAL_MAX_ANGLE, INITIAL_PATH, INITIAL_FILENAME, WAIT_TIME


flame = None
gonio = None
WAVELENGTHS = None
INTENSITIES = None

#LEN_WAVELENGTHS = 2028 # Just necessary if Mattias' saving scheme

TRACE_SPECTRA = [go.Scatter(x=[], y=[], name = 'counts', mode = 'lines')]
TRACE_SRI = [go.Scatter(x=[], y=[], name = 'sri', mode = 'markers',\
             xaxis = 'x2', yaxis = 'y2')]
RESET_TRACES = [go.Scatter(x=[], y=[], name = 'counts', mode = 'lines'),
                 go.Scatter(x=[], y=[], name = 'sri', mode = 'markers',\
             xaxis = 'x2', yaxis = 'y2')]
    
TRACES = [go.Scatter(x=[], y=[], name = 'counts', mode = 'lines'),
                 go.Scatter(x=[], y=[], name = 'sri', mode = 'markers',\
             xaxis = 'x2', yaxis = 'y2')]


LSPECTROMETERS = list_spectrometers()
LPORTS  = list_ports()
CURRENT_ANGLE = []
SRI = []

# Some useful funcions
def calculate_sri(wavelengths, intensity):
    wl_min = 380
    wl_max = 780
    ff = (wavelengths >= wl_min) & ( wavelengths <= wl_max)
    sri = np.trapz(intensity[ff], wavelengths[ff])
    return sri/1E6

# A second order symmetrical polynomial
pol2_sym = lambda x,a,c, x0: a*(x-x0)**2 + c

def write_to_file(etime, angle, data, file, debug = False):
    t  = np.hstack((etime, angle, data))
    t = t.reshape((1, t.shape[0]))
    with open(file, 'a') as f:
        np.savetxt(f, t, fmt = '% 8.2f')
    if debug:
        print(f'INFO: Data saved at \n\t{file:s}')

#%%    

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

plot_data = TRACES

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
        html.Span('Arduino COM port:'),
        dcc.Dropdown(id  = 'dropdown-arduino',
            options = [{'label' : name, 'value': name} for name in LPORTS],
            value = ARDUINO_PORT if ARDUINO_PORT in LPORTS else None,
            placeholder = 'No detected ports',
            style = {'width' : '200'},
            searchable = False
            )            
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
               buttonText = 'left',
               n_clicks = 0,
               ),
         daq.StopButton(id='button-move-right',
               disabled = True,
#               title = 'Moves the gonio 0.1125° right',
               buttonText = 'right',
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
                  size = '5')
              ]),
          html.Div(['Step angle (°): ',
              dcc.Input(
                  id = "input-step-angle",
                  type = 'number',
                  value = INITIAL_STEP,
                  size = '5')
              ]),
          html.Div(['Max. angle (°): ',
              dcc.Input(
                  id = "input-max-angle",
                  type = 'number',
                  value = INITIAL_MAX_ANGLE,
                  size = '5')
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
               Output('button-move-left', 'disabled'),
               Output('button-move-right', 'disabled'),
               Output('button-move-shutter', 'disabled'),
               Output('button-clear', 'disabled'),
               Output('button-update', 'disabled'),
               Output('button-set-bkg', 'disabled'),
               Output('button-autozero', 'disabled')],
              [Input('power-button', 'on')],
              [State('dropdown-spectrometers', 'value'),
               State('dropdown-arduino', 'value'),
               State('integration-time', 'value')],
              prevent_initial_call = True)
def enable_buttons(on, resource_spectrometer,resource_gonio, integration_time):
    global gonio, flame, WAVELENGTHS
    n_buttons = 10
    
    try:
        if on:

            flame = SpectraMeasurement(resource_spectrometer, integration_time)
            flame.open()
            
            gonio = ArduinoMotorController(resource_gonio)
            
            print('INFO: Instrument is configured and ready')
            WAVELENGTHS = flame.get_wavelengths()
            sleep(0.250)
            
            buttons_state = False

        else:
            
            print('INFO: Instrument is off')
            if flame is not None:
                flame.close()   
            if gonio is not None:
                gonio.close()
            sleep(1)
            buttons_state = True

    except Exception as e:
        print(e)
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
    global flame, gonio, WAVELENGTHS, INTENSITIES, TRACES, SRI, CURRENT_ANGLE, RESET_TRACES

     # Determine which button has been clicked
    ctx = dash.callback_context

    if not ctx.triggered:
        button_id = 'No clicks yet'
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    
    
    if button_id == 'button-adquire':
        temp = flame.get_intensities()
        
        if flame.background is not None:
            SRI.append(calculate_sri(WAVELENGTHS, temp -  flame.background))
            figure['data'][1] = go.Scatter(x=list(range(len(SRI))), y=SRI, name = 'counts', mode = 'markers', xaxis = 'x2', yaxis = 'y2')

        figure['data'][0] = go.Scatter(x=WAVELENGTHS, y=temp, name = 'counts', mode = 'lines')    


    elif button_id == 'button-update':
                
        figure['data'] = TRACES
        
        figure['data'][1] = go.Scatter(x = CURRENT_ANGLE, y = SRI, name = 'counts', mode = 'markers', xaxis = 'x2', yaxis = 'y2')
        
    elif button_id == 'button-clear':
        print('INFO: Clearing the plot')
        SRI = []
        CURRENT_ANGLE = []
        RESET_TRACES = [go.Scatter(x=[], y=[], name = 'counts', mode = 'lines'),
                 go.Scatter(x=[], y=[], name = 'sri', mode = 'markers',\
             xaxis = 'x2', yaxis = 'y2')]
        TRACES = RESET_TRACES
        figure['data'] = RESET_TRACES
    else:
        pass
        
    return figure

@app.callback(Output('trigger-update-graph', 'children'),
              [Input('button-start', 'n_clicks')],
              [State('folder-input', 'value'),
               State('filename-input', 'value'),
               State('input-n-spectra','value'),
               State('input-max-angle','value'),
               State('input-step-angle','value'),
               State('integration-time','value')],
              prevent_initial_call = True)
def run_measurement(n, folder, filename, Nspectra, angle_max, angle_step, int_time):
    global WAVELENGTHS, INTENSITIES, WAIT_TIME, TRACES, SRI, CURRENT_ANGLE,RESET_TRACES
#    global LEN_WAVELENGTHS # Just if using Mattias' scheme
    TRACES = RESET_TRACES
    
    SRI = []
    CURRENT_ANGLE = []
    n_angles = int(angle_max*100) // int(angle_step*100) + 1
    # n_columns = n_angles * 2 - 1 + 4
    n_steps = 2 * (n_angles -1)

    # Open the port again, just in case, and configure with the current integration time
    flame.open()
    flame.config(int_time, n_spectra = Nspectra)
    
     # Timestamps for the header and filename
    itimestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
    timestamp = datetime.now().strftime("%Y-%m-%dT%Hh%Mm%Ss")
    start_time = time()
    
    path = pjoin(folder, timestamp + '_' + filename + '.dat')
    
    with open(path, 'a') as f:
        f.write(itimestamp + ' # Timestamp at the beginning of the measurement\n')
        f.write(f'{int_time:.0f} # Integration time in (ms)\n')
        f.write(f'{Nspectra:d} # Number of spectra taken\n')
              
    
    print('INFO: Measurement STARTED!')
    Beep(3000, 250)
    # Getting the wavelength vector
    WAVELENGTHS = flame.get_wavelengths()
#    data[:,0] = WAVELENGTHS
    # Saving the data in the new scheme
    write_to_file(time() - start_time, np.nan, WAVELENGTHS, path)
    
    # Take the dark spectra at zero, assuming shutter closed
    print('\tINFO: Taking dark spectra')
    temp = flame.get_averaged_intensities()
    flame.set_background(temp)
    TRACES[0] = go.Scatter(x = WAVELENGTHS, y = temp, name = 'dark', mode = 'lines')
    # Saving the data in the new scheme
    write_to_file(time()-start_time, np.nan, temp, path)       
                  
    # Open the shutter
    gonio.move_shutter()
    [Beep(2000,200) for i in range (5)] # Reminder of measurement starting
    sleep(WAIT_TIME)
    
    # Take 1st spectra at zero
    print('\tINFO: Taking spectra at 0°')
    temp = flame.get_averaged_intensities()
    TRACES.append(go.Scatter(x = WAVELENGTHS, y = temp, name = '0°', mode = 'lines'))
    SRI.append(calculate_sri(WAVELENGTHS, temp - flame.background)) 
    CURRENT_ANGLE.append(0.0)
    # Saving the data in the new scheme
    write_to_file(time()-start_time, 0.0, temp, path)    
    
    # Starting measurement. First, move to the last position
    out_angle = gonio.move_angle(round(-1 * angle_max, 4))
    
    sleep(WAIT_TIME*3) # Wait long enough for the movement to finish
    
    #Initialize the error made
    error = round(-1.0 * (angle_max - out_angle),4)
    total = 0
    current_angle = -out_angle

    k = 0 
    for k in range(n_steps):
        # Save the angle
#        first_row[0, k + 3] = current_angle
        print(f'\tINFO: Taking spectra at {current_angle:.1f}°')
        temp = flame.get_averaged_intensities()
        #Warning in case of saturation
        if np.any(temp > SATURATION_COUNTS): print('WARNING: Some values are saturating...')
        # Saving the data in the new scheme
        write_to_file(time()-start_time, current_angle, temp, path)
        
        # Plotting globals
        TRACES.append(go.Scatter(x = WAVELENGTHS, y = temp, name = f'{current_angle:.0f}°', mode = 'lines')
              )
        SRI.append(calculate_sri(WAVELENGTHS, temp - flame.background))
        CURRENT_ANGLE.append(current_angle)
        
        # Calculating next step
        next_step = angle_step + error
    # print(f'Step = {angle_step:.0f}, Out_angle = {out_angle:.2f}, Corrected_angle = {next_step:.2f}, Error made: {error:.2}')
        out_angle = gonio.move_angle(next_step)  
        error = (next_step - out_angle)
        total += abs(out_angle)
        current_angle += out_angle
        sleep(WAIT_TIME)
        
    # Take last angle spectra

    print(f'\tINFO: Taking spectra at {current_angle:.1f}°')
    temp = flame.get_averaged_intensities()
    # Saving the data in the new scheme
    write_to_file(time()-start_time, current_angle, temp, path)
    
    # Plotting globals
    TRACES.append(go.Scatter(x = WAVELENGTHS, y = temp, name = f'{current_angle:.0f}°', mode = 'lines'))
    SRI.append(calculate_sri(WAVELENGTHS, temp - flame.background))
    CURRENT_ANGLE.append(current_angle)

    
    # Going back to initial angle
    back_angle = round(-1 * abs(current_angle), 4)
    out_angle = gonio.move_angle(back_angle)
    current_angle -= out_angle

    
    sleep(WAIT_TIME * 3)
    
    # Taking last spectra at zero
    print('\tINFO: Taking last spectra at  0°')
    temp = flame.get_averaged_intensities()
    # Saving the data in the new scheme
    write_to_file(time()-start_time, current_angle, temp, path)
    
    # Plotting globals    
    TRACES.append(go.Scatter(x = WAVELENGTHS, y = temp, name = '0°', mode = 'lines'))
    SRI.append(calculate_sri(WAVELENGTHS, temp - flame.background))
    CURRENT_ANGLE.append(current_angle)

    # Close shutter
    gonio.move_shutter()
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
               Output('dropdown-arduino', 'options'),
                Output('dropdown-spectrometers', 'value'),
               Output('dropdown-arduino', 'value')],
              [Input('button-refresh-ports', 'n_clicks')],
              prevent_initial_call = True)
def refresh_ports(n_ports):
    global LSPECTROMETERS,LPORTS
    LSPECTROMETERS = list_spectrometers()
    LPORTS = list_ports()
    
    options_arduino = [{'label' : name, 'value': name} for name in LPORTS]
    options_spec = [{'label' : str(name), 'value': name.serial_number} for name in  LSPECTROMETERS]
    
    value_spec = None if LSPECTROMETERS == [] else  LSPECTROMETERS[0].serial_number
    value_arduino = ARDUINO_PORT if ARDUINO_PORT in LPORTS else None

    return options_spec, options_arduino, value_spec, value_arduino


@app.callback(Output('label-it', 'children'),
              [Input('integration-time', 'value')],
              prevent_initial_call = True)
def set_integration_time(value):
    
    if flame is not None:
        flame.open()
        flame.config(value)
        
    print(f'INFO: Integration time set to {value:.4g} ms')

    return f'Integration time is {value:.4g} ms'

@app.callback(Output('motor-movement', 'children'),
              [Input('button-move-left', 'n_clicks'),
              Input('button-move-right', 'n_clicks'),
              Input('button-move-shutter', 'n_clicks'),
              Input('button-set-bkg', 'n_clicks'),
              Input('button-autozero', 'n_clicks')],
              [State('input-step-angle','value')])
def gonio_and_spectra_functions(nleft, nright, nshutter, nbkg, nautozero, angle_step):
    global gonio, traces
    # Determine which button has been clicked
    ctx = dash.callback_context

    if not ctx.triggered:
        button_id = 'No clicks yet'
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if button_id == 'button-move-left':
        tangle = round(angle_step,4)
        print(f'INFO: Moving {tangle:.2f}° counter-clockwise.')
        gonio.move_angle(tangle)
    elif button_id == 'button-move-right':
        tangle = round(-1.0 * angle_step,4)
        print(f'INFO: Moving {tangle:.2f}° clockwise.')
        gonio.move_angle(tangle)
    elif button_id == 'button-move-shutter':
        gonio.move_shutter()
    elif button_id == 'button-set-bkg':
        print('INFO: Background spectra set')
        flame.set_background(flame.get_intensities())
    elif button_id == 'button-autozero':
        
        if CURRENT_ANGLE  != []:
            x = np.array(CURRENT_ANGLE)
            y = np.array(SRI)
            ymax = y.max()
            y /= ymax
            popt, _ = curve_fit(pol2_sym,x,y,p0 = [-1e-6,1,0])
            x0 = popt[2]
            x1 = np.linspace(x.min(),x.max(),101)
            y1 = pol2_sym(x1, *popt) * ymax
            TRACES.append(go.Scatter(x = x1, y = y1, name = 'fit', mode = 'lines', xaxis = 'x2', yaxis = 'y2'))
            gonio.move_angle(x0)
            print(f'INFO: Zero offset is {x0:.2f}°')
        else: 
            print('ERROR: No data to fit')
        
    else:
        pass
    return
        

if __name__ == '__main__':
   # Default values
    debug = True
    port = 8051
    user_reloader = False
    argv = sys.argv[1:]
    
    try:
        options, args = getopt.getopt(argv, "p:d:r:",
                                   ["port =",
                                    "debug =",
                                    "user_reloader = "])
        
        
        for name, value in options:
            if name in ['-d', '--debug']:
                if value.lower() in ['true', '1']:
                    debug = True
                else:
                    debug = False       
            elif name in ['-p', '--port']:
                port = value
            elif name in ['-r', '--user_reloader']:
                if value.lower() in ['true', '1']:
                    user_reloader = True
                else:
                    user_reloader = False
      
        app.run_server(debug = debug, port = port, use_reloader = user_reloader)
    
    except KeyboardInterrupt:
        print("Program terminated.")
    except Exception as e:
        print(e)
   
    if flame is not None:
        flame.close()   
    if gonio is not None:
        gonio.close()