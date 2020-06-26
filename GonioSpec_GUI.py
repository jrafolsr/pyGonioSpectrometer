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



N_CLICK_PREVIOUS = 0
flame = None
gonio = None
WAVELENGTHS = None
INTENSITIES = None
WAIT_TIME = 3
LEN_WAVELENGTHS = 2028
TRACE_SPECTRA = [go.Scatter(x=[], y=[], name = 'counts', mode = 'lines')]
TRACE_SRI = [go.Scatter(x=[], y=[], name = 'sri', mode = 'markers',\
             xaxis = 'x2', yaxis = 'y2')]
RESET_TRACES = [go.Scatter(x=[], y=[], name = 'counts', mode = 'lines'),
                 go.Scatter(x=[], y=[], name = 'sri', mode = 'markers',\
             xaxis = 'x2', yaxis = 'y2')]
    
TRACES = [go.Scatter(x=[], y=[], name = 'counts', mode = 'lines'),
                 go.Scatter(x=[], y=[], name = 'sri', mode = 'markers',\
             xaxis = 'x2', yaxis = 'y2')]

SATURATION_COUNTS = 65535
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

def write_to_file(etime, angle, data, file):
    t  = np.hstack((etime, angle, data))
    t = t.reshape((1, t.shape[0]))
    with open(file, 'a') as f:
       np.savetxt(f, t, fmt = '% 8.2f')
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
                   yaxis2 = dict(title =  "Spectral Radiant Intensity (a.u.)",\
                            # range = [0, None],
                            anchor="x2"))

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
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
                value = 'No devices available' if LSPECTROMETERS == [] else  LSPECTROMETERS[0].serial_number,
                style = {'width' : '200'},
                searchable = False
            ),
        html.Span('Arduino COM port:'),
        dcc.Dropdown(id  = 'dropdown-arduino',
            options = [{'label' : name, 'value': name} for name in LPORTS],
            value = 'No detected ports' if LPORTS == [] else LPORTS[0],
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
               buttonText = 'adquire',
               n_clicks = 0,
               ),
            daq.StopButton(id='button-set-bkg',
               disabled = True,
               buttonText = 'set bkg',
               n_clicks = 0,
               ),
            daq.StopButton(id='button-start',
               disabled = True,
               buttonText = 'start',
               n_clicks = 0,
               ),
            daq.StopButton(id='button-update',
               disabled = True,
               buttonText = 'update',
               n_clicks = 0,
               ),
            daq.StopButton(id='button-clear',
               disabled = True,
               buttonText = 'clear',
               n_clicks = 0,
               ),            
            daq.StopButton(id='button-move-left',
               disabled = True,
               buttonText = 'left',
               n_clicks = 0,
               ),
         daq.StopButton(id='button-move-right',
               disabled = True,
               buttonText = 'right',
               n_clicks = 0,
               ),
         daq.StopButton(id='button-move-shutter',
               disabled = True,
               buttonText = 'shutter',
               n_clicks = 0,
               ),
         daq.StopButton(id='button-refresh-ports',
               disabled = False,
               buttonText = 'refresh',
               n_clicks = 0,
               ),
         daq.StopButton(id='button-autozero',
               disabled = True,
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
              value = 1.000E2,
              style = {'display': 'inline-block'}
              )
            ]),
          html.Div(['Number of spectra: ',
              dcc.Input(
                  id = "input-n-spectra",
                  type = 'number',
                  value = 10,
                  size = '5')
              ]),
          html.Div(['Step angle (°): ',
              dcc.Input(
                  id = "input-step-angle",
                  type = 'number',
                  value = 10,
                  size = '5')
              ]),
          html.Div(['Max. angle (°): ',
              dcc.Input(
                  id = "input-max-angle",
                  type = 'number',
                  value = 80,
                  size = '5')
              ]),
          html.Div(['Folder: ',          
              dcc.Input(id="folder-input",
                        type="text",
                        placeholder="Folder",
                        value = r'C:\Users\OPEGLAB\Documents\data\goniospectrometer',
                        size = '40'),
              html.Span(id = 'folder-exist', children = '')
              ]),
          html.Div(['Filename: ',
              dcc.Input(id="filename-input",
                        type="text",
                        placeholder="Filename",
                        size = '40',
                        value = 'angular_sri')
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
#            gonio.disable_gonio()
            gonio.close()
            flame.close()
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
            SRI.append(calculate_sri(WAVELENGTHS, temp))
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
    global WAVELENGTHS, INTENSITIES, WAIT_TIME, LEN_WAVELENGTHS, TRACES, SRI, CURRENT_ANGLE,RESET_TRACES
    TRACES = RESET_TRACES
    
    SRI = []
    CURRENT_ANGLE = []
    n_angles = int(angle_max*100) // int(angle_step*100) + 1
    n_columns = n_angles * 2 - 1 + 4
    n_steps = 2 * (n_angles -1)
    
    # Data saving according to Mattias strategy (to be improved)
#    first_row = np.zeros((1, n_columns)) # INT_TIME, N_AV + ANGLES
#    first_row[0,0] = int_time
#    first_row[0,1] = Nspectra
    # First columns will be the wavelengths, second the dark spectra, the rest the angles
#    data = np.zeros((LEN_WAVELENGTHS, n_columns)) # 2028 is the length of the output vector
    
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
              
    
    
    Beep(3000, 250)
    # Getting the wavelength vector
    WAVELENGTHS = flame.get_wavelengths()
#    data[:,0] = WAVELENGTHS
    # Saving the data in the new scheme
    write_to_file(time() - start_time, np.nan, WAVELENGTHS, path)
    
    # Take the dark spectra at zero, assuming shutter closed
    print('INFO: Taking dark spectra....')
    temp = flame.get_averaged_intensities()
    flame.set_background(temp)
    TRACES[0] = go.Scatter(x = WAVELENGTHS, y = temp, name = 'dark', mode = 'lines')
    # Saving the data in the new scheme
    write_to_file(time()-start_time, np.nan, temp, path)
    # Saving the data in Mattias's scheme   
#    data[:,1] = temp            
                  
    # Open the shutter
    gonio.move_shutter()
    Beep(3000, WAIT_TIME)
    
    # Take 1st spectra at zero
    print('INFO: Taking 0° spectra....')
    temp = flame.get_averaged_intensities()
    TRACES.append(go.Scatter(x = WAVELENGTHS, y = temp, name = '0°', mode = 'lines'))
    SRI.append(calculate_sri(WAVELENGTHS, temp - flame.background)) 
    CURRENT_ANGLE.append(0.0)
    # Saving the data in the new scheme
    write_to_file(time()-start_time, 0.0, temp, path)    
    # Saving the data in Mattias's scheme   
#    data[:,2] = temp   
#    first_row[0,2] = 0.0
    
    # Starting measurement. First, move to the last position
    out_angle = gonio.move_angle(-1 * angle_max)
    sleep(WAIT_TIME*2) # Wait long enough for the movement to finish
    
    #Initialize the error made
    error = (angle_max - out_angle)
    total = 0
    current_angle = -out_angle

    k = 0 
    for k in range(n_steps):
        # Save the angle
#        first_row[0, k + 3] = current_angle
        print('INFO: Taking spectra....')
        temp = flame.get_averaged_intensities()
        #Warning in case of saturation
        if np.any(temp > SATURATION_COUNTS): print('WARNING: Some values are saturating...')
        # Saving the data in the new scheme
        write_to_file(time()-start_time, current_angle, temp, path)
        # Saving the data in Mattias's scheme   
#        data[:, k + 3] = temp
        
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

    print('INFO: Taking spectra....')
    temp = flame.get_averaged_intensities()
    # Saving the data in the new scheme
    write_to_file(time()-start_time, current_angle, temp, path)
    # Saving the data in Mattias's scheme       
#    first_row[0,k + 4] = current_angle
#    data[:,k + 4] = temp
    
    # Plotting globals
    TRACES.append(go.Scatter(x = WAVELENGTHS, y = temp, name = f'{current_angle:.0f}°', mode = 'lines'))
    SRI.append(calculate_sri(WAVELENGTHS, temp - flame.background))
    CURRENT_ANGLE.append(current_angle)

    
    # Going back to initial angle
    back_angle = -1.0 * abs(current_angle)
    out_angle = gonio.move_angle(back_angle)
    current_angle -= out_angle

    
    sleep(WAIT_TIME)
    
    # Taking last spectra at zero
    print('INFO: Taking last 0° spectra....')
    temp = flame.get_averaged_intensities()
    # Saving the data in the new scheme
    write_to_file(time()-start_time, current_angle, temp, path)
    # Saving the data in Mattias's scheme        
#    first_row[0,k + 5] = current_angle
#    data[:, k + 5] = temp
    
    # Plotting globals    
    TRACES.append(go.Scatter(x = WAVELENGTHS, y = temp, name = f'0°', mode = 'lines'))
    SRI.append(calculate_sri(WAVELENGTHS, temp - flame.background))
    CURRENT_ANGLE.append(current_angle)

    # Close shutter
    gonio.move_shutter()
    
    # Saving the data with Mattias' scheme
#    ftimestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
#    data = np.vstack([first_row, data])
#    path2 = pjoin(folder, timestamp + '_'+ filename + '_Mattias.dat')
#    header = itimestamp + '\n' + ftimestamp
#    np.savetxt(path2, data, fmt = '% 8.2f', header= header)
#    print('INFO: Measurement finished data saved at (Mattias scheme)\n\t' + path2)
    
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
    
    value_spec = 'No devices available' if LSPECTROMETERS == [] else  LSPECTROMETERS[0].serial_number
    value_arduino = 'No detected ports' if LPORTS == [] else LPORTS[0]

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
              Input('button-autozero', 'n_clicks')])
def gonio_and_spectra_functions(nleft, nright, nshutter, nbkg, nautozero):
    global gonio, traces
    # Determine which button has been clicked
    ctx = dash.callback_context

    if not ctx.triggered:
        button_id = 'No clicks yet'
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if button_id == 'button-move-left':
        gonio.move_angle(1.8 / 8)
    elif button_id == 'button-move-right':
        gonio.move_angle(-1.8/8)
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
            print(f'ERROR: No data to fit')
        
    else:
        pass
    return
        

if __name__ == '__main__':
    try:
        app.run_server(debug=True, port = 8051)
    except KeyboardInterrupt as e:
        print(e)
    except Exception as e:
        print(e)
        
    if flame is not None:
        flame.close()   
    if gonio is not None:
        gonio.close()