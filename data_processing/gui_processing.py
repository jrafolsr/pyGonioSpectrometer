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
import plotly.graph_objs as go
import numpy as np
import getopt
import sys
from pathlib import Path
import shutil
from pyGonioSpectrometer.data_processing import __file__ as calibration_dir
from pyGonioSpectrometer.data_processing  import process_goniodata, process_L0, load_sri_file,interpolate_expdata, gci,load_simdata, error_landscape
import seaborn as sns
## Globals

# Calibration files
calibration_dir = Path(calibration_dir).parent / 'calibration_files'
calibrations ={'arduino_gonio_1': dict(path_IRF = calibration_dir / 'IRF_FlameBlueFiber_wLens2945K.txt',\
                                       abs_calfactor = 7.10e6),\
               'raspberry_gonio2':dict(path_IRF = calibration_dir / 'IRF_Flame2OrangeFiber_wLens-F260SMA-A+100um-slit_PMA12-reference.txt',\
                                       abs_calfactor = 6.08E6),\
               'raspberry_gonio2_old':dict(path_IRF = calibration_dir / 'IRF_Flame2OrangeFiber_wLens2945K.txt',\
                                       abs_calfactor = 1.44E6)}
    
class ProcessingFolder(object):
    def __init__(self):
        self.processed_folder = 'processed-data'
        self.fits_folder = 'fits'
        self.plots_folder = 'plots'
        self.raw_folder = 'raw-data'
        self.folder = Path.home()
        self.files = []
        self.filesL0 = []
        self.iv_files = []
        self.iv_file = None
        self.files_processed = []
        
        self.processed_flag = False
        
        self.fitting_flag = False
        
        self.calibration = None
        self.angle_offset = 0.0
        self.current = 1.0
        self.thickness = None
        self.t0 = None
        self.correct_time_drift = True

        self.vtimes = []
        
        self.file_EZP_time = None
        self.file_time_luminance = None
        
    def check_folder(self, folder):
        self.__init__()
        
        self.folder = Path(folder)

        # Check if the folder contains proper files with *times-series* in the name
        files = list(self.folder.glob('*_time-series*.dat')) +\
            list((self.folder / self.raw_folder).glob('*_time-series*.dat'))
        if len(files) == 0:
            return False
        else:
            # I assume it is the right folder if it contains at least one file with *time-series* in the name
            # Create folders with output of processed data, plots and fits
            for subfolder in [self.raw_folder, self.processed_folder, self.fits_folder, self.plots_folder]:
                t = self.folder /  subfolder
                t.mkdir() if not t.exists() else None  
            
            files = list((self.folder / self.raw_folder).glob('*_time-series.dat'))
            filesL0 = list((self.folder / self.raw_folder).glob('*_time-series*.dat'))
            
            if len(files) == 0:
                files = folder.glob('*_time-series*.dat')
                for f in files:
                    shutil.move(f, self.folder/ self.raw_folder / f.name)
                files = list((self.folder / self.raw_folder).glob('*_time-series.dat'))
                filesL0 = list((self.folder / self.raw_folder).glob('*_time-series*.dat'))
        
            # Sort the files by name (and therfore time), do not assumed that glob returns them sorted (e.g. not in Linux at least)
            self.files = sorted(files)
            self.filesL0 = sorted(filesL0)          
            
            # Check if there is an existing processed files
            self.files_processed = sorted(list((self.folder / self.processed_folder).glob('*.sri')))
            self.files_processed_integrated = sorted(list((self.folder / self.processed_folder).glob('*.integrated')))
            if len(self.files_processed) > 0:
                self.processed_flag = True
            
            # Define existing or files to be created files in the case the folder contains processed data
            self.file_time_luminance = self.folder / 'time-luminance.dat'
            
            # Check if there exists a previous fitting
            self.file_EZP_time =  self.folder / self.fits_folder / 'time-evolution_emission-zone.dat'
            
            if self.file_EZP_time.exists():
                self.fitting_flag = True

            
            return True
    
    def process_data(self):
        if self.calibration == None or self.iv_file == None:
            return 'Define a calibration and/or iv_files!\n'
        else:

            # in deg, put 'auto' if you think it was misaligned
            current = self.current /1000 if self.current != None else None  # mA
            
            path_IRF = self.calibration['path_IRF']
            abs_calfactor = self.calibration['abs_calfactor']

            # Initial timestamp at the starting of the measurement
            self.t0 = np.loadtxt(self.iv_file, max_rows = 1, dtype = np.datetime64)
            
            # print(f'{len(files):d} gonio files')
            print(f'Processing {len(self.filesL0):d} files...')
    
            vtimes, vluminances = process_L0(self.filesL0, self.t0, folder = self.folder, path_IRF =  path_IRF, abs_calfactor = abs_calfactor)
            header = 'Luminance obtained from the forward spectrum\n'
            header += f'IRF_file: {path_IRF.name:s}\n'
            header += f'correct_time_drift: {self.correct_time_drift}\n'
            header += f'abs_calfactor: {abs_calfactor:4.2e}\n'
            header += f'angle_offset = {self.angle_offset}\n' if self.angle_offset == 'auto' else f'angle_offset = {self.angle_offset:.2f} deg\n'
            header +=  f't0 = {self.t0}\n'
            header += 'EllapsedTime(s)\tLuminance(cd/m2)'
            np.savetxt(self.folder / self.file_time_luminance, np.vstack([vtimes[1:], vluminances[1:]]).T, fmt = '%10.4f', header = header)
        
        # ----------------- Goniometer data processing --------------------------
            for i, file in enumerate(self.files):            
                _ = process_goniodata(file, current = current,\
                                            angle_offset= self.angle_offset,\
                                            plot = False,\
                                            correct_time_drift = self.correct_time_drift,\
                                            verbose =  False,\
                                            folder = str(self.folder / self.processed_folder),
                                            path_IRF = path_IRF,
                                            abs_calfactor  = abs_calfactor)
            
            self.file_time_luminance = self.folder / 'time-luminance.dat'
            self.files_processed = sorted(list((self.folder / self.processed_folder).glob('*.sri')))
            self.files_processed_integrated = sorted(list((self.folder / self.processed_folder).glob('*.integrated'))) 
            
            print('Done!')
            
            return 'Data processed.\n'
        
    def create_times_vector(self):
        """Creates a time vector based on the processed gobio files"""
        # Initial timestamp at the starting of the measurement
        self.t0 = np.loadtxt(self.iv_file, max_rows = 1, dtype = np.datetime64)
        files = self.files_processed_integrated
        
        vtimes = np.zeros((len(files), ))
        
        for i, file in enumerate(files):
            # Ellapsed time since start of the whole measurement
            t1 = np.loadtxt(file, max_rows = 1, dtype = np.datetime64)
            rTimes = np.loadtxt(file, skiprows = 12, usecols = 0)
            eTime = rTimes[(rTimes.shape[0] + 1) // 2]
            
            # Add the time of the goniometer at the middle of the measurement, in us
            t2 = t1 + int(eTime * 1e6) # The time needs to be add un us
    
            # Calculate the relative time between the starting of the measurement (t0) and the current measruement
            vtimes[i] = np.float64(t2 - self.t0) / 1e6
        self.vtimes = vtimes
    
    def fit_data_time_series(self, error_landscape_file, ethickness = None):       
        # Load the error landscape file
        simEL = load_simdata(error_landscape_file, wl_limit=(450,800), angle_max=70)
        
        #Create and additional folder for the error maps
        error_landscape_dir = self.folder / fits_folder / 'error_landscapes'
        error_landscape_dir.mkdir() if not  error_landscape_dir.exists() else None
        
        output = np.zeros((len(self.files_processed),7))
        # Thickness
        print(f'Fittings assuming a thickness of {self.thickness:.0f} nm')
        
        for i, file in enumerate(self.files_processed):
            # Calculates the EZP based on the ethickness range
            if ethickness != None:
                error_low, pos_min_low, _, _ = error_landscape(file, self.thickness - max(ethickness,5), simEL, folder = error_landscape_dir)
                error_high, pos_min_high, _, _ = error_landscape(file, self.thickness + max(ethickness,5), simEL, folder = error_landscape_dir)
                error_min_high = error_high.min()
                error_min_low = error_low.min()
            else:
               pos_min_high, pos_min_low, error_min_high, error_min_low = [np.nan]*4
               
            error, pos_min, _, _ = error_landscape(file, self.thickness, simEL, folder = error_landscape_dir)
            # output[i,:] = [vtimes[i], pos_min, error.min(), 0, 0, 0, 0]
            output[i,:] = [self.vtimes[i], pos_min, error.min(), pos_min_low, error_min_low, pos_min_high, error_min_high]
        
        if ethickness == None:
            ethickness = np.nan
            
        # Saving the data    
        header = str(self.t0) +'\n'
        header += f'{self.thickness:.0f} +/- {ethickness:.0f}\n'
        header += ('{:^10s}\t'*3).format('Rel.Time(s)','EZP(bestFit)', 'RMSE')
        np.savetxt(self.file_EZP_time, output, header = header, fmt = '%10.6g')    
        
        self.fitting_flag = True
        
    
    def compare_with_simulation(self, file, error_landscape_file, position, etime):
        # Load the error landscape file
        simEL = load_simdata(error_landscape_file, wl_limit=(450,800))
        wavelengths, angles, sri = load_sri_file(file)
    
        wl_sim = simEL['wl']
        angles_sim = simEL['angles']
        dAL_sim = simEL['dAL']
        ipos_sim = simEL['ipos']
        simData = simEL['data']
        
        iNormExpSRI = interpolate_expdata(sri, wavelengths, angles, wavelengths, angles_sim)
        
        # Find the simulation data according to the input thickness   
        ithickness = gci(self.thickness, dAL_sim)
            
        # So, we choose the column according to the experimental thickness:
        thickness_sim = dAL_sim[ithickness] # simulated thickness
        # Find the indices of positions in the ipos_sim vector
        iposition = gci(position, ipos_sim)
        
        NormSimSRI = simData[iposition,ithickness] # Columns are the simulated thicknesses
        
        data2save = np.zeros((len(wavelengths) + 1, len(angles_sim) + 1)) * np.nan
        data2save[0, 1:] = angles_sim
        data2save[1:, 0] = wavelengths
        data2save[1:, 1:] = iNormExpSRI
        
        header = f'thickness_exp = {p.thickness:.1f} nm\nthickness_sim = {thickness_sim:.0f} nm\nEZP = {position:.2f}'
        
        sufix = f'_t={etime:.0f}s_dsim={thickness_sim:.0f}nm_EZP={position:.2f}_'
        
        np.savetxt(self.folder / self.fits_folder / (file.stem + sufix + 'measured-spectra.dat'), data2save, header=header, fmt  = '%.4f')
        
        data2save = np.zeros((len(wl_sim) + 1, len(angles_sim) + 1)) * np.nan
        data2save[0, 1:] = angles_sim
        data2save[1:, 0] = wl_sim
        data2save[1:, 1:] = NormSimSRI
                
        np.savetxt(self.folder / self.fits_folder / (file.stem + sufix + 'simulated-spectra.dat'), data2save, header=header, fmt  = '%.4f')
        
        return angles_sim, wavelengths, iNormExpSRI, wl_sim, NormSimSRI, thickness_sim

p = ProcessingFolder()

# Output folders
processed_folder = 'processed-data'
fits_folder = 'fits'
plots_folder = 'plots'
raw_folder = 'raw-data'

#%%    

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

plot_data = [go.Scatter(x=[], y=[], mode = 'lines+markers'),
          go.Scatter(x=[], y=[], name = '', mode = 'lines',\
             xaxis = 'x2', yaxis = 'y2'),
          go.Scatterpolar(r=[], theta=[], mode = 'lines'),
          go.Scatter(x = [], y = [], mode = 'lines+markers',\
             xaxis = 'x', yaxis = 'y3'),\
          go.Scatter(x = [], y = [], mode = 'lines+markers',\
             xaxis = 'x4', yaxis = 'y4')]

plot_layout = dict(height=800,#, width=1200,
                    margin =  {'l': 60, 'r': 60, 'b': 60, 't': 20},\
                   legend = dict(yanchor="top",
                                y=0.99,
                                xanchor="right",
                                x=0.99),\
                   xaxis = dict(title =  "Time (s)",
                            anchor = 'y',
                            range = [0, None],
                            domain=[0, 1]),\
                   yaxis = dict(range  =  [0, None],
                            anchor = 'x',
                            title =  "Voltage (V)",
                            domain = [0.65, 1]),
                   xaxis2 = dict(title =  "Wavelengths (nm)",\
                            anchor = 'y2',
                            range = [400, 800],
                            domain=[0, 0.45]),
                   yaxis2 = dict(title =  "Spectral radiant Intensity (a.u.)",\
                            domain = [0.20, 0.55] ,
                            # range = [0, None],
                            anchor="x2"),
                   yaxis3 = dict(anchor =  'x',
                               overlaying = 'y',
                               side = 'right',
                               title = {'text': 'Luminance (cd/m2'}),  
                   polar =  dict(title = 'Angular profile',\
                                 domain = {'x': [0.55, 1.0],
                                           'y': [0.20, 0.55]},
                                 sector  =  (0,180),
                                 angularaxis = dict(rotation= 90)),
                   xaxis4 = dict(title =  "Emission zone position",\
                                 anchor = 'y4',
                                range = [0, 1],
                                domain=[0, 1]),
                   yaxis4 = dict(title =  "RMSE",\
                            domain = [0, 0.10] ,
                            # range = [0, None],
                            anchor="x4"),
                   annotations = [{'font': {'size': 18},
                                    'showarrow': False,
                                    'text': 'Transient',
                                    'x': 0.5,
                                    'xanchor': 'center',
                                    'xref': 'paper',
                                    'y': 1.0,
                                    'yanchor': 'bottom',
                                    'yref': 'paper'},
                                 {'font': {'size': 16},
                                    'showarrow': False,
                                    'text': '',
                                    'x': 0.225,
                                    'xanchor': 'center',
                                    'xref': 'paper',
                                    'y': 0.56,
                                    'yanchor': 'bottom',
                                    'yref': 'paper'} ]
                   # xaxis3 = dict(title =  "Wavelengths (nm)",\
                   #          range = [-90, 90],
                   #          domain=[0.55, 1],
                   #          type = 'scatterpolar'),
                   # yaxis3 = dict(title =  "Spectral radiant Intensity (a.u.)",\
                   #          domain = [0, 0.375] ,\
                   #          # range = [0, None],
                   #          anchor="x3",
                   #          type = 'scatterpolar')
                   )

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.title = 'Data plotter'
app.layout = html.Div(children =  [
        html.Div(id='live-update-text', className = 'row', children  = [
        html.Div(className = 'column left', children = [       
            html.H4('Goniospectrometer data plotter', style = {'width': '100%', 'display': 'inline-block','vertical-align':'middle'}),
            html.Span([], id = 'trigger-update-graph', hidden  = True),
            dcc.Graph(id='live-update-graph', 
                figure= { "data": plot_data,
                          "layout": plot_layout
                          }
                ),
            html.Span('Time slider'),
            dcc.Slider(
                id = 'slider-time',
                disabled = True,
                min = 0,
                max = 0,
                value = 0,
                step = 1,
                updatemode='drag'
            ),
            ],    
        ),
        
        html.Div(className =  'column middle', children = [
         daq.StopButton(id='button-active-plot',
               disabled = False,
               buttonText = 'plot',
               n_clicks = 0,
               ),
         daq.StopButton(id='button-compare',
               disabled = True,
               buttonText = 'compare',
               n_clicks = 0,
               ),
         daq.StopButton(id='button-reset',
               disabled = False,
               buttonText = 'reset',
               n_clicks = 0,
               ),
         html.Span(id = 'motor-movement', children = [], hidden = True)
        ]),  
        html.Div(id  = 'div-inputs', className = 'column right', children = [
        html.Div(id = 'input-process-container', className = 'container', children = [
                html.Div([
                    html.Div([
                        daq.StopButton(id='button-process-folder',
                           disabled = False,
                           buttonText = 'process',
                           n_clicks = 0),
                        
                        'Folder: ',          
                        dcc.Input(id="folder-input",
                                  type="text",
                                  placeholder="Folder",
                                  value = str(Path.home()),
                                  size = '40',
                                  debounce = True)
                        ], style = {'width' : '60%', 'display' : 'inline-block','vertical-align' : 'top'}),
                    html.Div([
                          'Input current (mA):',              
                          dcc.Input(
                                id = "current-input",
                                type = 'number',
                                value = 1,
                                size = '5',
                                debounce = True),
                          html.Br(),
                          ' Angle offset (°):',              
                          dcc.Input(
                                id = "angle-offset-input",
                                type = 'number',
                                value = 0.0,
                                size = '5',
                                debounce = True),
                          dcc.RadioItems(
                                  id = 'radio-angle-offset',
                                  options=[
                                      {'label': 'Fix offset', 'value': 0},
                                      {'label': 'Auto offset', 'value': 'auto'}
                                  ],
                                  value = 0
                              )  
                      ], style = {'width' : '40%', 'display' : 'inline-block','vertical-align' : 'middle'}),
                    ]),

              html.Div( [
                  html.Div( children = [
                      'Available raw files:',
                        dcc.Textarea(
                          id='textarea-files',
                          value='The available files to process show be listed here.',
                          style={'width': '100%', 'height': 100},
                          readOnly = True
                          )
                          ], style = {'display': 'inline-block', 'width': '45%', 'padding-right' : '5px'}),
                html.Div( children = [
                        '(Already) available processed files:',
                        dcc.Textarea(
                          id='textarea-processed-files',
                          value='If there are existing processed files they will be showed here.',
                          style={'width': '100%', 'height': 100},
                          readOnly = True
                          )
                    ], style = {'display': 'inline-block', 'width': '45%', 'padding-left' : '5px'}
                    )
                  ], style = {'padding-top' : '10px'}),
                  'Pick IV file:',
                  dcc.Dropdown(id  = 'dropdown-IV-files',
                    options = [],
                    value = None,
                    placeholder = 'No IV files available',
                    style = {'width' : '200'},
                    searchable = False
                ),
                  'Pick spectrometer calibration file:',
                  dcc.Dropdown(id  = 'dropdown-calibration',
                    options = [{'label' : key, 'value' : key} for key in calibrations],
                    value = None,
                    placeholder = 'Pick the calibration file',
                    style = {'width' : '200'},
                    searchable = False
                ),
                    
                  ]),
            html.Div( id = 'fits-input-container', className = 'container', children = [
                daq.StopButton(id='button-fit',
                       disabled = False,
                       buttonText = 'fit',
                       n_clicks = 0,
                       style = {'display' : 'inline-block'}),
                      'Thickness (nm): ',
                      dcc.Input(
                          id = "input-thickness",
                          type = 'number',
                          value = None,
                          size = '5',
                          debounce = True),
                      html.P('Simulation file:'),
                      dcc.Dropdown(id  = 'dropdown-error-landscape',
                        options = [{'label' : 'ErrorLandscape_newPL', 'value' : str(calibration_dir /'ErrorLandscape_newPL.mat')}],
                        value = str(calibration_dir /'ErrorLandscape_newPL.mat'),
                        placeholder = 'Pick the error landscape file',
                        style = {'width' : '200'},
                        searchable = False
                        ),
              ]),
            html.Div( id = 'debugger-container', className = 'container', children = [
                  html.H6('Debugger'),
                  dcc.Textarea(
                    id='textarea-logger',
                    value='',
                    style={'width': '100%', 'height': 100},
                    readOnly = True
                    ),
              ])
            ])
       ])
    ])


@app.callback(Output('live-update-graph', 'figure'),
              [Input('button-active-plot', 'n_clicks'),
               Input('button-compare', 'n_clicks'),
               Input('button-reset', 'n_clicks'),
               Input('folder-input', 'value'),
               Input('slider-time', 'value')],
              [State('live-update-graph', 'figure'),
               State('dropdown-error-landscape', 'value')],
              prevent_initial_call =  True)

def update_plot(n_clicks, n_clicks2,n_clicks3, folder_button, ifile, figure, error_landscape_file):


    # Determine which button has been clicked
    ctx = dash.callback_context

    if not ctx.triggered:
        button_id = 'No clicks yet'
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if button_id in['button-reset', 'folder-input'] :
        figure['data'] = plot_data
        figure['layout']['annotations'][0]['text'] = 'Transient'
        figure['layout']['annotations'][1]['text'] = ''
        figure['layout']['shapes'] = []
        return figure
    
    if p.iv_file == None or p.files_processed == []:
        return figure
    
    
    file = p.files_processed[ifile] 
    
    # Plot voltage
    try:
        timeV, voltageV = np.loadtxt(p.iv_file, usecols = (1,3), unpack = True, skiprows = 3)
    except Exception as e:
        print('There is some error with the iv file, check it.', e)
        timeV, voltageV = [], []
        
    #Plot luminance
    timeL, luminanceL = np.loadtxt(p.file_time_luminance, usecols = (0,1), unpack = True)


    
    # Reset figure data
    figure['data'] = []
    
    figure['data'].append(go.Scatter(x = timeV, y = voltageV, name = 'voltage', mode = 'lines', xaxis = 'x', yaxis = 'y'))
    
    figure['data'].append(go.Scatter(x = timeL, y = luminanceL, name = 'luminance', mode = 'lines+markers',\
             xaxis = 'x', yaxis = 'y3'))
    figure['data'].append(go.Scatter(x = [], y = [], name = 'error', mode = 'lines+markers',\
             xaxis = 'x4', yaxis = 'y4') )  
    
        # Polar plot
    polar_angles, polar_sri = np.loadtxt(p.files_processed_integrated[ifile], usecols = (1,3), skiprows = 12, unpack=True)
    figure['data'].append(go.Scatterpolar(theta=polar_angles, r=polar_sri, mode = 'lines + markers', showlegend=False))
    lambertian_theta = np.linspace(-90,90, 100)
    lambertian_emitter = np.cos(lambertian_theta * np.pi /180)
    figure['data'].append(go.Scatterpolar(theta=lambertian_theta, r=lambertian_emitter, mode = 'lines', line_color = 'black', line_dash ='dash',showlegend=False))
    
    # Adding text indicating hte time you are plotting
    figure['layout']['annotations'][0]['text'] = f'Transient, time = {p.vtimes[ifile]:.0f} s'
    
    # Adding the vertical line
    figure['layout']['shapes'] =  [{'line': {'color': 'black', 'dash': 'dash', 'width': 1},
                        'type': 'line',
                        'x0': p.vtimes[ifile],
                        'x1': p.vtimes[ifile],
                        'xref': 'x',
                        'y0': 0.65,
                        'y1': 1,
                        'yref': 'paper'}]
    
    figure['layout']['annotations'][1]['text'] ='Click compare to plot the fitting'
    
    if button_id == 'button-active-plot' or button_id == 'slider-time' or p.thickness == None:
        wavelengths, angles, sri = load_sri_file(file)
        
        N = int((len(angles) + 1) / 2)
        plot_angles = angles[N:]
        plot_sri = sri[:, N:]
        plot_sri /= plot_sri[:,0].max()
        # Plot spectra
        c =  sns.cubehelix_palette(N, start=.5, rot=-.75, reverse = True)
        
        N = len(plot_angles)
        
        for i, angle in enumerate(plot_angles):
            offset = 0.125 * (N - 1 - i)
            color = 'rgba({:.4f},{:.4f},{:.4f},1)'.format(*c[i])
            figure['data'].append(go.Scatter(x = wavelengths, y = offset + plot_sri[:,i], name = f'{plot_angles[i]:.1f}°', mode = 'lines',\
                 xaxis = 'x2', yaxis = 'y2', line_color = color, showlegend=False))
        
    elif button_id == 'button-compare' and p.fitting_flag:
        if p.thickness == None:
            return figure      
        
        _, ezp_positions = np.loadtxt(p.file_EZP_time, usecols=(0,1), unpack=True)
        
        best_position = ezp_positions[ifile]

        angles, wavelengths, iNormExpSRI, wl_sim, NormSimSRI, thickness_sim = p.compare_with_simulation(file, error_landscape_file, best_position, p.vtimes[ifile])
        N = len(angles)
        
        c =  sns.cubehelix_palette(len(angles), start=.5, rot=-.75, reverse = True)
        
        for i, angle in enumerate(angles):
            offset = 0.125 * (N - 1 - i)
            color = 'rgba({:.4f},{:.4f},{:.4f},1)'.format(*c[i])
            figure['data'].append(go.Scatter(x = wavelengths, y = offset + iNormExpSRI[:,i], name = f'{angles[i]:.1f}°', mode = 'lines',\
                 xaxis = 'x2', yaxis = 'y2', line_color = color, showlegend=False))
            figure['data'].append(go.Scatter(x = wl_sim, y = offset + NormSimSRI[:,i], mode = 'lines',\
                 xaxis = 'x2', yaxis = 'y2', line_color = color,line_dash = 'dash', showlegend=False))
            
        epositions, error = np.loadtxt(p.folder / p.fits_folder / Path('error_landscapes') / (file.stem + '.el'), unpack=True, usecols = (0,1))
        
        # Add error plot
        figure['data'].append(go.Scatter(x = epositions, y = error, name = 'error', mode = 'lines+markers',\
             xaxis = 'x4', yaxis = 'y4', showlegend=False))  
            
        figure['layout']['annotations'][1]['text'] = f'Best EZP fit = {best_position:.2f}, dsim = {thickness_sim:.0f} nm'
        
        # Adding the vertical line
        figure['layout']['shapes'].append({'line': {'color': 'black', 'dash': 'dash', 'width': 1},
                            'type': 'line',
                            'x0': best_position,
                            'x1': best_position,
                            'xref': 'x4',
                            'y0': 0,
                            'y1': .1,
                            'yref': 'paper'})
        
    return figure

@app.callback([Output('folder-input', 'value'),
               Output('input-thickness', 'value'),
               Output('dropdown-calibration', 'value')],
              [Input('button-reset', 'n_clicks')])
def reset_all(n_clicks):
    folder = str(Path.home())
    p.__init__()
    thickness = None
    calibration = None
    return folder, thickness, calibration

@app.callback([Output('textarea-logger', 'value'),
               Output('button-compare', 'disabled'),
               Output('slider-time', 'disabled'),
               Output('slider-time', 'max'),
               ],
              [Input('current-input', 'value'),
               Input('input-thickness', 'value'),
               Input('dropdown-calibration', 'value'),
               Input('dropdown-IV-files', 'value'),
               Input('angle-offset-input', 'value'),
               Input('button-process-folder', 'n_clicks'),
               Input('button-active-plot', 'n_clicks'),
               Input('button-fit', 'n_clicks'),\
               Input('button-compare', 'n_clicks'),
               Input('button-reset', 'n_clicks')
               ],
              [State('textarea-logger', 'value'),
               State('dropdown-error-landscape', 'value'),
               State('button-compare', 'disabled'),
               State('slider-time', 'disabled'),
               State('slider-time', 'max'),
               ])
def update_processing_parameters(current, thickness, calibration, iv_file, angle_offset,\
                                 n_clicks, n_clicks2, n_clicks3,n_clicks4, n_click5,\
                                 text, error_landscape_file, disable_comparing, disable_slider, max_slider):
        # Determine which button has been clicked
    ctx = dash.callback_context

    if not ctx.triggered:
        button_id = 'No clicks yet'
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if button_id == 'current-input':
        p.current = current
        text += f'Input current of {current:.2f} mA.\n'
    
    elif button_id == 'dropdown-calibration':
        p.calibration = calibrations[calibration] if calibration != None else None 
        text += f'Calibration "{calibration}" selected.\n'
    
    elif button_id == 'dropdown-IV-files':
        if not iv_file == None:
            p.iv_file = p.folder  / iv_file
            if p.processed_flag:
                p.create_times_vector()
            text += f'IV-file set to "{iv_file}".\n'
    
    elif button_id == 'angle-offset-input':
        p.angle_offset = angle_offset
        text += f'Angle offset set to {angle_offset}.\n'
    
    elif button_id == 'button-process-folder':
        if p.files == []:
            text += 'The folder does not contain any files to process'
        else:
            text += p.process_data()
            p.create_times_vector()
    
    elif button_id == 'button-active-plot':
        if p.iv_file ==  None:
            text += 'Please select an IV-file from the dropdown menu.\n'
            
        N = len(p.files_processed)
        if N == 0:
            text += 'Please, process the files first.\n'
        
    
        max_slider = (N - 1) if N > 0 else 0
    
        disable_slider = True if  max_slider == 0 else False
        disable_comparing = False
        
    elif button_id == 'input-thickness':
        p.thickness = thickness
        if thickness != None:
            text += f'Input thickness of {thickness} nm.\n'
    
    elif button_id == 'button-fit':
        if p.thickness == None or p.t0 == None:
            text += 'Please define a thickness and/or an IV-file to proceed with the fitting.\n'
        else:
            text += 'Data fitted.'
            print(error_landscape_file)
            p.fit_data_time_series(error_landscape_file)
    
    elif button_id == 'button-compare':
        if not p.fitting_flag:
            text += 'Fitting is missing for the comparison.\n'
        elif p.thickness == None:
            text += 'Thickness missing for the comparison.\n'
        else:
            text += f'Comparing simulated and experimental data for {p.thickness:.0f} nm.\n'
    
    elif button_id == 'button-reset':
        text = ''
        disable_comparing = True
        disable_slider = True
        max_slider = 0
        
    return text, disable_comparing, disable_slider, max_slider


@app.callback([Output('textarea-files', 'value'),
              Output('dropdown-IV-files', 'options'),
              Output('dropdown-IV-files', 'value'),
              Output('textarea-processed-files', 'value')],
              [Input('folder-input', 'value')],)
def load_files(folder):
    folder = Path(folder)
    
    if not folder.exists():
        text = 'ERROR: The input folder does not exist'
        return text, [], None, text
    else:
        if p.check_folder(folder):
            if len(p.filesL0) == 0:
                text = 'No available files.'
            else:
                text = ''
                for f in p.filesL0:
                    text += f.name + '\n'
        
            options_iv_files = [{'label' : f.name, 'value' : f.name} for f in folder.glob('*.txt')]
            value_iv_file = options_iv_files[0]['value'] if len(options_iv_files) else None 

            if p.processed_flag:
                text_processed = ''
                for f in p.files_processed:
                    text_processed += f.name + '\n'
            else:
                text_processed = 'If there are existing processed files they will be showed here.'
        else:
            text = text_processed = 'Check if you input the correct folder\n(It does not seem to contain any *time-series* file.)'
            options_iv_files = []
            value_iv_file = None
        return text, options_iv_files, value_iv_file, text_processed
    
@app.callback([Output('angle-offset-input', 'disabled'),
               Output('angle-offset-input', 'value')],
              [Input('radio-angle-offset', 'value')])
def set_angle_offset_mode(value):
    disabled = False
    
    if value =='auto':
        disabled = True
    else:
        value = 0.0
    
    return disabled, value
    

if __name__ == '__main__':
   # Default values
    debug = True
    port = 8054
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
