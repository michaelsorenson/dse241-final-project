import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import json
import plotly.graph_objs as go
import pandas as pd
import numpy as np
import plotly.express as px
import random
import os
import sys
from dotenv import load_dotenv
load_dotenv()
px.set_mapbox_access_token(os.getenv('MAPBOX_ACCESS_TOKEN'))

# custom functions
import plotly_utils as pu

# Define variables for figure sizes
# left column
MAP_WIDTH = 600
MAP_HEIGHT = 700
CIRCOS_WIDTH = 600
CIRCOS_HEIGHT = 600
# right columns
SUNBURST_WIDTH = 800
SUNBURST_HEIGHT = 800
STACKED_BAR_WIDTH = 800
STACKED_BAR_HEIGHT = 600

map_figs, sunburst_figs, stacked_bar_figs = {}, {}, {}

# Load metadata
if not os.path.exists(os.path.join('data', 'NCOG_sample_log_DNA_stvx_meta_2014-2020_mod.tsv')):
    raise Exception('Must have the metadata file under data/NCOG_sample_log_DNA_stvx_meta_2014-2020_mod.tsv\n' +
                    'Even if you are loading the data from JSON or HTML divs, metadata is still required')
data_meta = pd.read_csv(os.path.join('data', 'NCOG_sample_log_DNA_stvx_meta_2014-2020_mod.tsv'), sep='\t')

# If any of the JSON figure data doesn't exist, compute it from the raw data
fig_data_fnames = ['map_figure_data.json', 'stacked_bar_figure_data.json', 'sunburst_figure_data.json']
if not all([os.path.exists(os.path.join('views', 'data', fig_data_file)) for fig_data_file in fig_data_fnames]):
    # Required raw data files
    req_files = [os.path.join('data', fname) for fname in ['NCOG_18sV4_asv_count_tax.tsv', 'NCOG_18sV9_asv_count_tax_S.tsv'
                                                           'NCOG_21_16S_redo2_asv_count_tax.tsv']]
    # Raw data exists, computing data
    if all([os.path.exists(req_file) for req_file in req_files]):
        print('JSON figure data files not found. Input data files found in data directory, computing figure data...')
        # Load data
        data_16S = pd.read_csv(os.path.join('data', 'NCOG_21_16S_redo2_asv_count_tax.tsv'), sep='\t')
        data_18Sv4 = pd.read_csv(os.path.join('data', 'NCOG_18sV4_asv_count_tax.tsv'), sep='\t')
        data_18Sv9 = pd.read_csv(os.path.join('data', 'NCOG_18sV9_asv_count_tax_S.tsv'), sep='\t')

        # Map figures
        map_fig_data = pu.compute_map_data(data_meta)
        with open(os.path.join('views', 'data', 'map_figure_data.json'), 'w') as mapfile:
            mapfile.write(json.dumps(map_fig_data))

        # Sunburst figures
        sunburst_fig_data = pu.compute_sunburst_fig_data(data_meta, data_16S, data_18Sv4, data_18Sv9,
                                                         path=['Phylum', 'Class', 'Order', 'Family'])
        with open(os.path.join('views', 'data', 'sunburst_figure_data.json'), 'w') as sunburstfile:
            sunburstfile.write(json.dumps(sunburst_fig_data))
        
        # Stacked bar figures
        stacked_bar_fig_data, colors_16S, colors_18S = pu.compute_stacked_bar_data(data_meta, data_16S, data_18Sv4, data_18Sv9)
        with open(os.path.join('views', 'data', 'stacked_bar_figure_data.json'), 'w') as stackedbarfile:
            stacked_bar_data = {}
            stacked_bar_data['colors_16S'] = colors_16S
            stacked_bar_data['colors_18S'] = colors_18S
            stacked_bar_data['stacked_bar_fig_data'] = stacked_bar_fig_data
            stackedbarfile.write(json.dumps(stacked_bar_data))

    # Raw data doesn't exist, throwing error
    else:
        missing_files = ['\t' + req_file for req_file in req_files if not os.path.exists(req_file)]
        error_str = 'JSON figure data files not found and input data files missing. Missing the following input data files:'
        error_str += '\n'.join(missing_files)
        raise Exception(error_str)
# Load data - if JSON figure data already exists, load figure data, else compute figure data from raw data
else:
    print('JSON figure data files found, loading...')

    # Load figures using JSON data files
    # Map figures
    with open(os.path.join('views', 'data', 'map_figure_data.json'), 'r') as mapfile:
        map_fig_data = json.load(mapfile)

    # Sunburst figures
    with open(os.path.join('views', 'data', 'sunburst_figure_data.json'), 'r') as sunburstfile:
        sunburst_fig_data = json.load(sunburstfile)

    # Stacked bar figures
    with open(os.path.join('views', 'data', 'stacked_bar_figure_data.json'), 'r') as stackedbarfile:
        stacked_bar_data = json.load(stackedbarfile)
    colors_16S = stacked_bar_data['colors_16S']
    colors_18S = stacked_bar_data['colors_18S']
    stacked_bar_fig_data = stacked_bar_data['stacked_bar_fig_data']

map_figs = pu.precompute_map_figs(map_fig_data, width=MAP_WIDTH, height=MAP_HEIGHT)
sunburst_figs = pu.precompute_sunburst_figs(sunburst_fig_data, path=['Phylum', 'Class', 'Order', 'Family'],
                                            width=SUNBURST_WIDTH, height=SUNBURST_HEIGHT)
stacked_bar_figs = pu.precompute_stacked_bar_figs(stacked_bar_fig_data, colors_16S, colors_18S,
                                                  width=STACKED_BAR_WIDTH, height=STACKED_BAR_HEIGHT)
circos_figs = {
    '16S & 18Sv4': html.Img(src=os.path.join('assets', 'NCOG_network_v4_chordDiagram_positive.png'), style={
        'template': 'plotly_dark',
        'plot_bgcolor': 'rgba(0, 0, 0, 0)',
        'paper_bgcolor': 'rgba(0, 0, 0, 0)',
        'height': CIRCOS_HEIGHT,
        'width': CIRCOS_WIDTH
    }),
    '16S & 18Sv9': html.Img(src=os.path.join('assets', 'NCOG_network_v9_chordDiagram_positive.png'), style={
        'template': 'plotly_dark',
        'plot_bgcolor': 'rgba(0, 0, 0, 0)',
        'paper_bgcolor': 'rgba(0, 0, 0, 0)',
        'height': CIRCOS_HEIGHT,
        'width': CIRCOS_WIDTH
    })
}

# Lat and Lon center for plots
cal_coast_center = dict(
    lat=np.mean([min(data_meta['Lat_Dec']), max(data_meta['Lat_Dec'])]),
    lon=np.mean([min(data_meta['Lon_Dec']), max(data_meta['Lon_Dec'])])
)

# Environmental variables and sample types
env_var_cols = ['T_degC', 'Salnty', 'O2ml_L', 'PO4ug', 'SiO3ug', 'NO3ug', 'NH3ug', 'ChlorA', 'IntC14', 'NCDepth']
sample_type_vals = data_meta['sample_type'].dropna().unique()

# Add X before metadata sampleids to match 16S and 18S sample ids
data_meta['sampleid'] = data_meta['sampleid'].apply(lambda x: 'X' + x)

# Load empty map, sunburst, and stacked bar to show on page load
empty_map_fig, empty_sunburst_fig, empty_stacked_bar_fig = pu.load_empty_figs(cal_coast_center, MAP_WIDTH, MAP_HEIGHT, SUNBURST_WIDTH,
                                                                              SUNBURST_HEIGHT, STACKED_BAR_WIDTH, STACKED_BAR_HEIGHT)

def build_app(app):
    global map_figs, sunburst_figs, stacked_bar_figs, circos_figs, MAP_WIDTH, MAP_HEIGHT, SUNBURST_WIDTH, SUNBURST_HEIGHT, CIRCOS_WIDTH, CIRCOS_HEIGHT
    progress_val = 0

    app.layout = dbc.Container([
        html.H1("Marine microbes and their environment in NCOG stations off the California coast", className="mt-4 mb-4"),
        dbc.Row([
            dbc.Col([
                dbc.Card(
                    dbc.CardBody([
                        html.Div([
                            html.H3('Loading Sunburst and Stacked Bar Plots')
                        ]),
                        html.Div([
                            dbc.Progress(value=progress_val, striped=True, color="warning", id='progress-loading'),
                            dcc.Interval(id='interval-loading', interval=500, n_intervals=0, disabled=False)
                        ])
                    ], id='loading-card')
                ),
                dbc.Card(
                    dbc.CardBody([
                        html.Div([
                            html.Label('Sample Type:'),
                            dcc.Dropdown(sample_type_vals, 'Surf', id='sample-type-dropdown')
                        ])
                    ])
                )
            ], width=11.5)
        ]),
        dbc.Row([
            dbc.Col([
                dbc.Card(
                    dbc.CardBody([
                        html.H4('NCOG Stations vs Environmental Variables'),
                        html.Div([
                            dbc.Row([
                                dbc.Col([
                                    html.Label('Environmental Variable (Color):', style={'font-size': '18px'}),
                                ], width=5),
                                dbc.Col([
                                    dcc.Dropdown(['Total Number of Samples'] + env_var_cols, 'Total Number of Samples', id='env-var-dropdown')
                                ]),
                            ]),
                        ], className="env-var-dropdown"),
                        html.Div([
                            dcc.Graph(figure=empty_map_fig, id='map-graph')
                        ], className="map-graph"),
                    ]),
                ),
                dbc.Card(
                    dbc.CardBody([
                        html.Div([
                            html.H4('FlashWeave Region-wide Network Co-occurrences with Bacillariophyta'),
                            html.Div([
                                dbc.Row([
                                    dbc.Col([
                                        html.Label('Network Dataset:', style={'font-size': '18px'}),
                                    ], width=3),
                                    dbc.Col([
                                        dcc.Dropdown(['16S & 18Sv4', '16S & 18Sv9'], '16S & 18Sv4', id='circos-dropdown')
                                    ]),
                                ]),
                            ], className="circos-dropdown"),
                        ]),
                        html.Div(id='circos-container', style={'margin-top': '2vw'}),
                    ])
                ),
            ], width=6),
            dbc.Col([
                dbc.Card(
                    dbc.CardBody([
                        html.H4('Selected Station Taxonomic Makeup'),
                        html.Div([
                            dbc.Row([
                                dbc.Col([
                                    html.Label('Sunburst Dataset:', style={'font-size': '18px'}),
                                ], width=3),
                                dbc.Col([
                                    dcc.Dropdown(['16S', '18Sv4', '18Sv9'], '16S', id='dataset-dropdown')
                                ]),
                            ]),
                        ], className="dataset-dropdown"),
                        html.Div([
                            html.P(''),
                            html.Button('Reset View', id='reset-button')
                        ]),
                        html.Div([
                            dcc.Graph(figure=empty_sunburst_fig, id='sunburst-graph')
                        ]),
                    ])
                ),
                dbc.Card(
                    dbc.CardBody([
                        html.Div([
                            dcc.Graph(figure=empty_stacked_bar_fig, id='stacked-bar-graph')
                        ])
                    ])
                )
            ], width=6)
        ])
    ], fluid=True)


    # Data loading callback
    @app.callback(
        Output('progress-loading', 'value'),
        [Input('interval-loading', 'n_intervals'),
         Input('progress-loading', 'value')]
    )
    def update_interval(n, progress_value):
        if progress_value >= 100:
            return
        # this doesn't actually track the loading progress of the figure loading but it gives the illusion LOL
        if len(sunburst_figs) > 0 and len(stacked_bar_figs) > 0:
            progress_value = 100
        elif len(sunburst_figs) > 0 and progress_value < 95:
            if random.random() <= 0.2:
                progress_value += 1
        elif progress_value < 75:
            if random.random() <= 0.2:
                progress_value += 1
        return progress_value
    

    # Callback to hide loading bar and refresh data figures once all the data is done loading
    @app.callback(
        [Output(component_id='loading-card', component_property='style'),
        Output('progress-loading', 'color'),
        Output('interval-loading', 'disabled')],
        [Input('progress-loading', 'value')]
    )
    def hide_interval_when_finished(progress_value):
        if progress_value >= 100:
            return {'display': 'none'}, 'info', True
        return {'display': 'block'}, 'warning', False


    # Map graph dropdown callback
    @app.callback(
        Output('map-graph', 'figure'),
        [Input('map-graph', 'clickData'),
        Input('sample-type-dropdown', 'value'),
        Input('env-var-dropdown', 'value'),
        Input('interval-loading', 'disabled')]
    )
    def update_map(click_data, dropdown_sample_type, dropdown_env_var, is_done_loading):
        if not is_done_loading:
            return empty_map_fig
        if click_data is None:
            station_id = data_meta['Sta_ID'].iloc[0]
        # Get station ID from hover data
        else:
            if 'customdata' not in click_data['points'][0]:
                station_id = data_meta['Sta_ID'].iloc[0]
            else:
                station_id = click_data['points'][0]['customdata'][0]
        cur_fig = go.Figure(map_figs[dropdown_sample_type][dropdown_env_var])
        cur_fig.update_layout(
            template='plotly_dark',
            plot_bgcolor= 'rgba(0, 0, 0, 0)',
            paper_bgcolor= 'rgba(0, 0, 0, 0)',
            width=MAP_WIDTH,
            height=MAP_HEIGHT
        )
        stations = cur_fig.data[0]['customdata']
        stations = [x[0] for x in stations]
        if station_id in stations:
            clicked_station_idx = stations.index(station_id)
            cur_fig.data[0]['marker']['size'][clicked_station_idx] = 5
        return cur_fig
        

    # Map graph click data callback for sunburst
    @app.callback(
        Output('sunburst-graph', 'figure'),
        [Input('sample-type-dropdown', 'value'),
        Input('map-graph', 'clickData'),
        Input('dataset-dropdown', 'value'),
        Input('reset-button', 'n_clicks'),
        Input('interval-loading', 'disabled')]
    )
    def update_sunburst(dropdown_sample_type, click_data, dropdown_dataset, n_clicks, is_done_loading):
        if not is_done_loading:
            return empty_sunburst_fig
        if click_data is None:
            station_id = data_meta['Sta_ID'].iloc[0]
        # Get station ID from hover data
        else:
            if 'customdata' not in click_data['points'][0]:
                station_id = data_meta['Sta_ID'].iloc[0]
            else:
                station_id = click_data['points'][0]['customdata'][0]
        return sunburst_figs[station_id][dropdown_sample_type][dropdown_dataset].update_layout(
            template='plotly_dark',
            plot_bgcolor= 'rgba(0, 0, 0, 0)',
            paper_bgcolor= 'rgba(0, 0, 0, 0)',
            height=SUNBURST_HEIGHT,
            width=SUNBURST_WIDTH
        )


    # Stacked bar callback
    @app.callback(
        Output('stacked-bar-graph', 'figure'),
        [Input('map-graph', 'clickData'),
        Input('dataset-dropdown', 'value'),
        Input('interval-loading', 'disabled')]
    )
    def update_stacked_bar(click_data, dropdown_dataset, is_done_loading):
        if not is_done_loading:
            return empty_stacked_bar_fig
        if click_data is None:
            station_id = data_meta['Sta_ID'].iloc[0]
        # Get station ID from hover data
        else:
            if 'customdata' not in click_data['points'][0]:
                station_id = data_meta['Sta_ID'].iloc[0]
            else:
                station_id = click_data['points'][0]['customdata'][0]
        return stacked_bar_figs[station_id][dropdown_dataset].update_layout(
            template='plotly_dark',
            plot_bgcolor= 'rgba(0, 0, 0, 0)',
            paper_bgcolor= 'rgba(0, 0, 0, 0)',
            height=STACKED_BAR_HEIGHT,
            width=STACKED_BAR_WIDTH
        )


    # Circos png callback
    @app.callback(
        Output('circos-container', 'children'),
        [Input('circos-dropdown', 'value'),
         Input('interval-loading', 'disabled')]
    )
    def update_circos(circos_dropdown, is_done_loading):
        if is_done_loading:
            return circos_figs[circos_dropdown]

if __name__ == '__main__':
    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SLATE])
    if len(sys.argv) > 1:
        port = sys.argv[1]
    else:
        port = 8050
    build_app(app)
    app.run_server(host='0.0.0.0', port=port)