import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
import pandas as pd
import numpy as np
import plotly.express as px
import os
import sys
from dotenv import load_dotenv
load_dotenv()
px.set_mapbox_access_token(os.getenv('MAPBOX_ACCESS_TOKEN'))

# custom functions
import plotly_utils as pu

# Load data
data_16S = pd.read_csv('NCOG_21_16S_redo2_asv_count_tax.tsv', sep='\t')
data_18Sv4 = pd.read_csv('NCOG_18sV4_asv_count_tax.tsv', sep='\t')
data_18Sv9 = pd.read_csv('NCOG_18sV9_asv_count_tax_S.tsv', sep='\t')
data_meta = pd.read_csv('NCOG_sample_log_DNA_stvx_meta_2014-2020_mod.tsv', sep='\t')

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

# precompute map figures for different sample types and environmental variables
map_figs = pu.precompute_map_figs(data_meta)

# precompute sunburst figures
sunburst_figs = pu.precompute_sunburst_figs(data_meta, data_16S, data_18Sv4, data_18Sv9)

# precompute stacked bar figures
stacked_bar_figs = pu.precompute_stacked_bar_figs(data_meta, data_16S, data_18Sv4, data_18Sv9)

def build_app(app):
    app.layout = dbc.Container([
        html.H1("Marine microbial life at California coastal stations", className="mt-4 mb-4"),
        dbc.Row([
            dbc.Col([
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
                        html.Div([
                            html.Label('Environmental Variable (Color):'),
                            dcc.Dropdown(['Total Number of Samples'] + env_var_cols, 'Total Number of Samples', id='env-var-dropdown')
                        ], className="env-var-dropdown")
                    ])
                ),
                dbc.Card(
                    dbc.CardBody([
                        html.Div([
                            dcc.Graph(id='map-graph')
                        ], className="map-graph"),
                    ]),
                ),
                dbc.Card(
                    dbc.CardBody([
                        html.Div([
                            html.Label('Circos Dataset:'),
                            dcc.Dropdown(['16S & 18Sv4', '16S & 18Sv9'], '16S & 18Sv4', id='circos-dropdown')
                        ]),
                        html.Div(id='circos-container', style={'margin-top': '2vw'}),
                    ])
                ),
            ], width=6),
            dbc.Col([
                dbc.Card(
                    dbc.CardBody([
                        html.Div([
                            html.Label('Dataset:'),
                            dcc.Dropdown(['16S', '18Sv4', '18Sv9'], '16S', id='dataset-dropdown')
                        ]),
                    ])
                ),
                dbc.Card(
                    dbc.CardBody([
                        html.Div([
                            html.Button('Reset View', id='reset-button')
                        ]),
                        html.Div([
                            dcc.Graph(id='sunburst-graph')
                        ]),
                    ])
                ),
                dbc.Card(
                    dbc.CardBody([
                        html.Div([
                            dcc.Graph(id='stacked-bar-graph')
                        ])
                    ])
                )
            ], width=6)
        ])
    ], fluid=True)


    # Map graph dropdown callback
    @app.callback(
        Output('map-graph', 'figure'),
        [Input('sample-type-dropdown', 'value'),
        Input('env-var-dropdown', 'value')]
    )
    def update_map(dropdown_sample_type, dropdown_env_var):
        return map_figs[dropdown_sample_type][dropdown_env_var].update_layout(
            template='plotly_dark',
            plot_bgcolor= 'rgba(0, 0, 0, 0)',
            paper_bgcolor= 'rgba(0, 0, 0, 0)',
            width=600,
            height=600
        )

    # Map graph click data callback for sunburst
    @app.callback(
        Output('sunburst-graph', 'figure'),
        [Input('sample-type-dropdown', 'value'),
        Input('map-graph', 'clickData'),
        Input('dataset-dropdown', 'value'),
        Input('reset-button', 'n_clicks')]
    )
    def update_sunburst(dropdown_sample_type, click_data, dropdown_dataset, n_clicks):
        if click_data is None:
            station_id = data_meta['Sta_ID'].iloc[0]
        # Get station ID from hover data
        else:
            if 'customdata' not in click_data['points'][0]:
                station_id = station_id = data_meta['Sta_ID'].iloc[0]
            else:
                station_id = click_data['points'][0]['customdata'][0]
        return sunburst_figs[station_id][dropdown_sample_type][dropdown_dataset].update_layout(
            template='plotly_dark',
            plot_bgcolor= 'rgba(0, 0, 0, 0)',
            paper_bgcolor= 'rgba(0, 0, 0, 0)',
            height=600,
            width=600
        )

    # Stacked bar callback
    @app.callback(
        Output('stacked-bar-graph', 'figure'),
        [Input('map-graph', 'clickData'),
        Input('dataset-dropdown', 'value')]
    )
    def update_stacked_bar(click_data, dropdown_dataset):
        if click_data is None:
            station_id = data_meta['Sta_ID'].iloc[0]
        # Get station ID from hover data
        else:
            if 'customdata' not in click_data['points'][0]:
                station_id = station_id = data_meta['Sta_ID'].iloc[0]
            else:
                station_id = click_data['points'][0]['customdata'][0]
        return stacked_bar_figs[station_id][dropdown_dataset].update_layout(
            template='plotly_dark',
            plot_bgcolor= 'rgba(0, 0, 0, 0)',
            paper_bgcolor= 'rgba(0, 0, 0, 0)',
            height=600,
            width=600
        )
    @app.callback(
        Output('circos-container', 'children'),
        [Input('circos-dropdown', 'value')]
    )
    def update_circos(circos_dropdown):
        if circos_dropdown == '16S & 18Sv4':
            src_path = 'assets/NCOG_network_v4_chordDiagram_positive.png'
        elif circos_dropdown == '16S & 18Sv9':
            src_path = 'assets/NCOG_network_v9_chordDiagram_positive.png'
        return html.Img(src=src_path, style={
            'template':'plotly_dark',
            'plot_bgcolor':'rgba(0, 0, 0, 0)',
            'paper_bgcolor':'rgba(0, 0, 0, 0)',
            'height':600,
            'width':600
        })

if __name__ == '__main__':
    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SLATE])
    if len(sys.argv) > 1:
        port = sys.argv[1]
    else:
        port = 8050
    build_app(app)
    app.run_server(host='0.0.0.0', port=port)
