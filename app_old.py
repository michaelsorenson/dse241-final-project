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

# Load data
data_16S = pd.read_csv('NCOG_21_16S_redo2_asv_count_tax.tsv', sep='\t')
data_18Sv4 = pd.read_csv('NCOG_18sV4_asv_count_tax.tsv', sep='\t')
data_18Sv9 = pd.read_csv('NCOG_18sV9_asv_count_tax_S.tsv', sep='\t')
data_meta = pd.read_csv('NCOG_sample_log_DNA_stvx_meta_2014-2020_mod.tsv', sep='\t')

# Add X before metadata sampleids to match 16S and 18S sample ids
data_meta['sampleid'] = data_meta['sampleid'].apply(lambda x: 'X' + x)

# Lat and Lon center for plots
cal_coast_center = dict(
    lat=np.mean([min(data_meta['Lat_Dec']), max(data_meta['Lat_Dec'])]),
    lon=np.mean([min(data_meta['Lon_Dec']), max(data_meta['Lon_Dec'])])
)

# Environmental variables and sample types
env_var_cols = ['T_degC', 'Salnty', 'O2ml_L', 'PO4ug', 'SiO3ug', 'NO3ug', 'NH3ug', 'ChlorA', 'IntC14', 'NCDepth']
sample_type_vals = data_meta['sample_type'].dropna().unique()

# precompute map figures for different sample types and environmental variables
map_figs = {sample_type: {} for sample_type in data_meta['sample_type'].dropna().unique()}
for sample_type in data_meta['sample_type'].dropna().unique():
    for env_var in env_var_cols:
        meta_subset = data_meta[data_meta['sample_type'] == sample_type]
        meta_subset = meta_subset[['Sta_ID', 'Lat_Dec', 'Lon_Dec', 'sample_type', env_var]].groupby('Sta_ID').agg({
            'Lat_Dec': 'min',
            'Lon_Dec': 'min',
            'sample_type': 'count',
            env_var: 'mean'
        }).rename({'sample_type': 'Number of Samples'}, axis=1).reset_index()
        hover_names = meta_subset['Sta_ID'].apply(lambda x: '<b>Station: </b>' + x)
        subset_fig = px.scatter_mapbox(meta_subset, lat='Lat_Dec', lon='Lon_Dec', center=cal_coast_center,
                                       color=env_var, hover_name=hover_names, hover_data='Number of Samples', #size="num_samples",
                                       color_continuous_scale='viridis', size_max=15, zoom=4.5, mapbox_style='outdoors',
                                       width=600, height=700, custom_data='Sta_ID')
        map_figs[sample_type][env_var] = subset_fig

# empty sunburst figure for if the station has no data with the selected sample type
empty_sunburst_data = {'Phylum': [], 'Class': [], 'Order': []}
for parent in ['Undetermined_1', 'Undetermined_2', 'Undetermined_3']:
    for child1 in ['Undetermined_1', 'Undetermined_2', 'Undetermined_3']:
        for child2 in ['Undetermined_1', 'Undetermined_2', 'Undetermined_3']:
            empty_sunburst_data['Phylum'].append(parent)
            empty_sunburst_data['Class'].append(child1)
            empty_sunburst_data['Order'].append(child2)
empty_sunburst_fig = px.sunburst(empty_sunburst_data, path=['Phylum', 'Class', 'Order'])

# precompute sunburst figures for every station id, sample type, dataset combination
counter = 0
total_computing = len(data_meta['Sta_ID'].dropna().unique()) * len(data_meta['sample_type'].dropna().unique()) * 3
print(f'Pre-computing {total_computing} Sunburst Figures')
sunburst_figs = {station_id: {} for station_id in data_meta['Sta_ID'].dropna().unique()}
taxonomies_all = {station_id: {} for station_id in data_meta['Sta_ID'].dropna().unique()}
for station_id in data_meta['Sta_ID'].dropna().unique():
    for sample_type in data_meta['sample_type'].dropna().unique():
        sunburst_figs[station_id][sample_type] = {}
        taxonomies_all[station_id][sample_type] = {}
        for dataset in ['16S', '18Sv4', '18Sv9']:
            cols_show_in_sunburst = ['Phylum', 'Class', 'Order']
            station_data = data_meta[(data_meta['Sta_ID'] == station_id) & (data_meta['sample_type'] == sample_type)]
            station_samples = station_data['sampleid'].tolist()
            if dataset =='16S':
                taxa_col_names = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

                # Merge with 16S dataframe to get taxonomy data for the samples
                asv_cols = pd.Series(data_16S.columns).isin(station_samples).values
                if np.sum(asv_cols) == 0: # If there are no samples for the sample_type at the station
                    # set the figure to an empty figure
                    title = f'No {sample_type} samples for station "{station_id}"'
                    fig = px.sunburst(empty_sunburst_data, path=['Phylum', 'Class', 'Order'], title=title)
                    sunburst_figs[station_id][sample_type][dataset] = fig
                    continue
                asv_cols[0] = True
                station_asvs = pd.concat([data_16S.loc[:,asv_cols], data_16S['silva_Taxon']], axis=1)

                # Get relative abundances
                values = station_asvs.drop(['Feature.ID', 'silva_Taxon'], axis=1).fillna(0).sum(axis=1)
                #values = values / values.sum()

                # Count occurrences of each taxonomy category
                taxonomies = station_asvs['silva_Taxon'].str.split('; ', expand=True)
                taxonomies.columns = taxa_col_names
                taxonomies = taxonomies.dropna(subset=cols_show_in_sunburst[0]).fillna('___Undetermined')[cols_show_in_sunburst]

                # get rid of the silva d__, p__, etc prefixes
                for col in taxonomies.columns:
                    taxonomies[col] = taxonomies[col].apply(lambda x: x[3:])

                # Get relative abundances
                taxonomies['abundance_values'] = values
                taxonomies = taxonomies[taxonomies['abundance_values'] != 0]

                # set title of plot
                title = '16S Silva Taxonomy, Station "' + station_id + '"'
            elif dataset == '18Sv4':
                taxa_col_names = ['Kingdom', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
                # Merge with 18Sv4 dataframe to get taxonomy data for the samples
                asv_cols = pd.Series(data_18Sv4.columns).isin(station_samples).values
                if np.sum(asv_cols) == 0: # If there are no samples for the sample_type at the station
                    # set the figure to an empty figure
                    title = f'No {sample_type} samples for station "{station_id}"'
                    fig = px.sunburst(empty_sunburst_data, path=['Phylum', 'Class', 'Order'], title=title)
                    sunburst_figs[station_id][sample_type][dataset] = fig
                    continue
                asv_cols[0] = True
                station_asvs = pd.concat([data_18Sv4.loc[:,asv_cols], data_18Sv4['pr2_Taxon']], axis=1)

                # Get relative abundances
                values = station_asvs.drop(['Feature.ID', 'pr2_Taxon'], axis=1).fillna(0).sum(axis=1)
                #values = values / values.sum()

                # Count occurrences of each taxonomy category
                taxonomies = station_asvs['pr2_Taxon'].str.split(';', expand=True)
                taxonomies = taxonomies.iloc[:, :8]
                taxonomies.columns = taxa_col_names
                taxonomies = taxonomies.dropna(subset='Phylum').fillna('Undetermined')[cols_show_in_sunburst]

                # Add relative abundances
                taxonomies['abundance_values'] = values
                taxonomies = taxonomies[taxonomies['abundance_values'] != 0]

                # set title of plot
                title = '18S v4 PR2 Taxonomy, Station "' + station_id + '"'

            elif dataset == '18Sv9':
                taxa_col_names = ['Kingdom', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

                # Merge with 18Sv9 dataframe to get taxonomy data for the samples
                asv_cols = pd.Series(data_18Sv9.columns).isin(station_samples).values
                if np.sum(asv_cols) == 0: # If there are no samples for the sample_type at the station
                    # set the figure to an empty figure
                    title = f'No {sample_type} samples for station "{station_id}"'
                    fig = px.sunburst(empty_sunburst_data, path=['Phylum', 'Class', 'Order'], title=title)
                    sunburst_figs[station_id][sample_type][dataset] = fig
                    continue
                asv_cols[0] = True
                station_asvs = pd.concat([data_18Sv9.loc[:,asv_cols], data_18Sv9['pr2_Taxon']], axis=1)

                # Get relative abundances
                values = station_asvs.drop(['Feature.ID', 'pr2_Taxon'], axis=1).fillna(0).sum(axis=1)
                #values = values / values.sum()

                # Count occurrences of each taxonomy category
                taxonomies = station_asvs['pr2_Taxon'].str.split(';', expand=True)
                taxonomies = taxonomies.iloc[:, :8]
                taxonomies.columns = taxa_col_names
                taxonomies = taxonomies.dropna(subset='Phylum').fillna('Undetermined')[cols_show_in_sunburst]

                # Add relative abundances
                taxonomies['abundance_values'] = values
                taxonomies = taxonomies[taxonomies['abundance_values'] != 0]

                # set title of plot
                title = '18S v9 PR2 Taxonomy, Station "' + station_id + '"'
            taxonomies['relative_abundance'] = taxonomies['abundance_values'] / taxonomies['abundance_values'].sum()
            taxonomies['relative_abundance'] = (taxonomies['relative_abundance'] * 100).round(2)
            fig = px.sunburst(taxonomies, path=['Phylum', 'Class', 'Order'], values='relative_abundance',
                              title=title, width=700, height=700)
            taxonomies_all[station_id][sample_type][dataset] = taxonomies
            sunburst_figs[station_id][sample_type][dataset] = fig
        counter += 3
        print('Finished: {0:.2f}%'.format(100 * counter / total_computing), end='\r')
print()

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.layout = dbc.Container([
    html.H1("Marine microbial life at California coastal stations", className="mt-4 mb-4"),
    dbc.Row([
        html.Div([
            html.Label('Sample Type:'),
            dcc.Dropdown(sample_type_vals, sample_type_vals[0], id='sample-type-dropdown')
        ])
    ]),
    dbc.Row([
        dbc.Col([
            html.Div([
                html.Label('Environmental Variable (Color):'),
                dcc.Dropdown(env_var_cols, 'NCDepth', id='env-var-dropdown')
            ], className="env-var-dropdown"),
            html.Div([
                dcc.Graph(id='map-graph')
            ], className="map-graph"),
        ], width=5),
        dbc.Col([
            html.Div([
                html.Label('Dataset:'),
                dcc.Dropdown(['16S', '18Sv4', '18Sv9'], '16S', id='dataset-dropdown')
            ]),
            html.Div([
                dcc.Graph(id='sunburst-graph')
            ]),
        ], width=5)
    ])
], fluid=True)


# Map graph dropdown callback
@app.callback(
    Output('map-graph', 'figure'),
    [Input('sample-type-dropdown', 'value'),
     Input('env-var-dropdown', 'value')]
)
def update_map(dropdown_sample_type, dropdown_env_var):
    return map_figs[dropdown_sample_type][dropdown_env_var]

# Map graph click data callback
@app.callback(
    Output('sunburst-graph', 'figure'),
    [Input('sample-type-dropdown', 'value'),
     Input('map-graph', 'clickData'),
     Input('dataset-dropdown', 'value')]
)
def update_sunburst(dropdown_sample_type, click_data, dropdown_dataset):
    if click_data is None:
        station_id = data_meta['Sta_ID'].iloc[0]
    # Get station ID from hover data
    else:
        if 'customdata' not in click_data['points'][0]:
            station_id = station_id = data_meta['Sta_ID'].iloc[0]
        else:
            station_id = click_data['points'][0]['customdata'][0]
    testing = (dropdown_sample_type, click_data, dropdown_dataset)
    return sunburst_figs[station_id][dropdown_sample_type][dropdown_dataset]

app.run_server(debug=True)
