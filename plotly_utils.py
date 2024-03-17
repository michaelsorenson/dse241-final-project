import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
import os
import pandas as pd
import numpy as np
import plotly.express as px
from dotenv import load_dotenv
load_dotenv()
px.set_mapbox_access_token(os.getenv('MAPBOX_ACCESS_TOKEN'))

def get_16S_top_phyla_colors(data_16S):
    taxa_16S = data_16S['silva_Taxon'].str.split(';', expand=True)
    taxa_16S.columns = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    abundance_cols = pd.Series(data_16S.columns).apply(lambda x: True if x[:3] == 'X20' else False).values
    abundance_df = data_16S.loc[:,abundance_cols]
    taxa_16S['abundance_values'] = (abundance_df / abundance_df.sum(axis=0)).sum(axis=1)
    taxa_16S = taxa_16S.dropna(subset='Phylum')[['Phylum', 'abundance_values']]
    taxa_16S['Phylum'] = taxa_16S['Phylum'].apply(lambda x: x[4:])
    top_10_phyla_16S = (
        taxa_16S.groupby('Phylum')['abundance_values'].sum()
        .sort_values(ascending=False).index.to_list()[:10]
    )

    colors_16S = {
        top_10_phyla_16S[i]: px.colors.qualitative.Plotly[i] for i in range(10)
    }
    colors_16S['Undetermined'] = '#C0C0C0' # light grey
    colors_16S['Other'] = '#71797E' # dark grey

    return (top_10_phyla_16S, colors_16S)


def get_18S_top_phyla_colors(data_18Sv4, data_18Sv9):
    # show these classes separately from their parent phyla
    special_18S_classes = ['Syndiniales', 'Bacillariophyta']
    # 18Sv4 ---------------------------------------------------------------------
    taxa_18Sv4 = data_18Sv4['pr2_Taxon'].str.split(';', expand=True).iloc[:, :8]
    taxa_18Sv4.columns = ['Kingdom', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    abundance_cols = pd.Series(data_18Sv4.columns).apply(lambda x: True if x[:3] == 'X20' else False).values
    abundance_df = data_18Sv4.loc[:,abundance_cols]
    taxa_18Sv4['abundance_values'] = (abundance_df / abundance_df.sum(axis=0)).sum(axis=1)
    taxa_18Sv4 = taxa_18Sv4.dropna(subset='Phylum')[['Phylum', 'abundance_values']]
    top_8_phyla_18Sv4 = (
        taxa_18Sv4.groupby('Phylum')['abundance_values'].sum()
        .sort_values(ascending=False).index.to_list()[:8]
    )

    colors_18S = {
        top_8_phyla_18Sv4[i]: px.colors.qualitative.Plotly[i] for i in range(8)
    }
    for i, special_class in enumerate(special_18S_classes):
        colors_18S[special_class] = px.colors.qualitative.Plotly[i + 8]
    colors_18S['Undetermined'] = '#C0C0C0' # light grey
    colors_18S['Other'] = '#71797E' # dark grey

    # 18Sv9 ---------------------------------------------------------------------
    taxa_18Sv9 = data_18Sv9['pr2_Taxon'].str.split(';', expand=True).iloc[:, :8]
    taxa_18Sv9.columns = ['Kingdom', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    abundance_cols = pd.Series(data_18Sv9.columns).apply(lambda x: True if x[:3] == 'X20' else False).values
    abundance_df = data_18Sv9.loc[:,abundance_cols]
    taxa_18Sv9['abundance_values'] = (abundance_df / abundance_df.sum(axis=0)).sum(axis=1)
    taxa_18Sv9 = taxa_18Sv9.dropna(subset='Phylum')[['Phylum', 'abundance_values']]
    top_8_phyla_18Sv9 = (
        taxa_18Sv9.groupby('Phylum')['abundance_values'].sum()
        .sort_values(ascending=False).index.to_list()[:8]
    )

    # colors_18Sv9 = {
    #     top_8_phyla_18Sv9[i]: px.colors.qualitative.Plotly[i] for i in range(8)
    # }
    # for i, (phylum, special_class) in enumerate(special_18S_classes):
    #     colors_18Sv9[special_class] = px.colors.qualitative.Plotly[i + 8]
    # colors_18Sv9['Undetermined'] = '#C0C0C0' # light grey
    # colors_18Sv9['Other'] = '#71797E' # dark grey
    if set(top_8_phyla_18Sv4) != set(top_8_phyla_18Sv9):
        print('Error: v4 and v9 do not have the same top 8 phyla. Are you sure the dataset is correct?')
    # colors_18Sv9 = colors_18Sv4 # the top 8 groups should match but I'll leave the code here just in case

    return (top_8_phyla_18Sv4, top_8_phyla_18Sv9, colors_18S)


def precompute_map_figs(data_meta):
    # Lat and Lon center for plots
    cal_coast_center = dict(
        lat=np.mean([min(data_meta['Lat_Dec']), max(data_meta['Lat_Dec'])]),
        lon=np.mean([min(data_meta['Lon_Dec']), max(data_meta['Lon_Dec'])])
    )

    # Environmental variables and sample types
    env_var_cols = ['T_degC', 'Salnty', 'O2ml_L', 'PO4ug', 'SiO3ug', 'NO3ug', 'NH3ug', 'ChlorA', 'IntC14', 'NCDepth']
    sample_type_vals = data_meta['sample_type'].dropna().unique()

    map_figs = {sample_type: {} for sample_type in sample_type_vals}
    for sample_type in sample_type_vals:
        for env_var in env_var_cols:
            meta_subset = data_meta[data_meta['sample_type'] == sample_type].copy()
            meta_subset = meta_subset[['Sta_ID', 'Lat_Dec', 'Lon_Dec', 'sample_type', env_var]].groupby('Sta_ID').agg({
                'Lat_Dec': 'min',
                'Lon_Dec': 'min',
                'sample_type': 'count',
                env_var: 'mean'
            }).rename({'sample_type': 'Number of Samples'}, axis=1).reset_index()
            hover_names = meta_subset['Sta_ID'].apply(lambda x: '<b>Station: </b>' + x)
            if env_var == 'NCDepth':
                color_continuous_scale='viridis_r'
            else:
                color_continuous_scale='viridis'
            subset_fig = px.scatter_mapbox(meta_subset, lat='Lat_Dec', lon='Lon_Dec', center=cal_coast_center,
                                        color=env_var, hover_name=hover_names, hover_data='Number of Samples', #size="num_samples",
                                        color_continuous_scale=color_continuous_scale, size_max=15, zoom=4.5, mapbox_style='outdoors',
                                        width=600, height=700, custom_data='Sta_ID')
            map_figs[sample_type][env_var] = subset_fig

        # Create a figure where environmental variable is "Total number of samples"
        meta_subset = data_meta[data_meta['sample_type'] == sample_type].copy()
        meta_subset = meta_subset[['Sta_ID', 'Lat_Dec', 'Lon_Dec', 'sample_type']].groupby('Sta_ID').agg({
            'Lat_Dec': 'min',
            'Lon_Dec': 'min',
            'sample_type': 'count'
        }).rename({'sample_type': 'Total Number of Samples'}, axis=1).reset_index()
        hover_names = meta_subset['Sta_ID'].apply(lambda x: '<b>Station: </b>' + x)
        subset_fig = px.scatter_mapbox(meta_subset, lat='Lat_Dec', lon='Lon_Dec', center=cal_coast_center,
                                    color='Total Number of Samples', hover_name=hover_names,
                                    color_continuous_scale='viridis', size_max=15, zoom=4.5, mapbox_style='outdoors',
                                    width=600, height=700, custom_data='Sta_ID')
        map_figs[sample_type]['Total Number of Samples'] = subset_fig
    return map_figs


def precompute_sunburst_figs(data_meta, data_16S, data_18Sv4, data_18Sv9):
    # empty sunburst data for if the station has no data with the selected sample type
    empty_sunburst_data = {'Phylum': [], 'Class': [], 'Order': []}
    for parent in ['Undetermined_1', 'Undetermined_2', 'Undetermined_3']:
        for child1 in ['Undetermined_1', 'Undetermined_2', 'Undetermined_3']:
            for child2 in ['Undetermined_1', 'Undetermined_2', 'Undetermined_3']:
                empty_sunburst_data['Phylum'].append(parent)
                empty_sunburst_data['Class'].append(child1)
                empty_sunburst_data['Order'].append(child2)

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
                    title = '16S Silva Taxonomy<br>(innermost) Phylum->Class->Order (outermost)<br>Station "' + station_id + '"'
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
                    title = '18S v4 PR2 Taxonomy<br>(innermost) Phylum->Class->Order (outermost)<br>Station "' + station_id + '"'

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
                    title = '18S v9 PR2 Taxonomy<br>(innermost) Phylum->Class->Order (outermost)<br>Station "' + station_id + '"'
                taxonomies['relative_abundance'] = taxonomies['abundance_values'] / taxonomies['abundance_values'].sum()
                taxonomies['relative_abundance'] = (taxonomies['relative_abundance'] * 100).round(2)
                fig = px.sunburst(taxonomies, path=['Phylum', 'Class', 'Order'], values='relative_abundance',
                                  title=title, width=700, height=700)
                fig.data[0].hovertemplate = '<b>%{label}</b><br>Relative Abundance: %{value}%'
                taxonomies_all[station_id][sample_type][dataset] = taxonomies
                sunburst_figs[station_id][sample_type][dataset] = fig
            counter += 3
            print('Finished: {0:.2f}%'.format(100 * counter / total_computing), end='\r')
    print()
    return sunburst_figs


def precompute_stacked_bar_figs(data_meta, data_16S, data_18Sv4, data_18Sv9):
    top_10_phyla_16S, colors_16S = get_16S_top_phyla_colors(data_16S)
    top_8_phyla_18Sv4, top_8_phyla_18Sv9, colors_18S = get_18S_top_phyla_colors(data_18Sv4, data_18Sv9)
    counter = 0
    total_computing = len(data_meta['Sta_ID'].dropna().unique())  * 3
    print(f'Pre-computing {total_computing} Stacked Bar Figures')
    stacked_bar_figs = {station_id: {} for station_id in data_meta['Sta_ID'].dropna().unique()}
    for station_id in data_meta['Sta_ID'].dropna().unique():
        for dataset in ['16S', '18Sv4', '18Sv9']:
            station_data = data_meta[(data_meta['Sta_ID'] == station_id)]
            station_samples = station_data['sampleid'].tolist()
            if dataset == '16S':
                taxa_col_names = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
                asv_cols = pd.Series(data_16S.columns).isin(station_samples).values
                if np.sum(asv_cols) == 0:
                    empty_fig = go.Figure(
                        data=[
                            go.Bar(name='Group 1', x=['No Samples'], y=[1], offsetgroup=1),
                            go.Bar(name='Group 2', x=['No Samples'], y=[1], offsetgroup=1, base=[1])
                        ],
                        layout=go.Layout(
                            title=f'Station {station_id} phyla distribution per sample',
                            yaxis_title='Relative Abundance', height=600, width=600
                        )
                    )
                    stacked_bar_figs[station_id][dataset] = empty_fig
                    continue
                
                asv_col_names = pd.Series(data_16S.columns[asv_cols])
                asv_col_sample_types = pd.Series([
                    station_data[station_data['sampleid'] == sampleid]['sample_type'].iloc[0]
                        for sampleid in asv_col_names
                ])
                asv_cols[0] = True
                station_asvs = pd.concat([data_16S.loc[:,asv_cols], data_16S['silva_Taxon']], axis=1)

                # Get relative abundances
                values = station_asvs.drop(['Feature.ID', 'silva_Taxon'], axis=1).fillna(0)
                values = values / values.sum()

                # Count occurrences of each taxonomy category
                taxonomies = station_asvs['silva_Taxon'].str.split('; ', expand=True)
                taxonomies.columns = taxa_col_names
                taxonomies = taxonomies.fillna('___Undetermined')[['Phylum']]

                # get rid of the silva d__, p__, etc prefixes
                taxonomies['Phylum'] = taxonomies['Phylum'].apply(lambda x: x[3:])

                # Get relative abundances of top phyla / special classes
                taxonomies = pd.concat([taxonomies, values], axis=1)
                phyla_counts = taxonomies.groupby('Phylum').sum()
                for phylum in top_10_phyla_16S:
                    if phylum not in phyla_counts.index:
                        phyla_counts.loc[phylum] = 0
                top_phyla_counts = phyla_counts.loc[top_10_phyla_16S]
                other_phyla_counts = phyla_counts[(
                    (phyla_counts.index != 'Undetermined') & (~phyla_counts.index.isin(top_10_phyla_16S))
                )].sum()
                top_phyla_counts.loc['Other'] = other_phyla_counts
                top_phyla_counts.loc['Undetermined'] = phyla_counts.loc['Undetermined']

                top_phyla_counts.columns = asv_col_sample_types + '_' + asv_col_names
                top_phyla_counts = top_phyla_counts[top_phyla_counts.columns.sort_values()]

                # Create stacked bar figure
                base = top_phyla_counts.loc[top_phyla_counts.index[0]].copy()
                stacked_bar_labels = top_phyla_counts.columns
                fig_data = [go.Bar(
                    name=top_phyla_counts.index[0],
                    x=stacked_bar_labels,
                    y=base,
                    customdata=base.copy(),
                    marker=go.bar.Marker(color=colors_16S[top_phyla_counts.index[0]]),
                    hovertemplate='Rel. Abundance:%{customdata:.3f}',
                    offsetgroup=1,
                )]
                for phyla in top_phyla_counts.index[1:]:
                    vals = top_phyla_counts.loc[phyla]
                    fig_data.append(go.Bar(
                        name=phyla,
                        x=stacked_bar_labels,
                        y=vals,
                        marker=go.bar.Marker(color=colors_16S[phyla]),
                        customdata=vals.copy(),
                        hovertemplate='Rel. Abundance:%{customdata:.3f}',
                        offsetgroup=1,
                        base=base
                    ))
                    base += vals
                stacked_bar_fig = go.Figure(
                    data=fig_data,
                    layout=go.Layout(
                        title=f'Phyla Distribution Per Sample, Station "{station_id}"',
                        yaxis_title='Relative abundance'
                    )
                )
                stacked_bar_figs[station_id][dataset] =  stacked_bar_fig
            # End 16S -----------------------------------------------------------------------------------------------

            if dataset == '18Sv4':
                taxa_col_names = ['Kingdom', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
                asv_cols = pd.Series(data_18Sv4.columns).isin(station_samples).values
                if np.sum(asv_cols) == 0:
                    empty_fig = go.Figure(
                        data=[
                            go.Bar(name='Group 1', x=['No Samples'], y=[1], offsetgroup=1),
                            go.Bar(name='Group 2', x=['No Samples'], y=[1], offsetgroup=1, base=[1])
                        ],
                        layout=go.Layout(
                            title=f'Station {station_id} phyla distribution per sample',
                            yaxis_title='Relative Abundance', height=600, width=600
                        )
                    )
                    stacked_bar_figs[station_id][dataset] = empty_fig
                    continue
                
                asv_col_names = pd.Series(data_18Sv4.columns[asv_cols])
                asv_col_sample_types = pd.Series([
                    station_data[station_data['sampleid'] == sampleid]['sample_type'].iloc[0]
                        for sampleid in asv_col_names
                ])
                asv_cols[0] = True
                station_asvs = pd.concat([data_18Sv4.loc[:,asv_cols], data_18Sv4['pr2_Taxon']], axis=1)

                # Get relative abundances
                values = station_asvs.drop(['Feature.ID', 'pr2_Taxon'], axis=1).fillna(0)
                values = values / values.sum()

                # Count occurrences of each taxonomy category
                taxonomies = station_asvs['pr2_Taxon'].str.split(';', expand=True).iloc[:,:8]
                taxonomies.columns = taxa_col_names
                taxonomies = taxonomies.fillna('Undetermined')[['Phylum', 'Class']]

                # Get relative abundances of top phyla / special classes
                taxonomies = pd.concat([taxonomies, values], axis=1)
                phyla_counts = taxonomies.drop('Class', axis=1).groupby('Phylum').sum()
                class_counts = taxonomies.drop('Phylum', axis=1).groupby('Class').sum()
                if 'Syndiniales' not in class_counts.index:
                    class_counts.loc['Syndiniales'] = 0
                if 'Bacillariophyta' not in class_counts.index:
                    class_counts.loc['Bacillariophyta'] = 0
                for phylum in top_8_phyla_18Sv4 + ['Dinoflagellata', 'Ochrophyta']:
                    if not phylum in phyla_counts.index:
                        phyla_counts.loc[phylum] = 0
                phyla_counts.loc['Dinoflagellata'] -= class_counts.loc['Syndiniales']
                phyla_counts.loc['Ochrophyta'] -= class_counts.loc['Bacillariophyta']

                class_counts.index.names = ['Phylum']
                # Top 10 phyla + special classes
                top_phyla_counts = pd.concat([
                    phyla_counts.loc[top_8_phyla_18Sv4],
                    class_counts.loc[['Syndiniales', 'Bacillariophyta']]
                ])
                other_phyla_counts = phyla_counts.loc[(
                    (phyla_counts.index != 'Undetermined') & (~phyla_counts.index.isin(top_8_phyla_18Sv4))
                )].sum()
                top_phyla_counts.loc['Other'] = other_phyla_counts
                top_phyla_counts.loc['Undetermined'] = phyla_counts.loc['Undetermined']
                
                top_phyla_counts.columns = asv_col_sample_types + '_' + asv_col_names
                top_phyla_counts = top_phyla_counts[top_phyla_counts.columns.sort_values()]
                
                # Create stacked bar figure
                base = top_phyla_counts.loc[top_phyla_counts.index[0]].copy()
                stacked_bar_labels = top_phyla_counts.columns

                fig_data = [go.Bar(
                    name=top_phyla_counts.index[0],
                    x=stacked_bar_labels,
                    y=base,
                    customdata=base.copy(),
                    marker=go.bar.Marker(color=colors_18S[top_phyla_counts.index[0]]),
                    hovertemplate='Rel. Abundance:%{customdata:.3f}',
                    offsetgroup=1,
                )]
                for phyla in top_phyla_counts.index[1:]:
                    vals = top_phyla_counts.loc[phyla]
                    fig_data.append(go.Bar(
                        name=phyla,
                        x=stacked_bar_labels,
                        y=vals,
                        marker=go.bar.Marker(color=colors_18S[phyla]),
                        customdata=vals.copy(),
                        hovertemplate='Rel. Abundance:%{customdata:.3f}',
                        offsetgroup=1,
                        base=base
                    ))
                    base += vals
                stacked_bar_fig = go.Figure(
                    data=fig_data,
                    layout=go.Layout(
                        title=f'Phyla Distribution Per Sample, Station "{station_id}"',
                        yaxis_title='Relative abundance'
                    )
                )
                stacked_bar_figs[station_id][dataset] =  stacked_bar_fig
            # End 18Sv4 ---------------------------------------------------------------------------------------------

            if dataset == '18Sv9':
                taxa_col_names = ['Kingdom', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
                asv_cols = pd.Series(data_18Sv9.columns).isin(station_samples).values
                if np.sum(asv_cols) == 0:
                    empty_fig = go.Figure(
                        data=[
                            go.Bar(name='Group 1', x=['No Samples'], y=[1], offsetgroup=1),
                            go.Bar(name='Group 2', x=['No Samples'], y=[1], offsetgroup=1, base=[1])
                        ],
                        layout=go.Layout(
                            title=f'Station {station_id} phyla distribution per sample',
                            yaxis_title='Relative Abundance', height=600, width=600
                        )
                    )
                    stacked_bar_figs[station_id][dataset] = empty_fig
                    continue
                asv_col_names = pd.Series(data_18Sv9.columns[asv_cols])
                asv_col_sample_types = pd.Series([
                    station_data[station_data['sampleid'] == sampleid]['sample_type'].iloc[0]
                        for sampleid in asv_col_names
                ])
                asv_cols[0] = True
                station_asvs = pd.concat([data_18Sv9.loc[:,asv_cols], data_18Sv9['pr2_Taxon']], axis=1)
                
                # Get relative abundances
                values = station_asvs.drop(['Feature.ID', 'pr2_Taxon'], axis=1).fillna(0)
                values = values / values.sum()

                # Count occurrences of each taxonomy category
                taxonomies = station_asvs['pr2_Taxon'].str.split(';', expand=True).iloc[:,:8]
                taxonomies.columns = taxa_col_names
                taxonomies = taxonomies.fillna('Undetermined')[['Phylum', 'Class']]

                # Get relative abundances of top phyla / special classes
                taxonomies = pd.concat([taxonomies, values], axis=1)
                phyla_counts = taxonomies.drop('Class', axis=1).groupby('Phylum').sum()
                class_counts = taxonomies.drop('Phylum', axis=1).groupby('Class').sum()
                for phylum in top_8_phyla_18Sv9 + ['Dinoflagellata', 'Ochrophyta']:
                    if phylum not in phyla_counts.index:
                        phyla_counts.loc[phylum] = 0
                for class_name in ['Syndiniales', 'Bacillariophyta']:
                    if class_name not in class_counts.index:
                        class_counts.loc[class_name] = 0
                phyla_counts.loc['Dinoflagellata'] -= class_counts.loc['Syndiniales']
                phyla_counts.loc['Ochrophyta'] -= class_counts.loc['Bacillariophyta']
                
                class_counts.index.names = ['Phylum']
                # Top 10 phyla + special classes
                top_phyla_counts = pd.concat([
                    phyla_counts.loc[top_8_phyla_18Sv9],
                    class_counts.loc[['Syndiniales', 'Bacillariophyta']]
                ])
                other_phyla_counts = phyla_counts[(
                    (phyla_counts.index != 'Undetermined') & (~phyla_counts.index.isin(top_8_phyla_18Sv9))
                )].sum()
                top_phyla_counts.loc['Other'] = other_phyla_counts
                top_phyla_counts.loc['Undetermined'] = phyla_counts.loc['Undetermined']

                if not all(top_phyla_counts.columns == asv_col_names):
                    print('HERE FKK')
                
                top_phyla_counts.columns = asv_col_sample_types + '_' + asv_col_names
                top_phyla_counts = top_phyla_counts[top_phyla_counts.columns.sort_values()]
                
                # Create stacked bar figure
                base = top_phyla_counts.loc[top_phyla_counts.index[0]].copy()
                stacked_bar_labels = top_phyla_counts.columns

                # Create stacked bar figure
                base = top_phyla_counts.loc[top_phyla_counts.index[0]].copy()
                stacked_bar_labels = asv_col_names + asv_col_sample_types
                fig_data = [go.Bar(
                    name=top_phyla_counts.index[0],
                    x=stacked_bar_labels,
                    y=base,
                    customdata=base.copy(),
                    marker=go.bar.Marker(color=colors_18S[top_phyla_counts.index[0]]),
                    hovertemplate='Rel. Abundance:%{customdata:.3f}',
                    offsetgroup=1,
                )]
                for phyla in top_phyla_counts.index[1:]:
                    vals = top_phyla_counts.loc[phyla]
                    fig_data.append(go.Bar(
                        name=phyla,
                        x=stacked_bar_labels,
                        y=vals,
                        marker=go.bar.Marker(color=colors_18S[phyla]),
                        customdata=vals.copy(),
                        hovertemplate='Rel. Abundance:%{customdata:.3f}',
                        offsetgroup=1,
                        base=base
                    ))
                    base += vals
                stacked_bar_fig = go.Figure(
                    data=fig_data,
                    layout=go.Layout(
                        title=f'Phyla Distribution Per Sample, Station "{station_id}"',
                        yaxis_title='Relative abundance'
                    )
                )
                stacked_bar_figs[station_id][dataset] =  stacked_bar_fig
            # End 18Sv9 ---------------------------------------------------------------------------------------------
        # Finished all three datasets
        counter += 3
        print('Finished: {0:.2f}%'.format(100 * counter / total_computing), end='\r')
    return stacked_bar_figs