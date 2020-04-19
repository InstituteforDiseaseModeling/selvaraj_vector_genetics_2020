import pandas as pd
import numpy as np
import os
import json
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import rcParams
# install package pyshp
import shapefile

mpl.rcParams['pdf.fonttype'] = 42
rcParams.update({'font.size': 21})


def features_within_polygon(shapefile_path, lngmin, lngmax, latmin, latmax):
    """
    Plot man-made features like roads or natural features like rivers
    :param shapefile_path: path to shapefile describing feature
    :param lngmin: westernmost longitude of polygon to plot
    :param lngmax: easternmost longitude of polygon to plot
    :param latmin: southernmost longitude of polygon to plot
    :param latmax: northernmost longitude of polygon to plot
    :return: features: coordinates of features
    """

    sf = shapefile.Reader(shapefile_path)
    features = []
    for shape in sf.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        x2 = [i for i in x if lngmin < i < lngmax]
        y = [i[1] for i in shape.shape.points[:]]
        y2 = [i for i in y if latmin < i < latmax]
        if x2 and y2:
            features.append(shape.shape.points[:])

    return features


def get_roads_rivers(grid_data):

    # Paths to roads and rivers shape file, and to full Burkina Faso raster file
    river_shapefile_path = os.path.join('Burkina', 'Burkina shapefiles', 'BFA_wat', 'BFA_water_lines_dcw')

    road_shapefile_path = os.path.join('Burkina', 'Burkina shapefiles', 'BFA_rds', 'BFA_roads')

    lngmin = min(grid_data['lon']) - 0.02
    lngmax = max(grid_data['lon']) + 0.02
    latmin = min(grid_data['lat']) - 0.04
    latmax = max(grid_data['lat']) + 0.02
    minmaxes = (lngmin, lngmax, latmin, latmax)
    rivers = features_within_polygon(river_shapefile_path, lngmin=lngmin, lngmax=lngmax, latmax=latmax, latmin=latmin)
    roads = features_within_polygon(road_shapefile_path, lngmin=lngmin, lngmax=lngmax, latmax=latmax, latmin=latmin)

    return minmaxes, rivers, roads


def get_grid_data():

    dir = 'Burkina_Faso/Grid'
    grid_data = pd.read_csv(os.path.join(dir, 'vector_genetics_full_pop_grid.csv'))
    grid_data = grid_data[['lat', 'lon', 'node_label']]
    grid_data.rename(columns={'node_label': 'NodeID'}, inplace=True)
    grid_data['NodeID'] += 1

    direc = 'input/'
    demographics_dir = os.path.join(direc, 'VectorGeneticsSpatial', 'Burkina_Faso', 'Demographics')
    demographics_file = 'vector_genetics_demographics.json'
    demographics_filepath = os.path.join(demographics_dir, demographics_file)

    nodes = []
    pop = []
    lat = []
    lon = []
    with open(demographics_filepath) as f:
        data = json.load(f)
        for node in data['Nodes']:
            nodes.append(node['NodeID'])
            pop.append(node['NodeAttributes']['InitialPopulation'])
            lat.append(node['NodeAttributes']['Latitude'])
            lon.append(node['NodeAttributes']['Longitude'])
    df_release = pd.DataFrame({'NodeID': nodes, 'pop': pop, 'lat': lat, 'lon': lon})

    return grid_data, df_release


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def get_relevant_dataframe(df, df_nets_temp, df_no_int_temp, transmission, net_cov):

    df_GM_temp = df.loc[(df['Transmission_To_Human'] == transmission) &
                        (df['ITN_Coverage'] == 0.0)
                        ]
    df_VC_GM_temp = df.loc[(df['Transmission_To_Human'] == transmission) &
                           (df['ITN_Coverage'] == net_cov)
                           ]
    df_temp = pd.concat([df_GM_temp, df_VC_GM_temp])
    df_temp = pd.concat([df_temp, df_nets_temp])
    df_temp = pd.concat([df_temp, df_no_int_temp])

    return df_temp


def plotting_shapefiles(times, df, cmap):

    grid_data, df_release_master = get_grid_data()
    minmaxes, rivers, roads = get_roads_rivers(grid_data)

    transmission_to_humans = list(df['Transmission_To_Human'].unique())
    transmission_to_humans.remove(1.0)
    itn_coverages = list(df['ITN_Coverage'].unique())
    itn_coverages.remove(0.0)

    for nc, net_cov in enumerate(itn_coverages):
        df_nets_temp = df.loc[(df['Transmission_To_Human'] == 1.0) & (df['ITN_Coverage'] == net_cov)]
        df_no_int_temp = df.loc[(df['Transmission_To_Human'] == 1.0) & (df['ITN_Coverage'] == 0.0)]
        df_release = df_release_master.sort_values('pop', ascending=False)[:6]
        for tr, transmission in enumerate(transmission_to_humans):

            df_temp = get_relevant_dataframe(df, df_nets_temp, df_no_int_temp, transmission, net_cov)

            for t, time in enumerate(times):

                fig = plt.figure(t)
                fig.set_size_inches((28, 7))

                expts = ['No_interventions', 'VC_only', 'GM_only', 'VC_and_GM']

                for j, exp in enumerate(expts):

                    df_frac = df_temp.loc[(df_temp['Experiment'] == exp) & (df_temp['Time']==time)]

                    axes = [fig.add_subplot(1, len(expts), j + 1)]
                    ax1 = axes[0]

                    for river in rivers:
                        x = [i[0] for i in river]
                        y = [i[1] for i in river]
                        ax1.plot(x, y, zorder=1, color='#7faddd')

                    for road in roads:
                        x = [i[0] for i in road]
                        y = [i[1] for i in road]
                        ax1.plot(x, y, zorder=1, color='#969696')

                    ax1.set_ylim([minmaxes[2], minmaxes[3]])
                    ax1.set_xlim([minmaxes[0], minmaxes[1]])

                    ##############################################################################################

                    norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
                    m = cm.ScalarMappable(norm=norm, cmap=cmap)

                    df_frac = df_frac[['Time', 'NodeID', 'GM', 'Experiment']]
                    df_frac = pd.merge(df_frac, df_release_master, on='NodeID')
                    df_frac['Colors'] = df_frac['GM'].apply(lambda x: m.to_rgba(x))
                    ax1.scatter(df_frac['lon'], df_frac['lat'], zorder=2, marker='s', s=100, c=df_frac['Colors'],
                                edgecolor='silver')
                    if time == 180:
                        sc = ax1.scatter(df_release['lon'], df_release['lat'], zorder=3, marker='s', s=100, c='k',
                                         edgecolor='r')
                    else:
                        sc = ax1.scatter(df_release['lon'], df_release['lat'], zorder=3, marker='s', s=100, c='k',
                                         edgecolor='#7FFF00')
                    sc.set_facecolor("none")

                    # plot scale
                    ax1.plot([-1.57267, -1.52673], [12.0, 12.0], 'k', lw=2)
                    ax1.text(-1.57267, 12.005, r'5 km', fontsize=21)
                    ax1.set_aspect('equal', adjustable='box')
                    ax1.axes.get_xaxis().set_visible(False)
                    ax1.axes.get_yaxis().set_visible(False)

                    if exp == 'VC_and_GM':
                        m.set_array(np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]))
                        cax = fig.add_axes([0.82, 0.54, 0.02, 0.35])
                        cbar = fig.colorbar(m, cax=cax, ticks=[0.0, 0.5, 1.0])
                        cbar.ax.tick_params(labelsize=21)
                        plt.text(5.8, 0.5, 'GM', verticalalignment='center', horizontalalignment='center',
                                 fontsize=21)

                    if exp == 'No_interventions':
                        ax1.set_title('No interventions', fontsize=21)
                    elif exp == 'VC_only':
                        ax1.set_title('VC', fontsize=21)
                    elif exp == 'GM_only':
                        ax1.set_title('GM', fontsize=21)
                    else:
                        ax1.set_title('VC and GM', fontsize=21)

                plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9, left=0.15, wspace=0.25)

                fig.suptitle('Netcov = %0.1f, Tran = %0.2f'
                                         % (net_cov, transmission))
                fig_file = os.path.join(fig_dir,  'Spatial_GM_%i.pdf'
                                         % time)
                plt.savefig(fig_file, type='pdf')
                plt.close('all')


if __name__ == '__main__':

    exp_name = 'Establishment_Burkina'
    channel = 'GM'

    data_dir = '/_data'
    fig_dir = '/_figures'

    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    times = [180, 1094, 2189]

    # File from analyzers
    file_reduced = os.path.join(data_dir, 'simulation_data.csv')

    df = pd.read_csv(file_reduced)
    if 'Unnamed: 0' in df.columns:
        df = df.drop('Unnamed: 0', axis=1)
    if 'Unnamed: 0.1' in df.columns:
        df = df.drop('Unnamed: 0.1', axis=1)

    #     # Color map
    cmap = cm.plasma
    cmap = truncate_colormap(cmap, 0.0, 1.0)
    plotting_shapefiles(times, df, cmap)
