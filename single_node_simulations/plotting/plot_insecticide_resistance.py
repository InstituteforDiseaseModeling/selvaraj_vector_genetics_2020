import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

plot_inset = 1
plot_GM = 0

exp_name = 'insecticide_resistance_single_node'
data_path = os.path.join('_data', exp_name)
fig_path = os.path.join('_data', exp_name)
os.makedirs(fig_path, exist_ok=True)

df_est = pd.read_csv(os.path.join(data_path, 'insecticide_resistance_establishment.csv'))
df_inset = pd.read_csv(os.path.join(data_path, 'insecticide_resistance_inset.csv'))

channels = ['Annual EIR']
groupby_columns = ['Mutation_Rate2']
colors = ['#2cbcb2', '#5e495a', '#d44e28', '#56445D', 'r']

if plot_inset:
    for channel in channels:
        if any(x in channel for x in channels):
            for group, gdf in df_inset.groupby(groupby_columns):
                fig = plt.figure()
                fig.set_size_inches((9.951, 6.72))
                for i, (start, startdf) in enumerate(gdf.groupby('Label')):
                    df_label = df_inset[df_inset['Label'] == start]
                    startdf = startdf[:6*365]
                    eir = sum(startdf['Annual EIR'][-365:])
                    clinical_cases = sum(startdf['Clinical Cases'])
                    time = list(range(len(startdf)))
                    plt.plot(time, startdf[channel],
                             label='%s, EIR=%i, Cases=%i' % (start, eir, clinical_cases),
                             color=colors[i], linewidth=2)
                    plt.fill_between(time,
                                     startdf[channel] - startdf[channel+'_std'],
                                     startdf[channel] + startdf[channel + '_std'],
                                     color=colors[i],
                                     linewidth=0,
                                     alpha=0.3)

                xticks = [t for t in range(0, 7*365, 365)]
                xticklabels = [labels for labels in range(7)]
                plt.xticks(ticks=xticks, labels=xticklabels)
                plt.xlim([0, 6*365])
                if 'Prevalence' in channel:
                    plt.ylim([0, 0.5])
                else:
                    plt.ylim([0, 2.5])
                plt.ylabel(channel)
                plt.xlabel('Years')
                plt.title('%s' % start)
                plt.legend(bbox_to_anchor=(1.0, -0.2), borderaxespad=0)
                plt.subplots_adjust(bottom=0.7)
                plt.tight_layout(rect=[0, 0, 0.75, 1])
                plt.savefig(os.path.join(fig_path, '%s_recessive_all.pdf' % channel),
                            bbox_inches="tight")
                plt.close('all')

if plot_GM:
    for start, startdf in df_est.groupby('Label'):
        startdf_inset = df_inset[df_inset['Label'] == start]
        for group, gdf in startdf.groupby(groupby_columns):
            gdf_inset = startdf_inset[(startdf_inset['Mutation_Rate2'] == group)
                                      ]
            eir = sum(gdf_inset['Annual EIR'][-365:])
            clinical_cases = sum(gdf_inset['Clinical Cases'])
            fig = plt.figure()
            fig.set_size_inches((9.951, 6.72))
            gdf = gdf[:6*365]
            time = list(range(len(gdf)))
            for channel in ['X-a0:X-a0', 'X-a0:X-a1', 'X-a1:X-a1']:
                plt.plot(time, gdf[channel], label=channel)
                plt.fill_between(time,
                                 gdf[channel] - gdf[channel + '_std'],
                                 gdf[channel] + gdf[channel + '_std'],
                                 linewidth=0,
                                 alpha=0.3)
            xticks = [t for t in range(0, 7*365, 365)]
            xticklabels = [labels for labels in range(7)]
            plt.xticks(ticks=xticks, labels=xticklabels)
            yticks = [0.0, 0.5, 1.0]
            yticklabels = [0.0, 0.5, 1.0]
            plt.yticks(ticks=yticks, labels=yticklabels)
            ax = plt.gca()
            majorLocator = MultipleLocator(365)
            majorFormatter = FormatStrFormatter('%d')
            minorLocator = MultipleLocator(365 / 12.0)
            ax.xaxis.set_major_locator(majorLocator)
            ax.xaxis.set_major_formatter(majorFormatter)
            ax.xaxis.set_minor_locator(minorLocator)
            plt.xlim([0, 6*365])
            plt.ylim([0, 1])
            plt.title('%s, \n EIR=%i, Cases=%i'
                      % (start, eir, clinical_cases))
            plt.ylabel('Establishment fraction')
            plt.xlabel('Years')
            plt.legend(bbox_to_anchor=(0.75, -0.2), borderaxespad=0)
            plt.subplots_adjust(bottom=0.7)
            plt.tight_layout(rect=[0, 0, 0.75, 1])
            plt.savefig(os.path.join(fig_path,
                                     "recessive_%s.pdf" % start), bbox_inches="tight")
            plt.close('all')
