import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

mpl.rcParams['pdf.fonttype'] = 42

exp_name = 'vector_genetics_single'
data_path = os.path.join('_data', exp_name)
fig_path = os.path.join('_data', exp_name)
os.makedirs(fig_path, exist_ok=True)

file = os.path.join(data_path, 'vector_genetics_data.csv')

df = pd.read_csv(file)
if 'Unnamed: 0' in df.columns:
    df = df.drop('Unnamed: 0', axis=1)

majorLocator = MultipleLocator(365)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(365/12.0)

fig, ax = plt.subplots()
colors = ['b', 'r', 'g', 'y']
columns = ['X-a0:X-a0', 'X-a0:X-a1', 'X-a1:X-a1']
for i, column in enumerate(columns):
    ax.plot(df['Time'], df[column], label=column, color=colors[i])
    ax.fill_between(df['Time'],
                    df[column] - df[column+'_std'],
                    df[column] + df[column+'_std'],
                    color=colors[i], linewidth=0,
                    alpha=0.3)
plt.legend()
plt.xlim([0, 6*365])
plt.ylim([-0.02, 1.02])
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)

# for the minor ticks, use no labels; default NullFormatter
ax.xaxis.set_minor_locator(minorLocator)
labels = [x for x in range(2017, 2027)]
ax.xaxis.set_ticklabels(labels)

fig_file = os.path.join(fig_path, 'Establishment_HW.png')
plt.savefig(fig_file)
plt.show()
