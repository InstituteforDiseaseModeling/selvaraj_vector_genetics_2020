import pandas as pd


def construct_spatial_output_df(rawdata, channel, timesteps=[]):
    n_nodes = rawdata['n_nodes']
    n_tstep = rawdata['n_tstep']
    if 'start' in rawdata:
        start = rawdata['start']
        interval = rawdata['interval']
    else:
        start, interval = 0, 1
    nodeids = rawdata['nodeids']
    data = rawdata['data']

    all_timesteps = range(start, (start + n_tstep) * interval, interval)

    df = pd.DataFrame({channel: [item for sublist in data for item in sublist],
                       'time': [item for sublist in [[x] * n_nodes for x in all_timesteps] for item in sublist],
                       'node': [item for sublist in [list(nodeids) * len(all_timesteps)] for item in sublist]})
    if not timesteps:
        return df

    timesteps = sorted(list(set(all_timesteps).intersection(timesteps)))
    return df[df['time'].isin(timesteps)]
