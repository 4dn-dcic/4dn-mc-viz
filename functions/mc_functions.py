import csv
import gzip
import numpy as np
#from matplotlib import pyplot as plt
import time
import urllib.request as use_urllib
import h5py
import pandas as pd
#from dcicutils import jh_utils


def get_chrom_sizes(accession):
    chromsizes_metadata = jh_utils.find_valid_file_or_extra_file(accession, format=None)

    with use_urllib.urlopen(chromsizes_metadata['full_href']) as f:
        chrom_sizes = [(row[0], row[1]) for row in [
            item.split('\t') for item in f.read().decode('utf-8').split('\n') if item
        ]]
    return chrom_sizes


def get_bins(chrom_sizes, binsize):
    bins_list = [(i, j, j + binsize) if int(item[1]) > j + binsize else (i, j, int(item[1]))
                 for i, item in enumerate(chrom_sizes) for j in range(0, int(item[1]), binsize)]
    return {bin_tuple[:2]: bins_list.index(bin_tuple) for bin_tuple in bins_list}


def create_clusters_h5(clusters_acc, binsize, chromsize_acc, outfile, lines=10000000):
    pwr = np.log10(binsize)
    if pwr != int(pwr):
        raise ValueError('binsize parameter should be a power of 10')
        return
    print('getting chrom sizes')
    chromsizes = get_chrom_sizes(chromsize_acc)
    print('getting bins')
    bin_dict = get_bins(chromsizes, binsize)
    file_info = jh_utils.find_valid_file_or_extra_file(clusters_acc, format=None)
    with use_urllib.urlopen(file_info['full_href']) as f:
        with gzip.open(f, 'rt') as infile:
            print('parsing clusters')
            chrom_dict = {item[0]: i for i, item in enumerate(chromsizes)}
            cluster_list = []
            for i in range(lines):
                line = [item for item in infile.readline().split()
                        if 'random' not in item and 'chrUn' not in item
                        and 'chrM' not in item and 'chrEBV' not in item]
                locs = list(set(
                    (item.split(':')[0], item.split(':')[1][:-1*int(pwr)] + str(binsize)[1:]) for item in line[1:]
                ))
                if len(locs) < 2:
                    continue
                for loc in locs:
                    bin_id = bin_dict[(chrom_dict[loc[0]], int(loc[1]))]
                    cluster_list.append((bin_id, i, len(locs)))
    cluster_list = sorted(cluster_list)
    with h5py.File(outfile, 'w') as f:
        print('writing h5 file')
        print('writing chroms')
        chroms = f.create_group("chroms")
        name = chroms.create_dataset("name", data=np.array([item[0] for item in chromsizes], dtype='S64'))
        length = chroms.create_dataset("length", data=np.array([item[1] for item in chromsizes], dtype=np.int32))
        bins = f.create_group("bins")
        print('writing bins')
        bin_list = [
            (i, j, j + binsize) if int(item[1]) > j + binsize else (i, j, int(item[1]))
            for i, item in enumerate(chromsizes) for j in range(0, int(item[1]), binsize)
        ]
        chrom = bins.create_dataset("chrom", data=np.array([item[0] for item in bin_list]), dtype=np.int32)
        start = bins.create_dataset("start", data=np.array([item[1] for item in bin_list]), dtype=np.int32)
        end = bins.create_dataset("end", data=np.array([item[2] for item in bin_list]), dtype=np.int32)
        clusters = f.create_group("clusters")
        clusters.attrs['max'] = max([item[1] for item in cluster_list])
        print('writing clusters')
        cluster_bins = [item[0] for item in cluster_list]
        cl_bin = clusters.create_dataset("bin", data=np.array(cluster_bins), dtype=np.int32)
        cl_name = clusters.create_dataset(
            "cluster_name", data=np.array([item[1] for item in cluster_list], dtype=np.int32)
        )
        cl_len = clusters.create_dataset(
            "cluster_length", data=np.array([item[2] for item in cluster_list]), dtype=np.int32
        )
        print('writing offsets')
        bin_offset_list = np.zeros(len(bin_list) + 1, dtype=np.int32)
        curr_val = 0
        for start, length, value in zip(*rlencode(cluster_bins, 1000000)):
            bin_offset_list[curr_val:value + 1] = start
            curr_val = value + 1
        bin_offset_list[curr_val:] = len(cluster_bins)
        idx = f.create_group('index')
        bin_offset = idx.create_dataset('bin_offset', data=np.array(bin_offset_list), dtype=np.int32)


def rlencode(array, chunksize=None):
    """
    Run length encoding.
    Based on http://stackoverflow.com/a/32681075, which is based on the rle
    function from R.
    Parameters
    ----------
    x : 1D array_like
        Input array to encode
    dropna: bool, optional
        Drop all runs of NaNs.
    Returns
    -------
    start positions, run lengths, run values
    """
    where = np.flatnonzero
    array = np.asarray(array)
    n = len(array)
    if n == 0:
        return (np.array([], dtype=int),
                np.array([], dtype=int),
                np.array([], dtype=array.dtype))

    if chunksize is None:
        chunksize = n

    starts, values = [], []
    last_val = np.nan
    for i in range(0, n, chunksize):
        x = array[i:i+chunksize]
        locs = where(x[1:] != x[:-1]) + 1
        if x[0] != last_val:
            locs = np.r_[0, locs]
        starts.append(i + locs)
        values.append(x[locs])
        last_val = x[-1]
    starts = np.concatenate(starts)
    lengths = np.diff(np.r_[starts, n])
    values = np.concatenate(values)

    return starts, lengths, values


def filter_cluster_data(filename, bins):

    with h5py.File(filename, 'r') as f:
        data = {'bin': f['clusters']['bin'], 'len': f['clusters']['cluster_name']}
        df = pd.DataFrame(data)
        offsets = np.array(f['index']['bin_offset'][bins[0]:bins[0]+2])
        cluster_max = f['clusters'].attrs['max'] + 1

    idx = np.zeros(cluster_max, dtype='bool')
    clusters = df['len'].iloc[offsets[0]:offsets[1]]
    idx[clusters] = True
    df = df.iloc[idx[df['len'].values]]
    for b in bins[1:]:
        idx = np.zeros(cluster_max, dtype='bool')
        clusters = df['len'][df['bin'] == b]
        idx[clusters] = True
        df = df[idx[df['len'].values]]

    return df


# def generate_dist_chart_h5(filename, mode='counts'):
#     with h5py.File(filename, 'r') as f:
#         data = {'bin': f['clusters']['bin'], 'len': f['clusters']['cluster_length']}
#         df = pd.DataFrame(data)
#     nums = pd.Series(range(0, 31000))
#     df2 = pd.DataFrame(data={'idx': nums})
#     df2['2-10'] = df[df['len'] < 11].groupby('bin').count()
#     df2['11-100'] = df[(df['len'] > 10) & (df['len'] < 101)].groupby('bin').count()
#     df2['101-1000'] = df[(df['len'] > 100) & (df['len'] < 1001)].groupby('bin').count()
#     df2['1000+'] = df[df['len'] > 1000].groupby('bin').count()
#     df2 = df2.fillna(0).astype(int)
#
#     subset = df2.loc[:2489]
#     if mode == 'counts':
#         y = np.row_stack(tuple(list(subset.iloc[:,i]) for i in range(1,5)))
#
#     elif mode == 'pct':
#         cols = ['2-10', '11-100', '101-1000', '1000+']
#         data2 = subset[cols].div(subset[cols].sum(axis=1), axis=0).multiply(100)
#         y = np.row_stack(tuple(list(data2.iloc[:,i]) for i in range(0,4)))
#
#     x = list(subset['idx'])
#     y_stack = np.cumsum(y, axis=0)
#
#     fig = plt.figure(figsize=(18,4))
#     ax1 = fig.add_subplot(111)
#     # if not colors:
#     colors = ["#DDA0DD", "#1DACD6", "#6E5160", "#32CD32", "#191970", "#FFA600"]
#
#     ax1.fill_between(x, 0, y_stack[0,:], facecolor=colors[0], alpha=.7)
#     for i in range(1, 5-1):
#         ax1.fill_between(x, y_stack[i-1,:], y_stack[i,:], facecolor=colors[i])
#
#     #ticks = x[::int(len(x)/12)+1]
#     #labs = [t[:-6] + 'Mb' if t.endswith('000000') else t[:-3] + 'kb' for t in ticks]
#     #plt.xticks(ticks, [item[:-3] + 'kb' for item in ticks])
#     plt.show()
