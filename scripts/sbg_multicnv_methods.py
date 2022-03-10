import sbg_cveto_util as cveto
import pandas as pd
from collections import Counter
import itertools
import matplotlib.pyplot as plt
import seaborn as sb


# Read ploidy
def read_purple_ploidy(path):
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows,
                     header=header, comment=comment, sep=sep,
                     dtype={"ploidy": float, "purity": float})
    ploidy = df['ploidy'][0]
    ploidy_to_return = round(float(ploidy), 2)

    return ploidy_to_return


def read_controlfreec_ploidy(path):
    ploidy = 2
    with open(path) as openfile:
        for line in openfile:
            if 'Output_Ploidy' in line:
                parts = line.split('\t')
                ploidy = float(parts[1].split('\n')[0])

                return ploidy
    return ploidy


def read_sclust_ploidy(path):
    f = open(path, "r")
    header_line = f.readline()
    header_indexes = header_line.split('\t')
    value_line = f.readline()
    value_indexes = value_line.split('\t')
    data_dict = dict(zip(header_indexes, value_indexes))
    ploidy_value = (data_dict['ploidy'])
    ploidy = float(ploidy_value)
    ploidy = int(ploidy)
    return ploidy


def read_facets_ploidy(path):
    skiprows = 0
    header = 0
    sep = '\t'
    comment = '#'
    df = pd.read_csv(filepath_or_buffer=path, skiprows=skiprows,
                     header=header, comment=comment, sep=sep)
    ploidy = df['ploidy'][0]
    ploidy_to_return = round(float(ploidy), 2)

    return ploidy_to_return


# Additional methods

def load_file(cnv_file, name):
    caller = name
    df, ploidy = cveto.read_csv(cnv_file, caller, None)

    # Remove CHR, Chr, chr from chromosome name
    # This also ensures that X is X and not x
    df['chromosome'] = df['chromosome'].str.upper().str.replace('CHR', '')

    # Keep 1, 2, ... 22, X chromosomes
    valid_chroms = (
        '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y'.split()
    )
    cond = df['chromosome'].isin(valid_chroms)
    df = df[cond]

    return df, cveto.to_regions(df), ploidy


def delete_small_and_na_regions(df):
    df = df.drop(df[df.length == 1.0].index)
    df['length'] = pd.to_numeric(df.length, errors='coerce')
    df.dropna(inplace=True)
    return df


def calculate_stats_total(df, column_name_part, tumor, stats_type):
    groups = [col for col in df.columns if column_name_part in col]
    stats = df.groupby(groups).agg({'length': 'sum'})
    stats.to_csv('{}_{}.STATS.csv'.format(stats_type, tumor))
    return stats


def calculate_stats_per_file(stats, column):
    stats = stats.reset_index().groupby(column).sum()
    stats.to_csv('{}_.STATS_per_file.csv'.format(column))
    return stats


def delete_chr_start_end_columns(df):
    df.drop('chromosome', inplace=True, axis=1)
    df.drop('start', inplace=True, axis=1)
    df.drop('end', inplace=True, axis=1)


def delete_status_columns(df):
    for column in df.columns:
        if '_status' in column:
            df.drop(column, inplace=True, axis=1)


def adjust_values_to_nearest(df_old, df_new, all_values, column):
    caller_values = []
    for index, value in df_old[column].items():
        if value not in all_values:
            closest_value = min(all_values, key=lambda x: abs(x - value))
            df_new.at[index, column] = closest_value
            caller_values.append(closest_value)
        else:
            caller_values.append(value)
    return df_new, caller_values


def create_combination_of_columns_list(df):
    callers_to_combine = []
    for column in df.columns:
        if 'length' not in column:
            callers_to_combine.append(column)
    callers_to_combine_list = list(
        itertools.combinations(callers_to_combine, 2))
    return callers_to_combine_list


def create_heatmap(caller_one, caller_two, df, tumor, method):
    src = df.groupby([caller_one, caller_two]).agg({'length': 'sum'})
    src = src.unstack(level=0)
    src.columns = src.columns.get_level_values(1)
    src.fillna(0, inplace=True)
    fig, ax = plt.subplots(figsize=(11, 9))
    # plot heatmap
    if str(method) == 'STATUS':
        annot = True
        title = 'Comparison of statuses per length from two files'
    else:
        annot = False
        title = 'Comparison of cr values per length from two files'
    ay = sb.heatmap(src, cmap="Blues", linewidth=0.3,
                    cbar_kws={"shrink": .8, 'label': 'Length'}, annot=annot)
    ay.invert_yaxis()
    plt.xlabel(caller_one)
    plt.ylabel(caller_two)
    plt.title(title)
    plt.savefig(
        '{}_{}_{}_{}.heatmap.png'.format(caller_one, caller_two, method,
                                         tumor))
    plt.close(fig)
    return src


# Interpret


def get_callers_index(df):
    list_of_callers_index = []
    for i in range(0, len(df.columns)):
        if '_status' in df.columns[i]:
            list_of_callers_index.append(i)
    return list_of_callers_index


def precise_call(row, list_of_callers_index):
    """
    Return value is ambigous if not all callers agree on it.
    """
    all_calls = []
    for i in range(min(list_of_callers_index), max(list_of_callers_index) + 1):
        all_calls.append(row[i])
    c = Counter(all_calls)
    if len(c.most_common()) > 1:
        return 'ambiguous'
    else:
        value, count = c.most_common()[0]
        return value


def majority_call(row, list_of_callers_index):
    """
    Return value is the most common value. If two values are most common, result is ambiguous.
    """
    all_calls = []
    for i in range(min(list_of_callers_index), max(list_of_callers_index) + 1):
        all_calls.append(row[i])
    c = Counter(all_calls)
    if len(c.most_common()) > 1:
        if not (c.most_common()[0][1] == c.most_common()[1][1]):
            value, count = c.most_common()[0]
            return value
        else:
            return 'ambiguous'
    else:
        value, count = c.most_common()[0]
        return value


def delete_ambiguous_regions(df, tumor, interpret_method):
    df = df.drop(df[df.final_call == 'ambiguous'].index)
    if interpret_method == 'MAJORITY':
        add = 'majority_calls'
    elif interpret_method == 'PRECISE':
        add = 'precise_calls'
    else:
        add = 'undefined'
    df.to_csv('{}.{}.{}.FINAL_result.csv'.format(tumor, interpret_method, add))


def segment_regions(df_new, df_old, tumor, interpret_method):
    for row in df_old.iterrows():
        index_old = row[0]  # index of old df
        if df_new.empty:
            index_new = 0
            df_new.loc[df_old.index[index_new]] = df_old.iloc[index_old]
            continue
        if df_old['chromosome'].iloc[index_old] == df_new['chromosome'].iloc[
            index_new  # this will work because it will
            # always go through the first if
        ]:
            final_call_old = df_old['final_call'].iloc[index_old]
            final_call_new = df_new['final_call'].iloc[index_new]
            if final_call_old == final_call_new:
                df_new.at[df_new.index[index_new], 'end'] = (
                    df_old.iloc[index_old]['end']
                )
                df_new.at[df_new.index[index_new], 'length'] = (
                        df_old.iloc[index_old]['end'] - df_new.iloc[index_new][
                    'start']
                )
            else:
                df_new.loc[df_new.index.max() + 1] = df_old.iloc[index_old]
                index_new = index_new + 1
        else:
            df_new.loc[df_new.index.max() + 1] = df_old.iloc[index_old]
            index_new = index_new + 1
    df_new.to_csv(
        '{}_{}.seg_with_ambiguous.csv'.format(tumor, interpret_method))
    return df_new


# Plotly visualisations
def final_call(row, list_of_callers_index):
    all_calls = []
    for i in range(min(list_of_callers_index), max(list_of_callers_index) + 1):
        all_calls.append(row[i])
    all_calls.sort()
    return '-'.join(all_calls)
    
def calculate_metrics_for_benchmark(callers, df_calculate_metrics):
    """
    This should be run to calculate metrics when benchmarking multiple callers.
    """
    import json
    truth_combinations = []
    # Create combinations of all samples
    callers_to_combine_total = list(itertools.combinations(callers, 2))
    # Leave only combinations with truth set
    for callers_pair in callers_to_combine_total:
        for caller_name in callers_pair:
            if 'truth' in caller_name.lower():
                truth_combinations.append(list(callers_pair))

    metrics_dict = {}
    for truth_combination_pair in truth_combinations:
        # set truth set to be the first element in the list
        if 'truth' in truth_combination_pair[1].lower():
            truth_combination_pair[0], truth_combination_pair[1] = truth_combination_pair[1], truth_combination_pair[0]
        # Set correct column name
        for column in df_calculate_metrics.columns:
            if truth_combination_pair[0] in column:
                truth_column = ('_'.join(column.split('_')[:-1]))
            if truth_combination_pair[1] in column:
                caller_column = ('_'.join(column.split('_')[:-1]))
        # calculate metrics
        metrics_dict[truth_combination_pair[1]] = cveto.calculate_metrics(
            df_calculate_metrics, truth_column, caller_column
            )
        metrics_dict[truth_combination_pair[1]] = metrics_dict[truth_combination_pair[1]].to_dict()
    with open('benchmark_metrics.json', 'w') as metrics:
        json.dump(metrics_dict, metrics)
    df = pd.read_json('benchmark_metrics.json')
    df.to_csv('benchmark_metrics.csv')