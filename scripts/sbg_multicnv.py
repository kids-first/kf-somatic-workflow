import argparse
import json
import os
import logging
import numpy as np
import pandas as pd

import sbg_cveto_util as cveto
import sbg_multicnv_methods as multicnv
import sbg_multicnv_visualisations as plotly_vis

__author__ = 'Jelisaveta Ilic, jelisaveta.ilic@sevenbridges.com'
__created__ = '2021-02-08'
__copyright__ = 'SevenBridges 2021.'
__version__ = '1.0.0'
__desc__ = """Welcome to SevenBridges multiCNV tool! This tool combines 
multiCNV Combine, Segment and Merge tools in order to create the final results 
from multiple callers and tumors.."""


def parse_arguments():
    """
    Parse command-line arguments using argparse module.
    :return: Tuple of parsed arguments
    """
    parser = argparse.ArgumentParser(description=__desc__)
    parser.add_argument('-f', '--files',
                        nargs='+',
                        type=str,
                        help='Files as outputs of different callers / tumors',
                        required=True)
                        
    args = vars(parser.parse_args())

    return args['files']


def main():
    # collect user input
    files = parse_arguments()

    logger = cveto.init_logging('SBG MultiCNV')
    logger.info('START')
    logger.info(files)

    if len(files) == 0:
        logger.info("there is 0 files on input! please provide files")

    data_files = []
    ploidy_files = []
    dict_data = {}
    full_dict = {}
    full_dict_cn = {}
    full_dict_cr = {}
    callers_and_ploidy = {}
    callers = []
    tumors = []
    list_of_dicts = [full_dict, full_dict_cn, full_dict_cr]
    
    files_to_remove = []
    for filepath in files:
        filename = filepath.split('/')[-1]
        logger.info("Filename is: {}".format(filename))
        if filename.endswith('Store') or 'metrics' in filename:
            files_to_remove.append(filepath)
            continue
        if ('ploidy' in filename) and (filepath not in files_to_remove):
            ploidy_files.append(filepath)
        elif (filepath not in files_to_remove):
            data_files.append(filepath)
        else:
            logger.info("This file will not be used: {}".format(filepath))
    

    logger.info("Data files are: {}".format(data_files))
    logger.info("Ploidy files are: {}".format(ploidy_files))

    for data_file in data_files:

        filename = data_file.split('/')[-1]

        parts = filename.split('.')
        # Case id should be 41245 OR 22421
        patient = parts[0]

        # Sample should be rec1-lnl
        tumor = parts[1]
        # Caller should be freec or gatk
        caller = parts[2]
        #     print(caller)
        # Name should be 41245.rec1-lnl.2.gatk1
        name = '.'.join([patient, tumor, caller])
        callers.append(caller)
        tumors.append(tumor)

        dict_data[caller] = [data_file, None]

        for ploidy_file in ploidy_files:
            ploidy_filename = ploidy_file.split('/')[-1]
            # Check if it is its ploidy file
            if (caller in ploidy_filename) and (tumor in ploidy_filename):
                dict_data[caller] = [data_file, ploidy_file]

        df_per_file, regions_per_file, ploidy_per_file = multicnv.load_file(
            dict_data[caller], name)
        logging.info(f"Ploidy for {name} is {ploidy_per_file}.")
        callers_and_ploidy[name] = ploidy_per_file

        if name not in full_dict.keys():
            full_dict[name] = {}
        full_dict[name]['df'] = df_per_file
        full_dict[name]['df_per_file'] = df_per_file
        full_dict[name]['regions_per_file'] = regions_per_file
        full_dict[name]['regions'] = regions_per_file
        full_dict[name]['ploidy'] = ploidy_per_file

        # stavljam dict po fajlovima CN
        if name not in full_dict_cn.keys():
            full_dict_cn[name] = {}
        full_dict_cn[name]['df'] = df_per_file
        full_dict_cn[name]['df_per_file'] = df_per_file
        full_dict_cn[name]['regions_per_file'] = regions_per_file
        full_dict_cn[name]['regions'] = regions_per_file
        full_dict_cn[name]['ploidy'] = ploidy_per_file

        # CR
        if name not in full_dict_cr.keys():
            full_dict_cr[name] = {}
        full_dict_cr[name]['df'] = df_per_file
        full_dict_cr[name]['df_per_file'] = df_per_file
        full_dict_cr[name]['regions_per_file'] = regions_per_file
        full_dict_cr[name]['regions'] = regions_per_file
        full_dict_cn[name]['ploidy'] = ploidy_per_file

    names = []
    for name in full_dict:
        names.append(name)
        for i, name in enumerate(full_dict.keys()):
            if i == 0:
                regions = full_dict[name]['regions']
            else:
                regions = cveto.combine_regions(
                    regions, full_dict[name]['regions'], True
                )

        # STATUS
        df = cveto.to_dataframe(regions, full_dict.keys())
        df.fillna('neutral', inplace=True)
        full_dict[name]['df'] = df
        full_dict[name]['regions'] = regions
        groups = [col for col in df.columns if '_status' in col]
        stats = df.groupby(groups).agg({'length': 'sum'})
        full_dict[name]['stats'] = stats

        # CN dict
        df_cn = cveto.to_dataframe_with_cn(regions, full_dict.keys())
        df_cn.fillna(2,
                     inplace=True)  # ovde treba ploidy, videti kako da izvucem ploidy iz svakog caller-a
        full_dict_cn[name]['df'] = df_cn
        full_dict_cn[name]['regions'] = regions
        groups_cn = [col for col in df_cn.columns if '_status' in col]
        stats_cn = df_cn.groupby(groups_cn).agg({'length': 'sum'})
        full_dict_cn[name]['stats'] = stats_cn

        # CR dict

        df_cr = cveto.to_dataframe_with_cr(regions, callers_and_ploidy)
        df_cr.fillna(1,
                     inplace=True)  # ovde treba ploidy, videti kako da izvucem ploidy iz svakog caller-a
        full_dict_cr[name]['df'] = df_cr
        full_dict_cr[name]['regions'] = regions
        groups_cr = [col for col in df_cr.columns if '_status' in col]
        stats_cr = df_cr.groupby(groups_cr).agg({'length': 'sum'})
        full_dict_cr[name]['stats'] = stats_cr

    # Remove duplicate values in callers and tumors
    callers = list(dict.fromkeys(callers))
    tumors = list(dict.fromkeys(tumors))

    # Show full-dict
    full_dict['all'] = full_dict[name]['df']

    # All variables to lower case
    logger.info("Callers and ploidy dict is: {}".format(callers_and_ploidy))
    logger.info("Callers are: {}".format(callers))
    logger.info("Tumors are: {}".format(tumors))

    if len(callers) > len(tumors):
        title = tumors[0]
    else:
        title = callers[0]

    # first output
    with open('{}.ploidy_overview.json'.format(title), 'w') as ploidy_overview:
        json.dump(callers_and_ploidy, ploidy_overview)

    # full_dict['all'].to_csv('multiCNV_overview.csv')

    for name in names:
        full_dict[name]['df_per_file'].to_csv('{}.uniformed.csv'.format(name))

    logger.info("Done")

    # ALL FLOW - AUTOMATE MULTICALLER, MULTITUMOR
    mutual = 'multicnv'


    # Delete '1' regions
    df = multicnv.delete_small_and_na_regions(df)
    df_cn = multicnv.delete_small_and_na_regions(df_cn)
    df_cr = multicnv.delete_small_and_na_regions(df_cr)

    variables = []
    for column in df.columns:
        if 'status' in column:
            variables.append('_'.join(column.split('_')[:-1]))

    ######################### SAVE INITIAL DATAFRAMES STATUS, CN, CR #########################

    df.to_csv('combine_regions.status.csv')
    df_cn.to_csv('combine_regions.status_CN.csv')
    df_cr.to_csv('combine_regions.status_CR.csv')

    ######################### RENAME COLUMNS IN DATAFRAMES AND CALCULATE STATS #########################

    stats_dict = {}
    stats_total = multicnv.calculate_stats_total(df, '_status',
                                                 mutual,
                                                 'ALL')
    for variable in variables:
        stats_dict[variable] = multicnv.calculate_stats_per_file(
            stats_total, '{}_status'.format(variable))

    ######################### PLOTLY VISUALISATIONS #########################
    if len(data_files) == 3:
        plotly_vis.create_interactive_visualisations(df)

    ######################### CR HEATMAPS #########################

    df_del_columns = df_cr.copy()
    multicnv.delete_chr_start_end_columns(df_del_columns)

    # Calculate all values
    all_values = np.arange(0, 3, 0.1).tolist() + np.arange(3, 6,
                                                           0.2).tolist() + [6]
    all_values = np.around(all_values, decimals=1).tolist()

    variables_values_dict = {}
    for variable in variables:
        variables_values_dict[variable] = []

    df_uniform = df_del_columns.copy()

    # Adjusting values to nearest (1, 1.1 ...) to create heatmap
    for variable in variables:
        df_uniform, variables_values_dict[
            variable] = multicnv.adjust_values_to_nearest(
            df_del_columns, df_uniform, all_values,
            '{}_cr_status'.format(variable))

    # Set values which are in all_values and does not exist in df_uniform
    for value in all_values:
        for variable in variables:
            if value not in variables_values_dict[variable]:
                df_uniform = df_uniform.append(
                    {'{}_cr_status'.format(variable): value},
                    ignore_index=True)
    # Fill all NA values with 0
    df_uniform.fillna(0, inplace=True)

    # Generating all possible combinations of callers for CR df
    variables_to_combine_CR_list = multicnv.create_combination_of_columns_list(
        df_uniform)

    # Creating heatmaps for CR for all callers
    for variable_pairs in variables_to_combine_CR_list:
        multicnv.create_heatmap(variable_pairs[0], variable_pairs[1],
                                df_uniform,
                                mutual, 'CR')

    ######################### STATUS HEATMAPS #########################

    # Create heatmap for STATUS
    df_heatmap_status = df.copy()
    multicnv.delete_chr_start_end_columns(df_heatmap_status)
    df_calculate_metrics = df_heatmap_status.copy()

    df_heatmap_status = df_heatmap_status.replace(
        {'loss': 'A_loss', 'neutral': 'B_neutral', 'gain': 'C_gain'},
        regex=True)

    # Generating all possible combinations of callers for STATUS df
    variables_to_combine_STATUS_list = multicnv.create_combination_of_columns_list(
        df_heatmap_status)
    for variable_pairs in variables_to_combine_STATUS_list:
        multicnv.create_heatmap(variable_pairs[0], variable_pairs[1],
                                df_heatmap_status,
                                mutual, 'STATUS')

    ######################### INTERPRET - PRECISE #########################

    df_precise = df.copy()

    list_of_variables_index = multicnv.get_callers_index(df_precise)
    df_precise['final_call'] = df_precise.apply(
        lambda row: multicnv.precise_call(row, multicnv.get_callers_index(
            df_precise)), axis=1)
    multicnv.delete_status_columns(df_precise)
    df_precise.index = range(len(df_precise.index))
    df_precise.to_csv('{}.interpret_precise.csv'.format(mutual))

    ######################### SEGMENT - PRECISE #########################

    df_segment = pd.DataFrame(columns=df_precise.columns)
    df_segment_precise = multicnv.segment_regions(df_segment, df_precise,
                                                  mutual, 'PRECISE')
    multicnv.delete_ambiguous_regions(df_segment_precise, mutual,
                                      'PRECISE')
    stats_total_precise = multicnv.calculate_stats_total(df_segment_precise,
                                                         'final_call',
                                                         mutual, 'PRECISE')

    ######################### INTERPRET - MAJORITY #########################

    df_majority = df.copy()

    list_of_variables_index = multicnv.get_callers_index(df_majority)
    df_majority['final_call'] = df_majority.apply(
        lambda row: multicnv.majority_call(row, multicnv.get_callers_index(
            df_majority)), axis=1)
    multicnv.delete_status_columns(df_majority)
    df_majority.index = range(len(df_majority.index))
    df_majority.to_csv('{}.interpret_majority.csv'.format(mutual))

    ######################### SEGMENT - MAJORITY #########################
    df_segment = pd.DataFrame(columns=df_majority.columns)
    df_segment_majority = multicnv.segment_regions(df_segment, df_majority,
                                                   mutual,
                                                   'MAJORITY')
    multicnv.delete_ambiguous_regions(
        df_segment_majority,
        mutual,
        'MAJORITY'
    )
    stats_total_majority = multicnv.calculate_stats_total(df_segment_majority,
                                                          'final_call',
                                                          mutual, 'MAJORITY')
                                                          
    ######################### CALCULATE METRICS FOR TRUTH SET #################
    if 'truth' in callers:
        multicnv.calculate_metrics_for_benchmark(callers, df_calculate_metrics)


if __name__ == '__main__':
    main()
