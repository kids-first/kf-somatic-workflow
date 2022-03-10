import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.express as px

import sbg_multicnv_methods as multicnv
import sbg_cveto_util as cveto


def create_interactive_visualisations(df):
    logger = cveto.init_logging('SBG MultiCNV')
    dff = df.copy()

    # Create additional 'FINAL_status' column
    logger.info("Creating FINAL_status column")
    dff['FINAL_status'] = dff.apply(
        lambda row: multicnv.final_call(row, multicnv.get_callers_index(dff)),
        axis=1
    )

    # Extracting column names
    callers_columns = []
    callers_names = []
    pairs = {}

    for column_name in dff.columns:
        if '_status' in column_name:
            # Appending all columns with _status
            callers_columns.append(column_name)
            # Extracting variable part
            callers_name = column_name.split('_')[0]
            # callers_name = (column_name.split('.')[-1]).split('_')[0]
            # Appending variables / callers names
            callers_names.append(callers_name)
    logger.info(f"Callers / Tumors are: {callers_names[:-1]}")
    # print(callers_columns[:-1])
    # print(callers_names)

    callers_pair = {}
    for i in range(0, len(callers_names)):
        column_name = 'caller/tumor_' + str(i + 1)
        dff[column_name] = callers_names[i]
        callers_pair[column_name] = callers_names[i]  # example: caller_1:FREEC

    for column in callers_columns:
        for key in callers_pair.keys():
            if str(callers_pair[key]) in column:
                pairs[column] = key  # FREEC_status: caller_1

    df_vis = dff.copy()

    # Colors
    color_gain = '#FFFF00'  # yellow
    color_loss = '#F08080'  # red
    color_neutral = '#ADD8E6'  # blue

    color_gain_gain_gain = '#FFFF00'  # yellow
    color_gain_gain_loss = '#ff7f0e'  # safety orange
    color_gain_gain_neutral = '#FFFF00'  # yellow
    color_gain_loss_loss = '#ff7f0e'  # safety orange
    color_gain_neutral_neutral = '#FFFF00'  # yellow
    color_gain_loss_neutral = '#FF00FF'  # MAGENTA

    color_loss_loss_loss = '#F08080'  # red
    color_loss_loss_neutral = '#F08080'  # red
    color_loss_neutral_neutral = '#F08080'  # red

    color_neutral_neutral_neutral = '#ADD8E6'  # blue

    # Values
    value_gain = 1
    value_loss = -1
    value_neutral = 0

    value_gain_gain_gain = 3
    value_gain_gain_loss = 1
    value_gain_gain_neutral = 2
    value_gain_loss_loss = -1
    value_gain_neutral_neutral = 1
    value_gain_loss_neutral = 0

    value_loss_loss_loss = -3
    value_loss_loss_neutral = -2
    value_loss_neutral_neutral = -1

    value_neutral_neutral_neutral = 0

    # Dict for coloring
    final_col_dict = {
        'gain-gain-gain': {
            "color": color_gain_gain_gain, "value": value_gain_gain_gain
        },
        'gain-gain-loss': {
            "color": color_gain_gain_loss, "value": value_gain_gain_loss
        },
        'gain-gain-neutral': {
            "color": color_gain_gain_neutral, "value": value_gain_gain_neutral
        },
        'gain-loss-loss': {
            "color": color_gain_loss_loss, "value": value_gain_loss_loss
        },
        'gain-neutral-neutral': {
            "color": color_gain_neutral_neutral,
            "value": value_gain_neutral_neutral
        },
        'gain-loss-neutral': {
            "color": color_gain_loss_neutral, "value": value_gain_loss_neutral
        },
        'loss-loss-loss': {
            "color": color_loss_loss_loss, "value": value_loss_loss_loss
        },
        'loss-loss-neutral': {
            "color": color_loss_loss_neutral, "value": value_loss_loss_neutral
        },
        'loss-neutral-neutral': {
            "color": color_loss_neutral_neutral,
            "value": value_loss_neutral_neutral
        },
        'neutral-neutral-neutral': {
            "color": color_neutral_neutral_neutral,
            "value": value_neutral_neutral_neutral
        }
    }

    # Separate chromosomes into different dataframes
    chromosomes = dff['chromosome'].unique().tolist()
    dict_of_dfs = {}
    for chromosome in chromosomes:
        dict_of_dfs[f'dff_{chromosome}'] = dff[dff['chromosome'] == chromosome]

    # Creating plotly for CNV statuses
    logger.info("Creating plotly visualisation for CNV calls combined")
    dict_of_figures = {}
    fig_num = 1
    for dff in dict_of_dfs.keys():
        df = dict_of_dfs[dff]
        fig = make_subplots(rows=4, cols=1, subplot_titles=callers_names)

        row = 1
        for caller_column in callers_columns[:-1]:
            for (start, end, value, status, chromosome) in zip(
                    df["start"],
                    df["end"],
                    df["length"],
                    df[caller_column],
                    df['chromosome']
            ):
                if status == "gain":
                    line = dict(color=color_gain, width=4)
                    value = value_gain
                elif status == "neutral":
                    line = dict(color=color_neutral, width=4)
                    value = value_neutral
                else:
                    line = dict(color=color_loss, width=4)
                    value = value_loss
                fig.add_trace(go.Scatter(x=[start, end], y=[value, value],
                                         mode='lines', name=status, line=line,
                                         connectgaps=True), row=row, col=1)
                # chr_dict = count_traces_acc_chromosome(chr_dict, chromosome)

            row = row + 1

        for (start, end, value, status, chromosome) in zip(
                df["start"],
                df["end"],
                df["length"],
                df["FINAL_status"],
                df["chromosome"]
        ):
            for key_status in final_col_dict.keys():
                if status == key_status:
                    line = dict(color=final_col_dict[key_status]["color"],
                                width=4)
                    value = final_col_dict[key_status]["value"]

            fig.add_trace(go.Scatter(
                x=[start, end], y=[value, value], mode='lines', line=line,
                name=status, connectgaps=True
            ), row=4, col=1)
            # chr_dict = count_traces_acc_chromosome(chr_dict, chromosome)

        fig.update_layout(
            showlegend=False,
            title=f"Copy number regions - Chromosome {chromosome}"
        )
        dict_of_figures[f'fig_{fig_num}'] = fig
        fig_num = fig_num + 1

    # Saving figures
    with open('CNV_calls_combined.html', 'a') as f:
        for figure in dict_of_figures.keys():
            fig = dict_of_figures[figure]
            f.write(fig.to_html(full_html=True, include_plotlyjs='cdn'))

    logger.info("Saving HTML: CNV calls combine")

    # Statistics - sum of lengths
    logger.info("Creating plotly visualisation for sum of lengths")
    color_discrete_map = {
        'gain': '#FFFF00', 'neutral': '#ADD8E6', 'loss': '#F08080',
        'gain-gain-gain': '#FFFF00',
        'gain-gain-loss': '#FFA500',
        'gain-gain-neutral': '#90EE90',
        'gain-loss-loss': '#FF4500',
        'gain-neutral-neutral': '#2E8B57',
        'gain-loss-neutral': '#FF00FF',
        'loss-loss-loss': '#F08080',
        'loss-loss-neutral': '#DDA0DD',
        'loss-neutral-neutral': '#800080',
        'neutral-neutral-neutral': '#ADD8E6'
    }
    df = df_vis.copy()
    list_of_figs = []
    for key in pairs.keys():
        fig = px.histogram(
            df, x="length", y=pairs[key], color=key,
            hover_name="length", title='Sum of lengths', histfunc='sum',
            color_discrete_map=color_discrete_map, facet_col='chromosome',
            facet_col_wrap=5
        )
        list_of_figs.append(fig)

    # Saving figures
    with open('sum_of_lengths_figures.html', 'a') as f:
        for fig in list_of_figs:
            f.write(fig.to_html(full_html=True, include_plotlyjs='cdn'))
    logger.info("Saving HTML: SUM of lengths per chromosome")

    # Statistics - Total sum of lengths
    logger.info("Creating plotly visualisation for Total sum of lengths")
    df = df_vis.copy()
    list_of_figs = []
    for key in pairs.keys():
        fig = px.histogram(
            df, x="length", y=pairs[key], color=key,
            hover_name="length", title='Sum of lengths', histfunc='sum',
            color_discrete_map=color_discrete_map
        )
        list_of_figs.append(fig)

    # Saving figures
    with open('TOTAL_sum_of_lengths_figures.html', 'a') as f:
        for fig in list_of_figs:
            f.write(fig.to_html(full_html=True, include_plotlyjs='cdn'))
    logger.info("Saving HTML: Total SUM of lengths")
