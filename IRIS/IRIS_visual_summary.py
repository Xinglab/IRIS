from __future__ import print_function

import collections
import os
import sys

import matplotlib
matplotlib.use('agg')  # sets the plotting mode to non-interactive

import matplotlib.gridspec as gridspec
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

COLOR_BLACK = '#000000'
COLOR_BURNT_ORANGE = '#c04e01'
COLOR_CREAM = '#ffffc2'
COLOR_GREEN = '#15b01a'
COLOR_LIGHT_VIOLET = '#d6b4fc'
COLOR_OCHRE = '#bf9005'

COLOR_BY_PANEL_TYPE = {
    'output': COLOR_BLACK,
    'tissue_matched_normal': COLOR_GREEN,
    'tumor': COLOR_BURNT_ORANGE,
    'normal': COLOR_GREEN
}

Z_LOWEST = 0

# size in points
FONT_SIZE_MEDIUM = 11
FONT_SIZE_LARGE = 13


def points_to_pixels(points):
    # 1 point = (1/72.0) inches
    # 1 inch = 96.0 pixels
    inches = points / 72.0
    return inches * 96.0


def get_matrix_index(event_id, index_f_name):
    with open(index_f_name, 'rt') as f_h:
        for line in f_h:
            event, offset_str = line.strip().split('\t')
            if event == event_id:
                offset_int, error = parse_int(offset_str)
                if error:
                    return None, error

                return offset_int, None

    return None, '{} not in {}'.format(event_id, index_f_name)


def read_event_row(event_id, matrix_f_name, index_f_name):
    row = collections.OrderedDict()
    offset, error = get_matrix_index(event_id, index_f_name)
    fill_nan = False
    if error:
        print('{}; filling with NaN values'.format(error), file=sys.stderr)
        fill_nan = True

    with open(matrix_f_name, 'rt') as f_h:
        raw_headers = f_h.readline().strip().split('\t')
        header = np.asarray([
            x.split('.')[0] if x.startswith('SRR') else x for x in raw_headers
        ])

        if fill_nan:
            data = np.asarray(['NaN' for _ in header])
        else:
            f_h.seek(offset, 0)
            data = np.asarray(f_h.readline().strip().split('\t'))

        row = collections.OrderedDict(zip(header, data))

    return row, None


def read_panel_parameters(line):
    parameters = dict()
    if not line.strip():
        return parameters, None

    splits = line.split(' ')
    if len(splits) != 2:
        return parameters, 'expected 2 fields in {} but found {}'.format(
            line, len(splits))

    cutoffs_str, groups_str = splits
    cutoffs = cutoffs_str.split(',')
    groups = groups_str.split(',')
    floats = list()
    for cutoff in cutoffs:
        parsed_float, error = parse_float(cutoff)
        if error:
            return parameters, error

        floats.append(parsed_float)

    parameters['psi_pval_cutoff'] = floats[0]
    parameters['delta_psi_cutoff'] = floats[1]
    parameters['foc_cutoff'] = floats[2]
    parameters['sjc_pval_cutoff'] = floats[3]
    parameters['group_cutoff'] = floats[4]
    parameters['groups'] = groups
    return parameters, None


def read_parameters(parameter_f_name):
    parameters = dict()

    with open(parameter_f_name, 'rt') as f_h:
        lines = [line.strip() for line in f_h]

    if len(lines) != 10:
        return parameters, 'expected 10 lines in {} but found {}'.format(
            parameter_f_name, len(lines))

    parameters['out_prefix'] = lines[0]
    parameters['db_dir'] = lines[1]
    panel_params = list()
    for line in lines[2:5]:
        this_panel_params, error = read_panel_parameters(line)
        if error:
            return parameters, error

        panel_params.append(this_panel_params)

    parameters['tissue_matched_normal_panel'] = panel_params[0]
    parameters['tumor_panel'] = panel_params[1]
    parameters['normal_panel'] = panel_params[2]
    parameters['test_mode'] = lines[5]
    parameters['use_ratio'] = lines[6] == 'True'
    parameters['blacklist'] = lines[7]
    parameters['mappability_path'] = lines[8]
    parameters['ref_genome_path'] = lines[9]

    return parameters, None


def get_group_to_panel_type(parameters):
    panel_tmn = parameters['tissue_matched_normal_panel'].get('groups', list())
    panel_t = parameters['tumor_panel'].get('groups', list())
    panel_n = parameters['normal_panel'].get('groups', list())
    output_group = parameters['out_prefix']

    group_to_panel_type = collections.OrderedDict()
    group_to_panel_type[output_group] = 'output'

    for group in panel_tmn:
        group_to_panel_type[group] = 'tissue_matched_normal'

    for group in panel_t:
        group_to_panel_type[group] = 'tumor'

    for group in panel_n:
        group_to_panel_type[group] = 'normal'

    return group_to_panel_type


def get_matrix_file_names(groups, parameters, event_type):
    group_to_matrix_f_name = collections.OrderedDict()
    group_to_matrix_index_f_name = collections.OrderedDict()

    if not groups:
        return group_to_matrix_f_name, group_to_matrix_index_f_name, 'no groups specified'

    db_dir = parameters['db_dir']
    for group_name in groups:
        splicing_f_name = 'splicing_matrix.{}.cov10.{}.txt'.format(
            event_type, group_name)
        group_to_matrix_f_name[group_name] = os.path.join(
            db_dir, group_name, 'splicing_matrix', splicing_f_name)

    for group, matrix_f_name in group_to_matrix_f_name.items():
        if not os.path.isfile(matrix_f_name):
            error = 'no matrix file found for {}. expected {}'.format(
                group, matrix_f_name)
            return group_to_matrix_f_name, group_to_matrix_index_f_name, error

        index_f_name = '{}.idx'.format(matrix_f_name)
        if not os.path.isfile(index_f_name):
            error = 'no index file found for {}. expected {}'.format(
                group, matrix_f_name)
            return group_to_matrix_f_name, group_to_matrix_index_f_name, error

        group_to_matrix_index_f_name[group] = index_f_name

    return group_to_matrix_f_name, group_to_matrix_index_f_name, None


def get_events(parameters, screening_out_dir, event_type):
    events = dict()

    for variant in ['tier1', 'tier2tier3']:
        variant_events = list()
        events[variant] = variant_events
        file_name = '{}.{}.{}.txt'.format(parameters['out_prefix'], event_type,
                                          variant)
        file_path = os.path.join(screening_out_dir, file_name)
        if not os.path.isfile(file_path):
            return events, 'missing required file: {}'.format(file_path)

        event_id_index = None
        with open(file_path, 'rt') as f_h:
            for i, line in enumerate(f_h):
                splits = line.strip().split('\t')
                if i == 0:  # header
                    event_id_header = 'as_event'
                    if event_id_header not in splits:
                        return events, 'required header {} not found in {}'.format(
                            event_id_header, file_path)

                    event_id_index = splits.index(event_id_header)
                    continue

                if len(splits) <= event_id_index:
                    return events, 'no {} for line {} of {}'.format(
                        event_id_header, i, file_path)

                event_id = splits[event_id_index]
                variant_events.append(event_id)

    return events, None


def get_psi_data_by_event(group_to_matrix_f_name, group_to_matrix_index_f_name,
                          events):
    psi_data_by_event = collections.OrderedDict()
    tier2tier3_events = events['tier2tier3']
    if not tier2tier3_events:
        return psi_data_by_event, 'no tier2tier3 events'

    groups = group_to_matrix_f_name.keys()
    for event_id in tier2tier3_events:
        psi_by_group = collections.OrderedDict()
        all_psi = list()
        for group in groups:
            row, error = read_event_row(event_id,
                                        group_to_matrix_f_name[group],
                                        group_to_matrix_index_f_name[group])
            if error:
                return psi_data_by_event, error

            psi_strings = list(row.values())[8:]
            psi_floats = list()
            for psi_s in psi_strings:
                psi_f, error = parse_float(psi_s)
                if error:
                    return psi_data_by_event, error

                psi_floats.append(psi_f)

            psi_by_group[group] = psi_floats
            all_psi.extend(psi_floats)

        abs_change = max(all_psi) - min(all_psi)
        if abs_change < 0.05:
            continue

        psi_df = pd.DataFrame.from_dict(psi_by_group,
                                        orient='index').transpose()[groups]
        psi_data_by_event[event_id] = psi_df

    return psi_data_by_event, None


def add_or_verify_match(key, source_dict, dest_dict, parser):
    if key not in source_dict:
        return 'missing {}'.format(key)

    source_v = source_dict[key]
    parsed_source_v, error = parser(source_v)
    if error:
        return error

    if key in dest_dict:
        existing_v = dest_dict[key]
        if parsed_source_v == existing_v:
            return None

        return 'differing values for {}: {}, {}'.format(
            key, parsed_source_v, existing_v)

    dest_dict[key] = parsed_source_v
    return None


def read_tsv(f_name):
    rows = list()
    header = list()
    with open(f_name, 'rt') as f_h:
        for i, line in enumerate(f_h):
            tokens = line.strip().split('\t')
            if i == 0:
                header = tokens
                continue

            if len(tokens) != len(header):
                return header, rows, 'expected {} columns at line {} of {} but found {}'.format(
                    len(header), i, f_name, len(tokens))

            rows.append(dict(zip(header, tokens)))

    return header, rows, None


def get_hlas_with_binding_affinity(hla_headers, row, row_num, f_name):
    hlas_with_binding_affinity = list()
    for header in hla_headers:
        full_str = row[header]
        if full_str == '-':
            continue

        col_error_prefix = 'row {} in {} with column {}={}: '.format(
            row_num, f_name, header, full_str)
        semi_splits = full_str.split(';')
        for semi_split in semi_splits:
            pipe_splits = semi_split.split('|')
            if len(pipe_splits) != 2:
                return hlas_with_binding_affinity, '{}expected exactly 1 "|" in {}'.format(
                    col_error_prefix, semi_split)

            hla_str, binding_str = pipe_splits
            hla_prefix = 'HLA-'
            if not hla_str.startswith(hla_prefix):
                return hlas_with_binding_affinity, '{}expected a value starting with "{}"'.format(
                    col_error_prefix, hla_prefix)

            hla_sub_str = hla_str[len(hla_prefix):]
            binding_float, error = parse_float(binding_str)
            if error:
                return hlas_with_binding_affinity, '{}{}'.format(
                    col_error_prefix, error)

            hlas_with_binding_affinity.append((hla_sub_str, binding_float))

    return hlas_with_binding_affinity, None


def parse_int_ratio(s):
    splits = s.split('/')
    if len(splits) != 2:
        return None, 'expected exactly one "/" in {}'.format(s)

    ints = list()
    for split in splits:
        as_int, error = parse_int(split)
        if error:
            return None, error

        ints.append(as_int)

    return ints, None


def parse_float(s):
    try:
        as_float = float(s)
    except ValueError as e:
        return None, 'could not parse {} as float: {}'.format(s, e)

    return as_float, None


def parse_int(s):
    try:
        as_int = int(s)
    except ValueError as e:
        return None, 'could not parse {} as int: {}'.format(s, e)

    return as_int, None


def process_epitope_summary_row(row, row_num, summary, hla_headers, f_name,
                                has_prediction):
    event_id = row['as_event']
    event_summary = summary.get(event_id)
    if not event_summary:
        event_summary = dict()
        summary[event_id] = event_summary

    event_epitopes_hla_affinity_patients = event_summary.get('epitopes')
    if not event_epitopes_hla_affinity_patients:
        event_epitopes_hla_affinity_patients = list()
        event_summary['epitopes'] = event_epitopes_hla_affinity_patients

    row_error_prefix = 'row {} in {}: '.format(row_num, f_name)
    int_ratio_keys = [
        'tissue_matched_normal_panel', 'tumor_panel', 'normal_panel'
    ]
    for key in int_ratio_keys:
        error = add_or_verify_match(key, row, event_summary, parse_int_ratio)
        if error:
            return '{}{}'.format(row_error_prefix, error)

    error = add_or_verify_match('fc_of_tumor_isoform', row, event_summary,
                                parse_float)
    if error:
        return '{}{}'.format(row_error_prefix, error)

    if not has_prediction:
        return None

    hlas_with_binding_affinity, error = get_hlas_with_binding_affinity(
        hla_headers, row, row_num, f_name)
    if error:
        return error

    if not hlas_with_binding_affinity:
        return '{}no HLAs found'.format(row_error_prefix)

    hlas_with_binding_affinity.sort(key=lambda p: p[1])
    hla, binding_affinity = hlas_with_binding_affinity[0]
    event_epitopes_hla_affinity_patients.append({
        'epitope':
        row['epitope'],
        'hla':
        hla,
        'binding_affinity':
        binding_affinity,
        'num_patients':
        row['num_sample']
    })

    return None


def get_epitope_summary(screening_out_dir, has_prediction, parameters,
                        event_type):
    summary = dict()

    out_prefix = parameters['out_prefix']
    if has_prediction:
        summary_f_name = os.path.join(screening_out_dir,
                                      '{}.tier2tier3'.format(event_type),
                                      'epitope_summary.peptide-based.txt')
    else:
        summary_f_name = os.path.join(
            screening_out_dir,
            '{}.{}.tier2tier3.txt'.format(out_prefix, event_type))

    header, rows, error = read_tsv(summary_f_name)
    if error:
        return summary, error

    expected_columns = {
        'as_event', 'meanPSI', 'Q1PSI', 'Q3PSI', 'deltaPSI',
        'fc_of_tumor_isoform', 'tissue_matched_normal_panel', 'tumor_panel',
        'normal_panel', 'tag', 'mappability', 'mappability_tag'
    }
    if has_prediction:
        expected_columns = expected_columns.union({
            'epitope', 'junction_peptide_form', 'inclusion_form', 'num_hla',
            'num_sample', 'hla_types', 'canonical_match', 'uniqueness'
        })

    ignored_optional_columns = {'meanGeneExp', 'Q1GeneExp', 'Q3GeneExp'}
    header_set = set(header)
    missing_headers = expected_columns.difference(header_set)
    if missing_headers:
        return summary, 'missing headers in {}: {}'.format(
            summary_f_name, ', '.join(sorted(missing_headers)))

    ignored_columns = expected_columns.union(ignored_optional_columns)
    extra_headers = header_set.difference(ignored_columns)
    for i, row in enumerate(rows):
        error = process_epitope_summary_row(row, i, summary, extra_headers,
                                            summary_f_name, has_prediction)
        if error:
            return summary, error

    if has_prediction:
        for event_summary in summary.values():
            event_summary['epitopes'].sort(key=lambda p: p['binding_affinity'])

    return summary, None


def remove_spines(ax):
    for spine in ax.spines.values():
        spine.set_visible(False)


def remove_ticks_and_labels(ax):
    ax.set_xticks(list())
    ax.set_yticks(list())


def remove_spines_ticks_and_labels(ax):
    remove_spines(ax)
    remove_ticks_and_labels(ax)


def hide_ax(ax):
    ax.set_visible(False)


def hide_ax_except_color(ax, color):
    ax.clear()
    ax.set_zorder(Z_LOWEST)
    remove_spines_ticks_and_labels(ax)
    ax.set_facecolor(color)


def make_violin_plots(events, psi_data_by_event, groups, group_to_panel_type,
                      axes_by_name):
    gene_header_ax, gene_axes = axes_by_name['gene_name']
    violins_header_ax, violins_axes = axes_by_name['violins']
    violins_y_ticks_header_ax, violins_y_ticks_axes = axes_by_name[
        'violins_y_ticks']

    hide_ax_except_color(violins_y_ticks_header_ax, COLOR_CREAM)
    for y_tick_ax in violins_y_ticks_axes:
        hide_ax(y_tick_ax)

    make_gene_header(gene_header_ax)

    if len(events) != len(psi_data_by_event) or len(events) != len(
            gene_axes) or len(events) != len(violins_axes):
        return 'all inputs must have the same length'

    sns.set(style="white", color_codes=True)
    for i, event in enumerate(events):
        event = events[i]
        psi_data = psi_data_by_event.get(event)
        if psi_data is None:
            return 'no psi data for {}'.format(event)

        gene_ax = gene_axes[i]
        violins_ax = violins_axes[i]

        sns.violinplot(data=psi_data,
                       ax=violins_ax,
                       inner="box",
                       cut=0,
                       scale='width',
                       linewidth=1.5)

        violins_ax.set_yticks(np.arange(0, 1.1, 0.5))
        violins_ax.set_xticklabels(list())
        violins_ax.xaxis.set_tick_params(which='both', length=0)
        sns.despine(ax=violins_ax, offset=0, trim=False)

        remove_spines_ticks_and_labels(gene_ax)
        gene_name = event.split(':')[1]
        gene_ax.text(0.25,
                     0.5,
                     gene_name,
                     horizontalalignment='center',
                     verticalalignment='center',
                     transform=gene_ax.transAxes,
                     style='italic',
                     size=FONT_SIZE_LARGE,
                     color=COLOR_BLACK,
                     fontweight='bold')

    top_violin_ax = violins_axes[0]
    make_violins_header(violins_header_ax, top_violin_ax, groups,
                        group_to_panel_type)
    return None


def make_gene_header(gene_header_ax):
    gene_header_ax.text(0.25,
                        0,
                        'Gene of\nAS event',
                        horizontalalignment='center',
                        verticalalignment='bottom',
                        transform=gene_header_ax.transAxes,
                        size=FONT_SIZE_LARGE,
                        color=COLOR_BLACK,
                        fontweight='bold')

    remove_spines_ticks_and_labels(gene_header_ax)


def make_violins_header(violin_header_ax, top_violin_ax, groups,
                        group_to_panel_type):
    violin_xticks_data_coords = [(x, 0) for x in top_violin_ax.get_xticks()]
    violin_xticks_display_coords = top_violin_ax.transData.transform(
        violin_xticks_data_coords)
    violin_tick_label_y_display_coord = violin_header_ax.transAxes.transform_point(
        (0, 0))[1]
    violin_tick_label_display_coords = [(p[0],
                                         violin_tick_label_y_display_coord)
                                        for p in violin_xticks_display_coords]
    trans_display_to_axes = violin_header_ax.transAxes.inverted(
    ).transform_point
    for i, coord in enumerate(violin_tick_label_display_coords):
        x, y = trans_display_to_axes(coord)
        group = groups[i]
        color = COLOR_BY_PANEL_TYPE[group_to_panel_type[group]]
        violin_header_ax.text(x,
                              y,
                              group,
                              horizontalalignment='center',
                              verticalalignment='bottom',
                              rotation=90,
                              transform=violin_header_ax.transAxes,
                              size=FONT_SIZE_MEDIUM,
                              color=color,
                              fontweight='bold')

    violin_header_ax.text(0.5,
                          1,
                          'PSI by tissue or tumor',
                          horizontalalignment='center',
                          verticalalignment='top',
                          transform=violin_header_ax.transAxes,
                          size=FONT_SIZE_LARGE,
                          color=COLOR_BLACK,
                          fontweight='bold')
    text_height = points_to_pixels(FONT_SIZE_LARGE)
    top_text_y_max_display = violin_header_ax.transAxes.transform_point(
        (0, 1))[1]
    top_text_y_min_display = top_text_y_max_display - text_height
    top_text_y_min_axes = trans_display_to_axes((0, top_text_y_min_display))[1]
    underline_x_min_display = violin_xticks_display_coords[0][0]
    underline_x_max_display = violin_xticks_display_coords[-1][0]
    underline_x_min_axes = trans_display_to_axes(
        (underline_x_min_display, 0))[0]
    underline_x_max_axes = trans_display_to_axes(
        (underline_x_max_display, 0))[0]
    violin_header_ax.plot([underline_x_min_axes, underline_x_max_axes],
                          [top_text_y_min_axes] * 2,
                          color=COLOR_BLACK,
                          linestyle='solid',
                          transform=violin_header_ax.transAxes)

    remove_spines_ticks_and_labels(violin_header_ax)
    violin_header_ax.set_facecolor(COLOR_CREAM)


def create_grid_of_axes(fig, num_events, has_prediction):
    content_rows = num_events + 1
    spacer_rows = content_rows - 1
    content_rows_to_spacer = 9
    extra_header_room = content_rows_to_spacer
    grid_rows = (content_rows_to_spacer * content_rows) + spacer_rows
    grid_rows += extra_header_room

    grid_col_names_and_widths = [('gene_name', 3), ('violins_y_ticks', 2),
                                 ('violins', 13), ('tissue_matched', 1),
                                 ('fold_change', 2), ('tumor', 1),
                                 ('normal', 1)]
    if has_prediction:
        grid_col_names_and_widths.extend([('epitopes', 5), ('hlas', 3),
                                          ('num_patients', 2),
                                          ('affinity', 3)])

    grid_col_intervals_by_name = dict()
    grid_cols = 0
    for name, width in grid_col_names_and_widths:
        start_cols = grid_cols
        grid_cols += width
        grid_col_intervals_by_name[name] = (start_cols, grid_cols)

    grid = gridspec.GridSpec(grid_rows, grid_cols, wspace=0, hspace=0)

    filler_rows = list()
    axes_by_name = dict()
    for name, col_interval in grid_col_intervals_by_name.items():
        start_col, end_col = col_interval
        axes = list()

        row = 0
        for _ in range(0, content_rows):
            is_header = row == 0
            if not is_header:
                filler_rows.append(
                    fig.add_subplot(grid[row, start_col:end_col]))
                row += 1

            start_row = row
            row += content_rows_to_spacer
            if is_header:
                row += extra_header_room

            axes.append(fig.add_subplot(
                grid[start_row:row, start_col:end_col]))

        axes_by_name[name] = (axes[0], axes[1:])

    return axes_by_name, filler_rows


def make_shaded_dots_column(header_ax, dots_axes, header_text, events,
                            epitope_summary, summary_ratio_key, color):
    header_ax.text(0.5,
                   0,
                   'vs. ',
                   horizontalalignment='center',
                   verticalalignment='bottom',
                   rotation=90,
                   transform=header_ax.transAxes,
                   size=FONT_SIZE_MEDIUM,
                   color=COLOR_BLACK,
                   fontweight='bold')

    y_min_display = header_ax.transAxes.transform_point((0, 0))[1]
    pixels_per_large_char = points_to_pixels(FONT_SIZE_MEDIUM)
    # 'vs. ' is about two large characters in the default variable width font
    pixels_of_text = pixels_per_large_char * 2
    vs_space_y_max_display = y_min_display + pixels_of_text

    trans_display_to_axes = header_ax.transAxes.inverted().transform_point
    vs_space_y_max_axes = trans_display_to_axes((0, vs_space_y_max_display))[1]
    # put the header_text after 'vs. ' in a (possibly) different color
    header_ax.text(0.5,
                   vs_space_y_max_axes,
                   header_text,
                   horizontalalignment='center',
                   verticalalignment='bottom',
                   rotation=90,
                   transform=header_ax.transAxes,
                   size=FONT_SIZE_MEDIUM,
                   color=color,
                   fontweight='bold')

    remove_spines_ticks_and_labels(header_ax)
    header_ax.set_facecolor(COLOR_CREAM)

    for i, event in enumerate(events):
        event_summary = epitope_summary.get(event)
        if not event_summary:
            return 'no event {} in epitope_summary'.format(event)

        hits, total = event_summary[summary_ratio_key]
        if total == 0:
            percent = 0
        else:
            percent = hits / float(total)

        dots_ax = dots_axes[i]
        xmin_display, ymin_display = dots_ax.transAxes.transform_point((0, 0))
        xmax_display, ymax_display = dots_ax.transAxes.transform_point((1, 1))
        x_spread_display = xmax_display - xmin_display
        y_spread_display = ymax_display - ymin_display
        # dots_ax is a rectangle.
        # Normalize the axes scale so a circle is not distorted
        x_mid_data = 0.5
        y_mid_data = 0.5
        if x_spread_display > y_spread_display:
            new_scale = x_spread_display / float(y_spread_display)
            half_new_scale = new_scale / 2.0
            dots_ax.set_xlim((0, new_scale))
            x_mid_data = half_new_scale
        else:
            new_scale = y_spread_display / float(x_spread_display)
            half_new_scale = new_scale / 2.0
            dots_ax.set_ylim((0, new_scale))
            y_mid_data = half_new_scale

        circle = matplotlib.patches.Circle((x_mid_data, y_mid_data),
                                           radius=0.25,
                                           alpha=percent,
                                           color=COLOR_OCHRE,
                                           transform=dots_ax.transData)
        dots_ax.add_patch(circle)

        remove_spines_ticks_and_labels(dots_ax)

    return None


def show_fold_change_values(header_ax, fc_axes, events, epitope_summary):
    header_ax.text(0.5,
                   0,
                   'FC of tumor\nisoform',
                   horizontalalignment='center',
                   verticalalignment='bottom',
                   rotation=90,
                   transform=header_ax.transAxes,
                   size=FONT_SIZE_MEDIUM,
                   color=COLOR_BLACK,
                   fontweight='bold')
    remove_spines_ticks_and_labels(header_ax)
    header_ax.set_facecolor(COLOR_CREAM)

    for i, event in enumerate(events):
        event_summary = epitope_summary.get(event)
        if not event_summary:
            return 'no event {} in epitope_summary'.format(event)

        fc = event_summary['fc_of_tumor_isoform']
        fc_ax = fc_axes[i]
        fc_ax.text(0.5,
                   0.5,
                   '{:.1f}'.format(fc),
                   horizontalalignment='center',
                   verticalalignment='center',
                   transform=fc_ax.transAxes,
                   size=FONT_SIZE_MEDIUM)
        remove_spines_ticks_and_labels(fc_ax)

    return None


def make_top_underline_header(left_header_ax, right_header_ax, other_axes,
                              header_text):
    # need to plot on all axes using display coordinates to avoid clipping
    all_axes = [left_header_ax, right_header_ax] + other_axes

    x_min_display = left_header_ax.transAxes.transform_point((0, 0))[0]
    x_max_display = right_header_ax.transAxes.transform_point((1, 0))[0]
    x_mid_display = (x_min_display + x_max_display) / 2.0
    y_max_display = left_header_ax.transAxes.transform_point((0, 1))[1]
    x_range_display = x_max_display - x_min_display
    underline_margin_display = 0.02 * x_range_display
    underline_x_min_display = x_min_display + underline_margin_display
    underline_x_max_display = x_max_display - underline_margin_display
    for ax in all_axes:
        trans_display_to_axes = ax.transAxes.inverted().transform_point
        x_mid_axes, y_max_axes = trans_display_to_axes(
            (x_mid_display, y_max_display))
        ax.text(x_mid_axes,
                y_max_axes,
                header_text,
                horizontalalignment='center',
                verticalalignment='top',
                transform=ax.transAxes,
                size=FONT_SIZE_LARGE,
                color=COLOR_BLACK,
                fontweight='bold',
                clip_on=True)

    text_y_min_display = y_max_display - points_to_pixels(FONT_SIZE_LARGE)

    for ax in all_axes:
        trans_display_to_axes = ax.transAxes.inverted().transform_point
        underline_x_min_axes, text_y_min_axes = trans_display_to_axes(
            (underline_x_min_display, text_y_min_display))
        underline_x_max_axes = trans_display_to_axes(
            (underline_x_max_display, 0))[0]
        ax.plot([underline_x_min_axes, underline_x_max_axes],
                [text_y_min_axes] * 2,
                color=COLOR_BLACK,
                linestyle='solid',
                transform=ax.transAxes,
                clip_on=True)


def make_shaded_dots_and_show_fold_change(events, epitope_summary,
                                          axes_by_name):
    tmn_header_ax, tmn_axes = axes_by_name['tissue_matched']
    t_header_ax, t_axes = axes_by_name['tumor']
    n_header_ax, n_axes = axes_by_name['normal']
    fc_header_ax, fc_axes = axes_by_name['fold_change']

    error = make_shaded_dots_column(
        tmn_header_ax, tmn_axes, 'Tissue Matched', events, epitope_summary,
        'tissue_matched_normal_panel',
        COLOR_BY_PANEL_TYPE['tissue_matched_normal'])
    if error:
        return error

    error = make_shaded_dots_column(t_header_ax, t_axes, 'Tumor', events,
                                    epitope_summary, 'tumor_panel',
                                    COLOR_BY_PANEL_TYPE['tumor'])
    if error:
        return error

    error = make_shaded_dots_column(n_header_ax, n_axes, 'Normal', events,
                                    epitope_summary, 'normal_panel',
                                    COLOR_BY_PANEL_TYPE['normal'])
    if error:
        return error

    show_fold_change_values(fc_header_ax, fc_axes, events, epitope_summary)

    make_top_underline_header(tmn_header_ax, n_header_ax,
                              [t_header_ax, fc_header_ax], 'Summary')

    return None


def show_epitope_names(header_ax, epitope_axes, events, epitope_summary):
    header_ax.text(0.5,
                   0,
                   'Junction\nepitopes',
                   horizontalalignment='center',
                   verticalalignment='bottom',
                   transform=header_ax.transAxes,
                   size=FONT_SIZE_LARGE,
                   color=COLOR_BLACK,
                   fontweight='bold')
    remove_spines_ticks_and_labels(header_ax)
    header_ax.set_facecolor(COLOR_LIGHT_VIOLET)

    for i, event in enumerate(events):
        event_summary = epitope_summary.get(event)
        if not event_summary:
            return 'no event {} in epitope_summary'.format(event)

        epitopes = event_summary['epitopes']
        epitope_names = [e['epitope'] for e in epitopes[:2]]
        name_text = '\n'.join(epitope_names)
        epitope_ax = epitope_axes[i]
        epitope_ax.text(0.5,
                        0.5,
                        name_text,
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=epitope_ax.transAxes,
                        size=FONT_SIZE_MEDIUM)
        remove_spines_ticks_and_labels(epitope_ax)

    return None


def show_hla_names(header_ax, hla_axes, events, epitope_summary):
    header_ax.text(0.5,
                   0,
                   'Best\nHLA',
                   horizontalalignment='center',
                   verticalalignment='bottom',
                   transform=header_ax.transAxes,
                   size=FONT_SIZE_LARGE,
                   color=COLOR_BLACK,
                   fontweight='bold')
    remove_spines_ticks_and_labels(header_ax)
    header_ax.set_facecolor(COLOR_LIGHT_VIOLET)

    for i, event in enumerate(events):
        event_summary = epitope_summary.get(event)
        if not event_summary:
            return 'no event {} in epitope_summary'.format(event)

        epitopes = event_summary['epitopes']
        epitope_hlas = [e['hla'] for e in epitopes[:2]]
        hla_text = '\n'.join(epitope_hlas)
        hla_ax = hla_axes[i]
        hla_ax.text(0.5,
                    0.5,
                    hla_text,
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=hla_ax.transAxes,
                    size=FONT_SIZE_MEDIUM)
        remove_spines_ticks_and_labels(hla_ax)

    return None


def show_num_patients(header_ax, patient_axes, events, epitope_summary):
    header_ax.text(0.5,
                   0,
                   '# Pt.\nw/HLA',
                   horizontalalignment='center',
                   verticalalignment='bottom',
                   transform=header_ax.transAxes,
                   size=FONT_SIZE_LARGE,
                   color=COLOR_BLACK,
                   fontweight='bold')
    remove_spines_ticks_and_labels(header_ax)
    header_ax.set_facecolor(COLOR_LIGHT_VIOLET)

    for i, event in enumerate(events):
        event_summary = epitope_summary.get(event)
        if not event_summary:
            return 'no event {} in epitope_summary'.format(event)

        epitopes = event_summary['epitopes']
        epitope_patients = [e['num_patients'] for e in epitopes[:2]]
        patient_text = '\n'.join(epitope_patients)
        patient_ax = patient_axes[i]
        patient_ax.text(0.5,
                        0.5,
                        patient_text,
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=patient_ax.transAxes,
                        size=FONT_SIZE_MEDIUM)
        remove_spines_ticks_and_labels(patient_ax)

    return None


def show_binding_affinity(header_ax, affinity_axes, events, epitope_summary):
    header_ax.text(0.5,
                   0,
                   'IC$_{50}$\n(nM)',
                   horizontalalignment='center',
                   verticalalignment='bottom',
                   transform=header_ax.transAxes,
                   size=FONT_SIZE_LARGE,
                   color=COLOR_BLACK,
                   fontweight='bold')
    remove_spines_ticks_and_labels(header_ax)
    header_ax.set_facecolor(COLOR_LIGHT_VIOLET)

    for i, event in enumerate(events):
        event_summary = epitope_summary.get(event)
        if not event_summary:
            return 'no event {} in epitope_summary'.format(event)

        epitopes = event_summary['epitopes']
        affinities = [
            str(int(round(e['binding_affinity']))) for e in epitopes[:2]
        ]
        affinity_text = '\n'.join(affinities)
        affinity_ax = affinity_axes[i]
        affinity_ax.text(0.5,
                         0.5,
                         affinity_text,
                         horizontalalignment='center',
                         verticalalignment='center',
                         transform=affinity_ax.transAxes,
                         size=FONT_SIZE_MEDIUM)
        remove_spines_ticks_and_labels(affinity_ax)

    return None


def show_epitope_bindings(events, epitope_summary, axes_by_name):
    epitope_header_ax, epitope_axes = axes_by_name['epitopes']
    hla_header_ax, hla_axes = axes_by_name['hlas']
    patient_header_ax, patient_axes = axes_by_name['num_patients']
    affinity_header_ax, affinity_axes = axes_by_name['affinity']

    error = show_epitope_names(epitope_header_ax, epitope_axes, events,
                               epitope_summary)
    if error:
        return error

    error = show_hla_names(hla_header_ax, hla_axes, events, epitope_summary)
    if error:
        return error

    error = show_num_patients(patient_header_ax, patient_axes, events,
                              epitope_summary)
    if error:
        return error

    error = show_binding_affinity(affinity_header_ax, affinity_axes, events,
                                  epitope_summary)
    if error:
        return error

    make_top_underline_header(epitope_header_ax, affinity_header_ax,
                              [hla_header_ax, patient_header_ax],
                              'Predicted HLA-epitope binding')

    return None


def make_plots(psi_data_by_event, epitope_summary, group_to_panel_type,
               out_file_name, has_prediction):
    num_events = len(psi_data_by_event)
    if num_events == 0:
        return 'no events to plot'

    fig_width = 8
    if has_prediction:
        fig_width = 12

    # 1.5 inches per row.
    # The header is given 2 rows worth of height when creating the grid.
    fig_height = 1.5 * (num_events + 2)
    fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=False)

    events = list(psi_data_by_event.keys())
    groups = psi_data_by_event[events[0]].columns

    axes_by_name, filler_rows = create_grid_of_axes(fig, num_events,
                                                    has_prediction)
    for filler_row in filler_rows:
        hide_ax(filler_row)

    error = make_violin_plots(events, psi_data_by_event, groups,
                              group_to_panel_type, axes_by_name)
    if error:
        return error

    error = make_shaded_dots_and_show_fold_change(events, epitope_summary,
                                                  axes_by_name)
    if error:
        return error

    if has_prediction:
        error = show_epitope_bindings(events, epitope_summary, axes_by_name)
        if error:
            return error

    plt.savefig(out_file_name)
    plt.close(fig)
    return None


def filter_events(events, epitope_summary):
    variants = list(events.keys())
    for variant in variants:
        event_list = events[variant]
        filtered = list()
        for event in event_list:
            # filter to events that are in epitope summary
            if event in epitope_summary:
                filtered.append(event)
            # only keep up to 10 events
            if len(filtered) == 10:
                break

        events[variant] = filtered


def exit_with_error(error):
    print(error, file=sys.stderr)
    sys.exit(1)


def main(args):
    event_type = args.splicing_event_type
    has_prediction = not args.no_prediction

    parameters, error = read_parameters(args.parameter_fin)
    if error:
        exit_with_error(error)

    group_to_panel_type = get_group_to_panel_type(parameters)
    groups = list(group_to_panel_type.keys())

    group_to_matrix_f_name, group_to_matrix_index_f_name, error = get_matrix_file_names(
        groups, parameters, event_type)
    if error:
        exit_with_error(error)

    events, error = get_events(parameters, args.screening_out_dir, event_type)
    if error:
        exit_with_error(error)

    epitope_summary, error = get_epitope_summary(args.screening_out_dir,
                                                 has_prediction, parameters,
                                                 event_type)
    if error:
        exit_with_error(error)

    filter_events(events, epitope_summary)

    psi_data_by_event, error = get_psi_data_by_event(
        group_to_matrix_f_name, group_to_matrix_index_f_name, events)
    if error:
        exit_with_error(error)

    error = make_plots(psi_data_by_event, epitope_summary, group_to_panel_type,
                       args.out_file_name, has_prediction)
    if error:
        exit_with_error(error)
