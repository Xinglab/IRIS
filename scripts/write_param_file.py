import argparse
import os
import os.path


def parse_args():
    parser = argparse.ArgumentParser(
        description=('write the parameter file for IRIS screen'))
    parser.add_argument('--out-path',
                        required=True,
                        help='path to write parameter file')
    parser.add_argument('--group-name',
                        required=True,
                        help='name to use for sub directory in IRIS_data/db/')
    parser.add_argument('--iris-db',
                        required=True,
                        help='/path/to/IRIS_data/db')
    parser.add_argument(
        '--psi-p-value-cutoffs',
        required=True,
        help=('comma separated p-value cutoffs for PSI-based statistical tests'
              ' (tissue-matched normal, tumor, normal)'))
    parser.add_argument(
        '--sjc-p-value-cutoffs',
        required=True,
        help=('comma separated p-value cutoffs for SJC-based statistical tests'
              ' (tissue-matched normal, tumor, normal)'))
    parser.add_argument('--delta-psi-cutoffs',
                        required=True,
                        help=('comma separated minimum required delta PSIs'
                              ' (tissue-matched normal, tumor, normal)'))
    parser.add_argument('--fold-change-cutoffs',
                        required=True,
                        help=('comma separated minimum required fold changes'
                              ' (tissue-matched normal, tumor, normal)'))
    parser.add_argument(
        '--group-count-cutoffs',
        required=True,
        help=('comma separated minimum counts of reference groups that'
              ' need to meet other requirements'
              ' (tissue-matched normal, tumor, normal)'))
    parser.add_argument(
        '--reference-names-tissue-matched-normal',
        required=True,
        help='comma separated reference groups for tissue-matched normal')
    parser.add_argument('--reference-names-tumor',
                        required=True,
                        help='comma separated reference groups for tumor')
    parser.add_argument('--reference-names-normal',
                        required=True,
                        help='comma separated reference groups for normal')
    parser.add_argument('--comparison-mode',
                        required=True,
                        choices=['group', 'individual'],
                        help=('mode for statistical test'
                              ' (group requires at least 2 input samples)'))
    parser.add_argument('--statistical-test-type',
                        required=True,
                        choices=['parametric', 'nonparametric'],
                        help='type of statistical test')
    parser.add_argument('--use-ratio',
                        action='store_true',
                        help='use ratio instead of count for group cutoffs')
    parser.add_argument('--blocklist-file', help='list of AS events to remove')
    parser.add_argument('--mapability-bigwig',
                        help='allows evaluatio of splice region mapability')
    parser.add_argument('--reference-genome',
                        help='required for IRIS translate')

    args = parser.parse_args()
    check_file_exists(args.blocklist_file, parser)
    check_file_exists(args.mapability_bigwig, parser)
    check_file_exists(args.reference_genome, parser)

    args.psi_p_value_cutoffs = parse_floats(args.psi_p_value_cutoffs)
    args.sjc_p_value_cutoffs = parse_floats(args.sjc_p_value_cutoffs)
    args.delta_psi_cutoffs = parse_floats(args.delta_psi_cutoffs)
    args.fold_change_cutoffs = parse_floats(args.fold_change_cutoffs)
    args.group_count_cutoffs = parse_floats(args.group_count_cutoffs)

    if not (1 <= len(args.psi_p_value_cutoffs) <= 3):
        parser.error('must give 1 to 3 cutoffs')

    expected_len = len(args.psi_p_value_cutoffs)
    expected_non_none = [x is not None for x in args.psi_p_value_cutoffs]
    if not (expected_non_none[0] or
            ((expected_len == 3) and expected_non_none[2])):
        parser.error('must provide values for at least one of'
                     ' tissue-matched-normal or normal')

    for name, values in [('sjc_p_value_cutoffs', args.sjc_p_value_cutoffs),
                         ('delta_psi_cutoffs', args.delta_psi_cutoffs),
                         ('fold_change_cutoffs', args.fold_change_cutoffs),
                         ('group_count_cutoffs', args.group_count_cutoffs)]:
        if len(values) != expected_len:
            parser.error('{} has len {}, but expected {}'.format(
                name, len(values), expected_len))

        for i, value in enumerate(values):
            is_non_none = value is not None
            if is_non_none != expected_non_none[i]:
                expected = 'non-None' if expected_non_none[i] else 'None'
                parser.error('{} value {} was {} when {} was expected'.format(
                    name, i, value, expected))

    db_names = get_db_names(args.iris_db, parser)
    args.reference_names_tissue_matched_normal = parse_reference_names(
        args.reference_names_tissue_matched_normal, args.group_count_cutoffs,
        0, 'tissue-matched-normal', db_names, args.use_ratio, parser)
    args.reference_names_tumor = parse_reference_names(
        args.reference_names_tumor, args.group_count_cutoffs, 1, 'tumor',
        db_names, args.use_ratio, parser)
    args.reference_names_normal = parse_reference_names(
        args.reference_names_normal, args.group_count_cutoffs, 2, 'normal',
        db_names, args.use_ratio, parser)

    return args


def get_db_names(db_path, parser):
    if not os.path.exists(db_path):
        parser.error('{} does not exist'.format(db_path))
    if not os.path.isdir(db_path):
        parser.error('{} is not a directory'.format(db_path))

    db_names = list()
    dir_entries = os.listdir(db_path)
    for db_name in dir_entries:
        full_path = os.path.join(db_path, db_name)
        if os.path.isdir(full_path):
            db_names.append(db_name)

    return db_names


def parse_floats(floats_str):
    parts = floats_str.split(',')
    stripped = [x.strip() for x in parts]
    results = list()
    for string in stripped:
        if len(string) != 0:
            float_value = float(string)
            results.append(float_value)
        else:
            results.append(None)

    return results


def parse_reference_names(names_str, group_cutoffs, cutoff_index, group_name,
                          db_names, use_ratio, parser):
    if len(group_cutoffs) <= cutoff_index:
        parser.error(
            'missing list of reference groups for {}'.format(group_name))

    cutoff = group_cutoffs[cutoff_index]
    parts = names_str.split(',')
    stripped = [x.strip() for x in parts]
    if cutoff is None:
        return list()  # this group is skipped

    if use_ratio:
        if not (0 <= cutoff <= 1):
            parser.error('cutoff for {} was {} with use_ratio'.format(
                group_name, cutoff))
    elif cutoff > len(stripped):
        parser.error('{} cutoff is {}, but only {} references'.format(
            group_name, cutoff, len(stripped)))

    for name in stripped:
        if name not in db_names:
            parser.error('reference {} in {} not found in db/'.format(
                name, group_name))

    return stripped


def check_file_exists(file_path, parser):
    if file_path is None:
        return

    if not os.path.isfile(file_path):
        parser.error('{} does not exist'.format(file_path))


def write_file_line_or_empty_line(out_handle, maybe_file):
    if maybe_file:
        out_handle.write('{}\n'.format(maybe_file))
    else:
        out_handle.write('\n')


def write_param_file(args):
    with open(args.out_path, 'wt') as out_handle:
        out_handle.write('{}\n'.format(args.group_name))
        abs_db_path = os.path.abspath(args.iris_db)
        out_handle.write('{}\n'.format(abs_db_path))
        references_by_i = [
            args.reference_names_tissue_matched_normal,
            args.reference_names_tumor, args.reference_names_normal
        ]
        for i, psi_cutoff in enumerate(args.psi_p_value_cutoffs):
            if psi_cutoff is None:
                out_handle.write('\n')
                continue

            cutoffs = [
                psi_cutoff, args.delta_psi_cutoffs[i],
                args.fold_change_cutoffs[i], args.sjc_p_value_cutoffs[i],
                args.group_count_cutoffs[i]
            ]
            cutoffs = [str(x) for x in cutoffs]
            references = references_by_i[i]
            out_handle.write('{} {}\n'.format(','.join(cutoffs),
                                              ','.join(references)))

        out_handle.write('{} {}\n'.format(args.comparison_mode,
                                          args.statistical_test_type))
        out_handle.write('{}\n'.format('True' if args.use_ratio else 'False'))
        write_file_line_or_empty_line(out_handle, args.blocklist_file)
        write_file_line_or_empty_line(out_handle, args.mapability_bigwig)
        write_file_line_or_empty_line(out_handle, args.reference_genome)


def main():
    args = parse_args()
    write_param_file(args)


if __name__ == '__main__':
    main()
