import argparse
import os
import os.path


def parse_args():
    parser = argparse.ArgumentParser(
        description=('write input files for IRIS format'))
    parser.add_argument('--matrix-out',
                        required=True,
                        help='path to write the list of matrix directories')
    parser.add_argument('--sample-out',
                        required=True,
                        help='path to write the list of BAM lists')
    parser.add_argument('--summaries',
                        required=True,
                        nargs='+',
                        help='the summary files from the matrix directories')

    args = parser.parse_args()
    return args


def prepare_iris_format(matrix_out, sample_out, summaries):
    with open(matrix_out, 'wt') as matrix_out_handle:
        with open(sample_out, 'wt') as sample_out_handle:
            prepare_iris_format_with_handles(matrix_out_handle,
                                             sample_out_handle, summaries)


def prepare_iris_format_with_handles(matrix_out_handle, sample_out_handle,
                                     summaries):
    for summary in summaries:
        matrix_dir_path = os.path.dirname(summary)
        matrix_dir_name = os.path.basename(matrix_dir_path)
        matrix_dir_name_suffix = '.matrix'
        if not matrix_dir_name.endswith(matrix_dir_name_suffix):
            raise Exception('unexpected directory name for {}'.format(summary))

        matrix_dir_name_prefix = matrix_dir_name[:-len(matrix_dir_name_suffix)]
        matrix_dir_parent_dir_path = os.path.dirname(matrix_dir_path)
        sample_list_name = '{}_rmatspost_list.txt'.format(
            matrix_dir_name_prefix)
        sample_path = os.path.join(matrix_dir_parent_dir_path,
                                   sample_list_name)

        matrix_out_handle.write('{}\n'.format(matrix_dir_path))
        sample_out_handle.write('{}\n'.format(sample_path))


def main():
    args = parse_args()
    prepare_iris_format(args.matrix_out, args.sample_out, args.summaries)


if __name__ == '__main__':
    main()
