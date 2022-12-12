import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description=('write input file for IRIS sjc_matrix'))
    parser.add_argument('--sj-out',
                        required=True,
                        help='path to write the list of SJ count files')
    parser.add_argument('--sj-files',
                        required=True,
                        nargs='+',
                        help='the SJ count files from extract_sjc')

    args = parser.parse_args()
    return args


def prepare_iris_sjc_matrix(sj_out, sj_files):
    with open(sj_out, 'wt') as out_handle:
        for file_name in sj_files:
            out_handle.write('{}\n'.format(file_name))


def main():
    args = parse_args()
    prepare_iris_sjc_matrix(args.sj_out, args.sj_files)


if __name__ == '__main__':
    main()
