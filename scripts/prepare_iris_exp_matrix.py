import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description=('write input file for IRIS exp_matrix'))
    parser.add_argument('--out-manifest',
                        required=True,
                        help='path to write the list of gene expression files')
    parser.add_argument('--fpkm-files',
                        required=True,
                        nargs='+',
                        help='the fpkm files from cufflinks')

    args = parser.parse_args()
    return args


def prepare_iris_exp_matrix(out_manifest, fpkm_files):
    with open(out_manifest, 'wt') as out_handle:
        for file_name in fpkm_files:
            out_handle.write('{}\n'.format(file_name))


def main():
    args = parse_args()
    prepare_iris_exp_matrix(args.out_manifest, args.fpkm_files)


if __name__ == '__main__':
    main()
