import argparse
import os
import os.path


def parse_args():
    parser = argparse.ArgumentParser(
        description=('determine the set of read lengths based on'
                     ' the rmats output file names'))
    parser.add_argument(
        '--parent-dir',
        required=True,
        help='path of directory which contains 1 directory per read length')
    parser.add_argument('--run-name',
                        required=True,
                        help='prefix used to name output files')
    parser.add_argument('--out',
                        required=True,
                        help='path to write read lengths')

    args = parser.parse_args()
    return args


def check_read_lengths(parent_dir, run_name, out):
    file_names = os.listdir(parent_dir)
    prefix = '{}.RL'.format(run_name)
    read_lengths = list()
    for file_name in file_names:
        file_path = os.path.join(parent_dir, file_name)
        if not (os.path.isdir(file_path) and file_name.startswith(prefix)):
            continue

        suffix = file_name[len(prefix):]
        try:
            read_length = int(suffix)
        except ValueError:
            print('ignoring: {}'.format(file_path))
            continue

        read_lengths.append(suffix)

    if not read_lengths:
        raise Exception('no read lengths found in {}'.format(parent_dir))

    with open(out, 'wt') as out_handle:
        for read_length in read_lengths:
            out_handle.write('{}\n'.format(read_length))


def main():
    args = parse_args()
    check_read_lengths(args.parent_dir, args.run_name, args.out)


if __name__ == '__main__':
    main()
