import argparse
import os
import os.path


def parse_args():
    parser = argparse.ArgumentParser(
        description=('find which tasks were created by IRIS predict'))
    parser.add_argument('--out-list',
                        required=True,
                        help='path to write a file listing all created tasks')
    parser.add_argument(
        '--task-dir',
        required=True,
        help='directory where task files are expected to be found')
    parser.add_argument(
        '--splice-type',
        required=True,
        help='alternative splicing event type (expected in task file name)')

    args = parser.parse_args()
    return args


def count_iris_predict_tasks(out_list, task_dir, splice_type):
    file_names = os.listdir(task_dir)
    task_paths = list()
    base_prefix = 'pep2epitope_{}.'.format(splice_type)
    tiers = ['tier1', 'tier2tier3']
    prefixes = ['{}{}.'.format(base_prefix, tier) for tier in tiers]
    suffix = '.sh'
    for name in file_names:
        file_path = os.path.join(task_dir, name)
        if not name.endswith(suffix):
            continue

        for prefix in prefixes:
            if name.startswith(prefix):
                number_string = name[len(prefix):-len(suffix)]
                try:
                    int(number_string)
                except ValueError:
                    raise Exception('unexpected file: {}'.format(file_path))

                task_paths.append(file_path)
                continue

    if not task_paths:
        raise Exception('could not find any predict tasks')

    with open(out_list, 'wt') as out_handle:
        for task_path in task_paths:
            out_handle.write('{}\n'.format(task_path))


def main():
    args = parse_args()
    count_iris_predict_tasks(args.out_list, args.task_dir, args.splice_type)


if __name__ == '__main__':
    main()
