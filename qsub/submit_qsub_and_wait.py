import argparse
import subprocess
import time

import qsub


def local_subprocess_cmd_exec_func(tokens, _):
    process = subprocess.run(tokens, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if process.returncode != 0:
        print('returncode: {}, command: {}, stdout: {}, stderr: {}'.format(
            process.returncode, tokens, process.stdout, process.stderr))
        return None

    command_result = qsub.CommandResult()
    command_result.from_completed_process(process)
    return command_result


def submit_job(cmd):
    cmd = cmd.rstrip()
    print('executing: {}'.format(cmd))

    tokens = cmd.split(' ')
    command_result = local_subprocess_cmd_exec_func(tokens, None)
    if not command_result:
        print('error submitting job')
        return None

    job_name = ''
    return qsub.QsubJob(command_result, job_name, local_subprocess_cmd_exec_func)


def main():
    parser = argparse.ArgumentParser(
        description='execute qsub commands and wait for those jobs to complete')
    parser.add_argument('command_file', type=str, help='a file with 1 qsub command per line')
    parser.add_argument(
        '--poll-interval-seconds', type=int, default=30, help='how frequently to check job status')
    args = parser.parse_args()

    jobs = list()
    with open(args.command_file, 'rt') as f_handle:
        for cmd in f_handle:
            job = submit_job(cmd)
            if job:
                jobs.append(job)

    while jobs:
        time.sleep(args.poll_interval_seconds)
        print('checking {} job(s)'.format(len(jobs)))
        new_jobs = list()
        for job in jobs:
            if not job.is_finished():
                new_jobs.append(job)
            else:
                print('finished j_id: {}'.format(job.j_id))

        jobs = new_jobs


if __name__ == '__main__':
    main()
