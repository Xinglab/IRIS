import argparse
import subprocess
import time

SLEEP_SECONDS_BETWEEN_CHECKS = 30

class Job(object):
    def __init__(self, submit_process):
        self.args = submit_process.args
        self.submit_return_code = submit_process.returncode
        self.submit_out = submit_process.stdout
        self.submit_err = submit_process.stderr

    def is_finished(self):
        if self.submit_return_code != 0:
            print('is_finished: submit_return_code: {}'.format(self.submit_return_code))
            print('out: {}'.format(self.submit_out.decode()))
            print('err: {}'.format(self.submit_err.decode()))
            return True

        j_id = self.get_j_id()
        if j_id is None:
            return False

        return j_id_is_finished(j_id)

    # Your job-array 11173394.1-954:1 ("IRIS_pep2epitope") has been submitted
    def get_j_id(self):
        tokens = self.submit_out.decode().split(' ')
        for token in tokens:
            if token and token[0].isdigit():
                return token.split('.')[0]

        print('get_j_id: could not find j_id in: {}'.format(tokens))

def j_id_is_finished(j_id):
    poll_process = subprocess.run(['qstat', '-j', j_id, '-xml'],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if poll_process.returncode != 0:
        print('j_id_is_finished: poll_process.returncode: {}'.format(poll_process.returncode))
        print('out: {}'.format(poll_process.stdout.decode()))
        print('err: {}'.format(poll_process.stderr.decode()))
        return True # give up on job since cannot poll

    return '<unknown_jobs' in poll_process.stdout.decode()


def submit_job(cmd):
    cmd = cmd.rstrip()
    print('executing: {}'.format(cmd))
    return Job(subprocess.run(
        cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE))


def main():
    parser = argparse.ArgumentParser(
        description='execute qsub commands and wait for those jobs to complete')
    parser.add_argument('command_file', type=str,
                        help='a file with 1 qsub command per line')
    args = parser.parse_args()

    jobs = list()
    with open(args.command_file, 'rt') as f_handle:
        for cmd in f_handle:
            jobs.append(submit_job(cmd))

    while jobs:
        time.sleep(SLEEP_SECONDS_BETWEEN_CHECKS)
        print('checking {} job(s)'.format(len(jobs)))
        new_jobs = list()
        for job in jobs:
            if not job.is_finished():
                new_jobs.append(job)
            else:
                print('finished j_id: {}'.format(job.get_j_id()))

        jobs = new_jobs


if __name__ == '__main__':
    main()
