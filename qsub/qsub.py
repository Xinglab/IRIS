import os

import bs4


class CommandResult(object):
    def __init__(self):
        self.status_code = None
        self.out = None
        self.err = None

    def from_completed_process(self, completed_process):
        self.status_code = completed_process.returncode
        self.out = completed_process.stdout
        self.err = completed_process.stderr


class QsubJob(object):
    """
    cmd_execute_func(cmd_tokens, exec_func_ref_data) -> CommandResult

    cmd_execute_func and log_error_func must be defined at the top level of a module so that
    QsubJob can be pickled.

    The execute function allows either running a local qsub or qsub on a remote host.
    """
    def __init__(self,
                 command_result,
                 job_name,
                 cmd_execute_func,
                 exec_func_ref_data=None,
                 out_dir=None,
                 log_error_func=print):
        self.job_name = job_name
        self._cmd_execute_func = cmd_execute_func
        self._exec_func_ref_data = exec_func_ref_data
        self._log_error_func = log_error_func

        self._status = 'running'

        self.j_id = _extract_qsub_j_id(command_result.out.decode())
        self._set_output_file_names(out_dir)

    def _set_output_file_names(self, out_dir):
        if out_dir is None:
            self.qsub_out = None
            self.qsub_err = None
            return

        out_f_name = '{}.o{}'.format(self.job_name, self.j_id)
        err_f_name = '{}.e{}'.format(self.job_name, self.j_id)
        self.qsub_out = os.path.join(out_dir, out_f_name)
        self.qsub_err = os.path.join(out_dir, err_f_name)

    def _execute_cmd(self, cmd_tokens):
        return self._cmd_execute_func(cmd_tokens, self._exec_func_ref_data)

    def get_status(self):
        if self._status == 'finished':
            return self._status

        if self.is_finished():
            self._status = 'finished'
            return self._status

        return self._status

    def is_finished(self):
        if not self.j_id:
            self._log_error_func('no j_id when checking if qsub job is finished')
            return False

        qstat_command_tokens = ['qstat', '-j', self.j_id, '-xml']
        qstat_process = self._execute_cmd(qstat_command_tokens)
        if not qstat_process:
            return None

        qstat_soup = _make_soup(qstat_process.out.decode())
        if _qstat_output_has_job_details(qstat_soup, self.j_id):
            return False

        if _qstat_output_lists_job_as_unknown(qstat_soup, self.j_id):
            return True

        self._log_error_func('unexpected output from: {}\n{}'.format(qstat_command_tokens,
                                                                     qstat_soup))
        return True

    def get_stdout(self):
        if self.qsub_out is None:
            return None

        return self._cat_file(self.qsub_out)

    def get_stderr(self):
        if self.qsub_err is None:
            return None

        return self._cat_file(self.qsub_err)

    def _cat_file(self, f_name):
        cat_command_tokens = ['cat', f_name]
        process = self._execute_cmd(cat_command_tokens)
        if not process:
            return None

        return process.out


def _make_soup(output):
    return bs4.BeautifulSoup(output, 'html.parser')


def _extract_qsub_j_id(out):
    """
    Find j_id in out that looks like:
    Your job {j_id}.{array_details} ("{job_name}") has been submitted
    """
    tokens = out.split(' ')
    for token in tokens:
        if token and token[0].isdigit():
            return token.split('.')[0]

    return None


def _qstat_output_has_job_details(soup, j_id):
    """
    Return True if j_id is a JB_job_number in soup that looks like:
    <detailed_job_info ...>
      <djob_info>
        <element>
          <JB_job_number>11693229</JB_job_number>
    """
    djobs = soup.find_all('djob_info')
    if len(djobs) != 1:
        return False

    numbers = djobs[0].find_all('jb_job_number')
    for number in numbers:
        if number.string == j_id:
            return True

    return False


def _qstat_output_lists_job_as_unknown(soup, j_id):
    """
    Return True if j_id is an st_name in soup that looks like:
    <unknown_jobs ...>
      <element>
        <ST_name>123</ST_name>
    """
    unks = soup.find_all('unknown_jobs')
    if len(unks) != 1:
        return False

    names = unks[0].find_all('st_name')
    for name in names:
        if name.string == j_id:
            return True

    return False
