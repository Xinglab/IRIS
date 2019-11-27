import os
import sys
import tempfile
import time
import unittest

import submit_qsub_and_wait


def _write_lines(f_name, lines):
    with open(f_name, 'wt') as f_h:
        for line in lines:
            f_h.write('{}\n'.format(line))


class TestSubmitQsubAndWait(unittest.TestCase):
    def test(self):
        temp_f_name_1 = None
        temp_f_name_2 = None
        try:
            with tempfile.NamedTemporaryFile(delete=False) as temp_f_handle:
                temp_f_name_1 = temp_f_handle.name

            with tempfile.NamedTemporaryFile(delete=False) as temp_f_handle:
                temp_f_name_2 = temp_f_handle.name

            self._test(temp_f_name_1, temp_f_name_2, 1, 60)
            self._test(temp_f_name_1, temp_f_name_2, 30, 90)
        finally:
            if temp_f_name_1 is not None:
                os.remove(temp_f_name_1)

            if temp_f_name_2 is not None:
                os.remove(temp_f_name_2)

    def _test(self, f_name_1, f_name_2, sleep_seconds, max_seconds):
        lines_1 = [
            '#!/bin/bash',
            'sleep {}'.format(sleep_seconds),
        ]
        _write_lines(f_name_1, lines_1)
        lines_2 = [
            'qsub {}'.format(f_name_1),
            'qsub {}'.format(f_name_1),
        ]
        _write_lines(f_name_2, lines_2)

        sys.argv = ['submit_qsub_and_wait.py', f_name_2, '--poll-interval-seconds', '5']

        begin = time.time()
        submit_qsub_and_wait.main()
        end = time.time()

        elapsed_seconds = end - begin
        self.assertTrue(elapsed_seconds >= sleep_seconds)
        self.assertTrue(elapsed_seconds <= max_seconds)


if __name__ == '__main__':
    unittest.main(verbosity=2)
