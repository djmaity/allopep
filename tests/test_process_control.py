"""Exercise real subprocess deadlines, logs, and detached descendant cleanup."""

import json
import sys
import tempfile
import time
import unittest
from pathlib import Path

import psutil

from process_control import run_logged


class ProcessControlTests(unittest.TestCase):
    def test_logs_success_and_nonzero_exit(self):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / 'stage.log'
            run_logged([sys.executable, '-c', 'print("saved output")'],
                       log_file=log, deadline=time.monotonic() + 10)
            self.assertEqual(log.read_text(), 'saved output\n')
            with self.assertRaisesRegex(RuntimeError, 'code 7'):
                run_logged([sys.executable, '-c', 'print("failure detail"); raise SystemExit(7)'],
                           log_file=log, deadline=time.monotonic() + 10)
            self.assertIn('failure detail', log.read_text())

    def test_timeout_kills_parent_and_detached_descendant(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            script = '\n'.join([
                'import json, os, subprocess, sys, time',
                'child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"], '
                'start_new_session=True)',
                f'with open({str(root / "pids.json")!r}, "w") as handle:',
                '    json.dump([os.getpid(), child.pid], handle)',
                'print("running", flush=True)',
                'time.sleep(60)',
            ])
            with self.assertRaisesRegex(TimeoutError, 'timed out'):
                run_logged([sys.executable, '-c', script], log_file=root / 'stage.log',
                           deadline=time.monotonic() + 1)
            self.assertIn('running', (root / 'stage.log').read_text())
            for pid in json.loads((root / 'pids.json').read_text()):
                try:
                    process = psutil.Process(pid)
                    self.assertEqual(process.status(), psutil.STATUS_ZOMBIE)
                except psutil.NoSuchProcess:
                    pass

    def test_expired_deadline_does_not_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(TimeoutError, 'before launching'):
                run_logged(['/nonexistent'], log_file=Path(directory) / 'stage.log',
                           deadline=time.monotonic() - 1)


if __name__ == '__main__':
    unittest.main()
