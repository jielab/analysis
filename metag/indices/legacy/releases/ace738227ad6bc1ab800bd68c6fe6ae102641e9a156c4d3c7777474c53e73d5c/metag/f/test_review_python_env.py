"""Python selection regressions; no analysis jobs or library reads."""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SITE_PYTHON = Path('/work/sph-huangj/software/venvs/metag-20261115/bin/python')


class PythonEnvironment(unittest.TestCase):
    def selected(self, extra=None):
        env = dict(os.environ, submit='local', SUBMIT='local')
        for key in ('METAG_PYTHON', 'METAG_ENV_FILE', 'CONDA_PREFIX'):
            env.pop(key, None)
        env.update(extra or {})
        result = subprocess.run(
            ['bash', '-c', 'source "$1/f/env.sh"; printf "SELECTED=%s\\n" "$METAG_PYTHON"',
             'env-review', str(ROOT)], env=env, text=True, capture_output=True, check=True)
        return next(line.split('=', 1)[1] for line in result.stdout.splitlines()
                    if line.startswith('SELECTED='))

    def test_explicit_python_is_preserved(self):
        self.assertEqual(self.selected({'METAG_PYTHON': sys.executable}), sys.executable)

    def test_environment_file_python_is_preserved(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / 'worker.env'
            path.write_text('export METAG_PYTHON=' + shlex.quote(sys.executable) + '\n')
            self.assertEqual(self.selected({'METAG_ENV_FILE': str(path)}), sys.executable)

    def test_missing_site_python_does_not_block_fallback(self):
        selected = self.selected()
        if not os.access(SITE_PYTHON, os.X_OK):
            self.assertNotEqual(selected, str(SITE_PYTHON))
        self.assertIsNotNone(shutil.which(selected))


if __name__ == '__main__':
    unittest.main()
