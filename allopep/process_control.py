"""Run logged pipeline stages with a deadline and descendant cleanup."""

import os
import signal
import subprocess
import time
from pathlib import Path

import psutil


def _remember_children(pid, children):
    try:
        for child in psutil.Process(pid).children(recursive=True):
            children[(child.pid, child.create_time())] = child
    except psutil.Error:
        pass


def _stop_tree(process, children):
    # Ray may start descendants in a separate session. Remember those processes
    # while the parent is alive, and use Process identities to avoid PID reuse.
    _remember_children(process.pid, children)
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    for child in children.values():
        try:
            child.terminate()
        except psutil.Error:
            pass
    _, alive = psutil.wait_procs(list(children.values()), timeout=2)
    for child in alive:
        try:
            child.kill()
        except psutil.Error:
            pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def run_logged(command, *, log_file, deadline, cwd=None, env=None):
    """Run until an absolute monotonic deadline; retain output even on failure."""
    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    if time.monotonic() >= deadline:
        raise TimeoutError(f'Round deadline expired before launching {command[0]}')
    children = {}
    with log_file.open('w', encoding='utf-8') as handle:
        process = subprocess.Popen(
            command, cwd=cwd, env=env, stdout=handle,
            stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while True:
                _remember_children(process.pid, children)
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(f'Round timed out; see {log_file}')
                try:
                    returncode = process.wait(timeout=min(0.2, remaining))
                    break
                except subprocess.TimeoutExpired:
                    continue
        finally:
            _stop_tree(process, children)
    if returncode:
        # Logs may be large; read only their tail for the exception.
        with log_file.open('rb') as handle:
            handle.seek(max(0, log_file.stat().st_size - 3000))
            tail = handle.read().decode('utf-8', errors='replace')
        raise RuntimeError(f'{command[0]} exited with code {returncode}; '
                           f'see {log_file}\n{tail}')
