import subprocess
import multiprocessing
import queue
import signal
import time
import os
import sys
import warnings
import pytest

# Make the OFT Python package (<build_dir>/python) importable by all tests
oft_python_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python'))
if oft_python_dir not in sys.path:
    sys.path.append(oft_python_dir)


def env_flag(name):
    '''Return `True` if environment variable `name` is set to an enabled value ("1", "true", "yes", or "on")

    Unset, empty, and unrecognized values count as disabled, so a malformed setting never raises.
    '''
    return os.environ.get(name, '').strip().lower() in ('1', 'true', 'yes', 'on')


def run_command(command, cwd=None):
    pid = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
    outs, errs = pid.communicate()
    errcode = pid.poll()
    return outs, errs, errcode


def kill_process_group(pid):
    '''Kill the process group led by `pid` (started with `start_new_session=True`), if any of it is still running'''
    try:
        os.killpg(pid.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def run_OFT(command, nproc, timeout, return_stdout=False, expect_abort=None):
    '''Run OFT `command` on `nproc` processes with a `timeout` (in seconds) and check that it succeeded

    Returns `True` if the run succeeded (and its standard output if `return_stdout`). A run fails if it times
    out, exits with a non-zero error code, or reports an error or warning. If `expect_abort` is given, the run
    must instead exit with a non-zero error code and report `expect_abort` in its error output (eg. from
    `oft_abort`).
    '''
    if nproc > 1:
        if env_flag('OFT_HAVE_MPI'):
            command = "mpirun --map-by :OVERSUBSCRIBE -np {0} {1}".format(nproc, command)
        else:
            pytest.skip("Not compiled with MPI")
    # Run in a new session so a timeout can kill the whole process group (shell, mpirun, and ranks). The group
    # is then outside the terminal's foreground group, so an interrupt (Ctrl-C) must also kill it explicitly
    pid = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           start_new_session=True)
    # Increase timeout if testing in debug configuration
    if env_flag('OFT_DEBUG_TEST'):
        timeout *= 4
    # Wait for process to complete or timeout, reading the output as it arrives; waiting without reading
    # deadlocks once the output fills the pipe buffer (64 KB)
    timed_out = False
    try:
        outs, errs = pid.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        kill_process_group(pid)
        outs, errs = pid.communicate()
    except BaseException:
        kill_process_group(pid)
        raise
    errcode = pid.poll()
    #
    std_out = outs.decode()
    std_err = errs.decode()
    print("========== OFT STD OUTPUT ==========")
    print(std_out)
    print("========== OFT ERR OUTPUT ==========")
    print(std_err)
    print("ERRCODE = {0}".format(errcode))
    print("========== END OFT OUTPUT ==========")
    #
    if timed_out:
        print("FAILED: OFT timed out after {0} s".format(timeout))
        if return_stdout:
            return False, std_out
        else:
            return False
    if std_out.find('Not compiled with PETSc') > -1:
        pytest.skip("Not compiled with PETSc")
        if return_stdout:
            return True, std_out
        else:
            return True
    if std_out.find('SKIP TEST') > -1:
        pytest.skip("Configuration not supported")
        if return_stdout:
            return True, std_out
        else:
            return True
    if expect_abort is not None:
        result = True
        if errcode == 0:
            print("FAILED: OFT exited normally, but an error was expected!")
            result = False
        elif std_err.find(expect_abort) < 0:
            print("FAILED: expected OFT error not found!")
            result = False
        if return_stdout:
            return result, std_out
        else:
            return result
    if errcode != 0:
        if (errcode == 143) or (errcode == 15):
            warnings.warn("WARNING: OFT exited with error code 143 (external SIGTERM)")
        else:
            print("FAILED: OFT exited with non-zero error code!")
            if return_stdout:
                return False, std_out
            else:
                return False
    if std_out.find('ERROR:') > -1:
        print("FAILED: detected OFT error!")
        if return_stdout:
            return False, std_out
        else:
            return False
    if std_out.find('WARNING:') > -1:
        print("FAILED: detected OFT warning!")
        if return_stdout:
            return False, std_out
        else:
            return False
    if return_stdout:
        return True, std_out
    else:
        return True


def mp_run(target, args, cwd, timeout=180):
    '''Run `target(*args, mp_q)` in a separate process from `cwd` and return the first result put on `mp_q`

    Only one `OFT_env` may exist per process, so tests using the Python interface run each case in its own
    process. Returns `None` if the process does not finish within `timeout` seconds or produces no result.
    The caller's working directory is restored before returning.
    '''
    # Increase timeout if testing in debug configuration
    if env_flag('OFT_DEBUG_TEST'):
        timeout *= 4
    deadline = time.time() + timeout
    orig_dir = os.getcwd()
    os.chdir(cwd)  # The child process starts in the current directory
    try:
        mp_q = multiprocessing.Queue()
        p = multiprocessing.Process(target=target, args=args + (mp_q,))
        p.start()
        # Read the result while waiting: a child cannot exit until a large result has been read from the queue
        test_result = None
        while True:
            try:
                test_result = mp_q.get(timeout=0.5)
                break
            except queue.Empty:
                pass
            if not p.is_alive():
                # Exited without a result (eg. crash); allow for a result still in transit
                try:
                    test_result = mp_q.get(timeout=1)
                except queue.Empty:
                    print("Failed to get output")
                break
            if time.time() > deadline:
                break
        # Wait for the child to finish (eg. writing coverage data), still within the timeout
        p.join(max(deadline - time.time(), 0.0))
        if p.is_alive():
            print("Timeout reached")
            p.terminate()
            p.join()
            return None
        return test_result
    finally:
        os.chdir(orig_dir)


def approx_list(list, rel=None, abs=None):
    return [pytest.approx(val, rel=rel, abs=abs) if val is not None else None for val in list]