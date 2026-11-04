import subprocess
import multiprocessing
import time
import os
import sys
import warnings
import pytest

# Make the OFT Python package (<build_dir>/python) importable by all tests
oft_python_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python'))
if oft_python_dir not in sys.path:
    sys.path.append(oft_python_dir)

def run_command(command, cwd=None):
    pid = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
    outs, errs = pid.communicate()
    errcode = pid.poll()
    return outs, errs, errcode


def run_OFT(command, nproc, timeout, return_stdout=False):
    if nproc > 1:
        if int(os.environ.get('OFT_HAVE_MPI', 0)) == 1:
            command = "mpirun --map-by :OVERSUBSCRIBE -np {0} {1}".format(nproc, command)
        else:
            pytest.skip("Not compiled with MPI")
    pid = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Increase timeout if testing in debug configuration
    if int(os.environ.get('OFT_DEBUG_TEST', 0)) == 1:
        timeout *= 4
    # Wait for process to complete or timeout
    timeout *= 10
    while (pid.poll() is None) and (timeout > 0):
        time.sleep(0.1)
        timeout -= 1
    if timeout <= 0:
        pid.kill()
    outs, errs = pid.communicate()
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
    '''
    # Increase timeout if testing in debug configuration
    if int(os.environ.get('OFT_DEBUG_TEST', 0)) == 1:
        timeout *= 4
    os.chdir(cwd)
    mp_q = multiprocessing.Queue()
    p = multiprocessing.Process(target=target, args=args + (mp_q,))
    p.start()
    start = time.time()
    while time.time() - start <= timeout:
        if not p.is_alive():
            break
        time.sleep(.5)
    else: # Reached timeout
        print("Timeout reached")
        p.terminate()
        p.join()
        return None
    # Completed successfully?
    try:
        test_result = mp_q.get(timeout=5)
    except:
        print("Failed to get output")
        return None
    p.join()
    return test_result


def approx_list(list, rel=None, abs=None):
    return [pytest.approx(val, rel=rel, abs=abs) if val is not None else None for val in list]