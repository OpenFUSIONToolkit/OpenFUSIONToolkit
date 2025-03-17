import subprocess
import time
import os
import pytest

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
    if os.environ.get('OFT_DEBUG_TEST', 0):
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
        print("FAILED: OFT exited with non-zero error code!")
        if return_stdout:
            return False, std_out
        else:
            return False
    if std_out.find('ERROR') > -1:
        print("FAILED: detected OFT error!")
        if return_stdout:
            return False, std_out
        else:
            return False
    if std_out.find('WARNING') > -1:
        print("FAILED: detected OFT warning!")
        if return_stdout:
            return False, std_out
        else:
            return False
    if return_stdout:
        return True, std_out
    else:
        return True