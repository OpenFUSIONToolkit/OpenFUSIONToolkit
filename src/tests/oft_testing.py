import subprocess
import time
import os
import warnings
import pytest
import tempfile

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

    with tempfile.TemporaryFile() as out_file, tempfile.TemporaryFile() as err_file:
            pid = subprocess.Popen(
                command,
                shell=True,
                stdout=out_file,
                stderr=err_file
            )

            # Increase timeout if testing in debug configuration
            if os.environ.get('OFT_DEBUG_TEST', 0):
                timeout *= 4

            timeout *= 10
            while pid.poll() is None and timeout > 0:
                time.sleep(0.1)
                timeout -= 1

            if timeout <= 0:
                pid.kill()

            pid.wait()  # Ensure process is done

            # Read the outputs from temporary files
            out_file.seek(0)
            err_file.seek(0)
            std_out = out_file.read().decode(errors='replace')
            std_err = err_file.read().decode(errors='replace')
    errcode = pid.poll()
    #
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
        warnings.warn("WARNING: OFT emitted a warning during execution")
    if return_stdout:
        return True, std_out
    else:
        return True