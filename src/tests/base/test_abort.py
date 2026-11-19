import os
import re
import subprocess
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import env_flag, kill_process_group

# Basic template for input file
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
 omp_nthreads={0}
/

&test_abort_options
 test_id={1}
/
"""


def run_abort(test_id, nthreads=4, nproc=1, timeout=30):
    '''Run test_abort, which always ends through oft_abort, and return (exit code, combined output)'''
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nthreads, test_id))
    command = "./test_abort"
    if nproc > 1:
        if env_flag('OFT_HAVE_MPI'):
            command = "mpirun --map-by :OVERSUBSCRIBE -np {0} {1}".format(nproc, command)
        else:
            pytest.skip("Not compiled with MPI")
    pid = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           start_new_session=True)
    try:
        outs, errs = pid.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        kill_process_group(pid)
        pid.communicate()
        pytest.fail("oft_abort did not end the run within {0} s".format(timeout))
    except BaseException:
        kill_process_group(pid)
        raise
    # Garbled (non-UTF-8) output was one symptom of concurrent aborts, so decode leniently
    output = outs.decode(errors='replace') + errs.decode(errors='replace')
    print(output)
    print("ERRCODE = {0}".format(pid.returncode))
    assert output.find('Run continued after oft_abort') < 0
    return pid.returncode, output


# Every thread of a parallel loop calls oft_abort: one error block, non-zero exit
@pytest.mark.coverage
def test_abort_parallel_loop():
    errcode, output = run_abort(1)
    assert errcode == 99
    assert output.count('ERROR:') == 1


# oft_abort outside a parallel region
@pytest.mark.coverage
def test_abort_serial():
    errcode, output = run_abort(2)
    assert errcode == 99
    assert output.count('ERROR:') == 1


# Abort callback that calls oft_abort again must stop the run, not deadlock
@pytest.mark.coverage
def test_abort_reentrant():
    errcode, output = run_abort(3)
    assert errcode == 99
    assert output.count('ERROR:') == 1
    assert output.find('Abort called from abort callback') < 0


# Parallel loop abort on 2 MPI ranks: at most one error block per rank, non-zero exit
@pytest.mark.mpi
@pytest.mark.coverage
def test_abort_parallel_loop_mpi():
    errcode, output = run_abort(1, nthreads=2, nproc=2)
    assert errcode != 0
    rank_counts = [len(re.findall(r'\[\s*{0}\] ERROR:'.format(rank), output)) for rank in range(2)]
    assert 1 <= sum(rank_counts)
    assert max(rank_counts) <= 1
