import os
import shutil
import subprocess

def run_command(command):
    pid = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outs, errs = pid.communicate()
    errcode = pid.poll()
    return outs.decode(), errs.decode(), errcode

def get_prereqs(filename):
    preqs = []
    known_libs = {}
    output, error, errcode = run_command('readelf -d {0}'.format(filename))
    if errcode != 0:
        return preqs, known_libs
    #print(filename)
    for line in output.splitlines():
        if line.find('Shared library:') > 0:
            preqs.append(line.split('[')[1].split(']')[0])
    #for preq in preqs:
    #    print('  '+preq)
    ldd_output, error, errcode = run_command('ldd -d {0}'.format(filename))
    if errcode == 0:
        for line in ldd_output.splitlines():
            if line.find('=>') > 0:
                known_libs[line.split('=>')[0].strip()] = line.split('=>')[1].split('(')[0].strip()
    return preqs, known_libs

print('Copying required dynamic libraries')
output, error, errcode = run_command('ldconfig -p')
ld_confoutput = output.splitlines()

MKL_ROOT = None
for i in range(2):
    full_prequisites = {}
    known_resolutions = {}
    # Gather prerequisites
    patch_files = [f for f in os.listdir() if os.path.isfile(f)]
    for filename in patch_files:
        preqs, new_known = get_prereqs(filename)
        for preq in preqs:
            if preq not in full_prequisites:
                full_prequisites[preq] = None
        for known_lib in new_known:
            if known_lib not in known_resolutions:
                known_resolutions[known_lib] = new_known[known_lib]

    # Find system libraries that are not typically in a base install
    filter_list = [
        'libgomp.','libgfortran.','libquadmath.',
        'libmkl','libifport','libifcoremt','libimf','libsvml','libintlc','libirng','libiomp5'
    ]
    keep_keys = []
    for key in full_prequisites:
        for known_lib in filter_list:
            if key.startswith(known_lib):
                keep_keys.append(key)
    prequisites = {}
    for key in keep_keys:
        prequisites[key] = full_prequisites[key]
    # Find actual files
    for line in ld_confoutput:
        stripped_line = line.strip()
        for key in prequisites:
            if key in known_resolutions:
                prequisites[key] = known_resolutions[key]
            elif stripped_line.startswith(key):
                prequisites[key] = stripped_line.split('=>')[1].strip()
    for key in prequisites:
        path = prequisites[key]
        realpath = os.path.realpath(prequisites[key])
        if i == 1:
            if path == realpath:
                print('  {0} -> {1}'.format(key, path))
            else:
                print('  {0} -> {1} ({2})'.format(key, path, realpath))
            if (MKL_ROOT is None) and key.startswith('libmkl'):
                MKL_ROOT = os.path.dirname(realpath)
        filename = os.path.basename(realpath)
        symname = os.path.basename(key)
        if filename == symname:
            shutil.copy(path, key)
        else:
            shutil.copy(realpath, filename)
            try:
                os.remove(symname)
            except FileNotFoundError:
                pass
            os.symlink(filename, symname)

# Copy MKL computational kernels if linked with MKL
if MKL_ROOT is not None:
    print('Detected oneAPI MKL usage: {0}'.format(MKL_ROOT))
    mkl_comp_libs = ('libmkl_def', 'libmkl_avx2')
    mkl_libs = []
    for mkl_lib in os.listdir(MKL_ROOT):
        for comp_lib in mkl_comp_libs:
            if mkl_lib.startswith(comp_lib):
                mkl_libs.append(mkl_lib)
                break
    for mkl_lib in mkl_libs:
        path = os.path.join(MKL_ROOT,mkl_lib)
        realpath = os.path.realpath(path)
        filename = os.path.basename(realpath)
        symname = os.path.basename(mkl_lib)
        print('  Copying MKL kernel lib: {0}'.format(mkl_lib))
        if filename == symname:
            shutil.copy(path, mkl_lib)
        else:
            shutil.copy(realpath, filename)
            try:
                os.remove(symname)
            except FileNotFoundError:
                pass
            os.symlink(filename, symname)

# Patch files with $ORIGIN runpath
for filename in patch_files:
    run_command("patchelf --set-rpath '$ORIGIN' {0}".format(filename))
