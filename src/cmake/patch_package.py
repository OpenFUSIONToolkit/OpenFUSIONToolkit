import os
import shutil
import subprocess

def run_command(command):
    pid = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outs, errs = pid.communicate()
    errcode = pid.poll()
    return outs.decode(), errs.decode(), errcode

def get_prereqs(filename):
    output, error, errcode = run_command('readelf -d {0}'.format(filename))
    # print(filename)
    preqs = []
    for line in output.splitlines():
        if line.find('Shared library:') > 0:
            preqs.append(line.split('[')[1].split(']')[0])
    # for preq in preqs:
    #     print('  '+preq)
    return preqs

output, error, errcode = run_command('ldconfig -p')
ld_confoutput = output.splitlines()

for i in range(2):
    full_prequisites = {}
    # Gather prerequisites
    patch_files = [f for f in os.listdir() if os.path.isfile(f)]
    for filename in patch_files:
        preqs = get_prereqs(filename)
        for preq in preqs:
            if preq not in full_prequisites:
                full_prequisites[preq] = None

    # Find system libraries that are not typically in a base install
    filter_list = ['libgomp.','libgfortran.','libquadmath.']
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
            if stripped_line.startswith(key):
                prequisites[key] = stripped_line.split('=>')[1].strip()
    for key in prequisites:
        path = prequisites[key]
        realpath = os.path.realpath(prequisites[key])
        if i == 1:
            print('{0} -> {1} ({2})'.format(key, path, realpath))
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

# Patch files with $ORIGIN runpath
for filename in patch_files:
    run_command("patchelf --set-rpath '$ORIGIN' {0}".format(filename))
