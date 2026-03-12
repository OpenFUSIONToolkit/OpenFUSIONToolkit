#!/bin/bash

PATH=/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/oft_venv/bin:/home/runner/work/_temp/ghcca-node/node/bin:/snap/bin:/home/runner/.local/bin:/opt/pipx_bin:/home/runner/.cargo/bin:/home/runner/.config/composer/vendor/bin:/usr/local/.ghcup/bin:/home/runner/.dotnet/tools:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin OFT_HAVE_MPI=0 OFT_DEBUG_TEST=0 /home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/oft_venv/bin/pytest $*
