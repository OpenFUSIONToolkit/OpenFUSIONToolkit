#!/usr/bin/env python3
"""Cross-language session-IO helper for the TokaMaker.jl test suite.

Drives the Python ``OpenFUSIONToolkit.TokaMaker`` interface so the Julia test
can verify that native TokaMaker session files (``save_tokamaker`` /
``load_tokamaker``) round-trip between the two language bindings (both wrap the
same Fortran routines in liboftpy).

Usage:
    xlang_session.py save <mesh.h5> <out_session.h5> <out_psi.txt>
        Solve the shared spheromak case, save the session, and dump psi.
    xlang_session.py load <mesh.h5> <in_session.h5> <out_psi.txt>
        Set up the same mesh, load the session from file, and dump psi.

``psi`` is written one raw (un-normalized) value per line. The mesh is loaded
from the shared fixture so both languages solve/load on an identical mesh.
"""
import sys

import numpy as np

from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.TokaMaker import TokaMaker
from OpenFUSIONToolkit.TokaMaker.meshing import load_gs_mesh

FE_ORDER = 2


def _build(mesh_file):
    pts, lc, reg, coil_dict, cond_dict = load_gs_mesh(mesh_file)
    env = OFT_env(nthreads=1)
    gs = TokaMaker(env)
    gs.setup_mesh(pts, lc, reg)
    gs.setup_regions(cond_dict=cond_dict, coil_dict=coil_dict)
    gs.settings.free_boundary = False
    gs.setup(order=FE_ORDER)
    return gs


def _solve_spheromak(gs):
    gs.set_profiles(ffp_prof={'type': 'linterp', 'x': [0.0, 1.0], 'y': [1.0, 0.0]},
                    pp_prof={'type': 'flat'})
    gs.p_scale = 0.0
    gs.settings.nl_tol = 1.0e-12
    gs.settings.maxits = 100
    gs.settings.urf = 0.0
    gs.update_settings()
    gs.init_psi()
    gs.solve()


def main():
    mode, mesh_file, session_file, psi_file = sys.argv[1:5]
    gs = _build(mesh_file)
    if mode == 'save':
        _solve_spheromak(gs)
        gs._tMaker_equil.save_TokaMaker(session_file)
    elif mode == 'load':
        gs.replace_eq(source_file=session_file)
    else:
        raise SystemExit("mode must be 'save' or 'load'")
    psi = gs.get_psi(False)
    np.savetxt(psi_file, np.asarray(psi))


if __name__ == '__main__':
    main()
