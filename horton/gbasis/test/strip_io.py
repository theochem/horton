#!/usr/bin/env python

from horton import *
import numpy as np
from glob import glob

log.set_level(log.silent)
np.set_printoptions(precision=16)


def grep_file():
    """Searches a test file for IOData and saves the data"""
    filenames = []
    skipped = []
    for filename in glob("test_*.py"):
        print("Processing: ", filename)
        with open(filename) as fh:
            for i in fh:
                if ".json" in i:
                    skipped.append((filename, i))
                    continue
                if "test/" in i:
                    fn, san_fn = santize_fn(i)
                    filenames.append((fn, san_fn))

    filenames = set(filenames)
    skipped = set(skipped)

    for fn, san_fn in filenames:
        if "fchk" in san_fn:
            flog = fn[:-4] + "log"
            try:
                mol = IOData.from_file(context.get_fn(fn), context.get_fn(flog))
                save_ints(mol, san_fn)
            except (IOError):
                pass

    for i in (save_gobasis_params, save_dms, save_exps, save_moldata, save_quads, save_dipoles):
        for fn, san_fn in filenames:
            mol = IOData.from_file(context.get_fn(fn))
            try:
                i(mol, san_fn)
            except AttributeError:
                pass
        print "-" * 80

    print("Processed:")
    for f in filenames:
        print(f)
    print("Skipped:")
    for s in skipped:
        print(s)


def santize_fn(line):
    if "\'" in line:
        fn = line.split("'")[1]
    else:
        fn = line.split('"')[1]
    san_fn = fn.replace(".", "_").replace("-", "_").split("/")[1]
    return fn, san_fn


def save_ints(mol, san_fn):
    for i in ("olp", "kin", "na", "er", "two_mo"):
        try:
            np.save("ints/{}/{}".format(i, san_fn), getattr(mol, i))
        except AttributeError:
            pass


def save_gobasis_params(mol, san_fn):
    obasis = mol.obasis
    print san_fn + " = ",
    print((obasis.centers, obasis.shell_map, obasis.nprims, obasis.shell_types, obasis.alphas, obasis.con_coeffs))


def save_dms(mol, san_fn):
    dm = mol.get_dm_full()
    np.save("dms/" + san_fn, dm)


def save_exps(mol, san_fn):
    orba = mol.orb_alpha
    np.save("orbs_a/coeffs/" + san_fn, orba.coeffs)
    np.save("orbs_a/occs/" + san_fn, orba.occupations)
    np.save("orbs_a/dms/" + san_fn, orba.to_dm())
    orbb = mol.orb_beta
    np.save("orbs_b/coeffs/" + san_fn, orbb.coeffs)
    np.save("orbs_b/occs/" + san_fn, orbb.occupations)
    np.save("orbs_b/dms/" + san_fn, orbb.to_dm())


def save_moldata(mol, san_fn):
    d = {"coordinates": mol.coordinates,
         "numbers": mol.numbers,
         "pseudo_numbers": mol.pseudo_numbers}
    print san_fn + "=",
    print(d)


def save_quads(mol, san_fn):
    quad = mol.quadrupole_moment
    np.save("quads/" + san_fn, quad)


def save_dipoles(mol, san_fn):
    dipole = mol.dipole_moment
    np.save("dipoles/" + san_fn, dipole)


grep_file()
