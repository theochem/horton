#!/usr/bin/env python
import json

from os import mkdir

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
                # if ".json" in i:
                #     skipped.append((filename, i))
                #     continue
                if "test/" in i:
                    fn, san_fn = santize_fn(i)
                    filenames.append((fn, "cached/" + san_fn))

    filenames = set(filenames)
    filenames2 = set()
    skipped = set(skipped)

    for fn, san_fn in filenames:
        try:
            mkdir(san_fn)
        except OSError:
            pass

        if "json" in san_fn:
            save_json(context.get_fn(fn), san_fn)
            continue

        if "fchk" in san_fn:
            flog = fn[:-4] + "log"
            try:
                mol = IOData.from_file(context.get_fn(fn), context.get_fn(flog))
                save_ints(mol, san_fn)
            except IOError:
                pass

        filenames2.add((fn, san_fn))

    for i in (save_gobasis_params, save_dms, save_exps, save_moldata, save_quads, save_dipoles):
        for fn, san_fn in filenames2:
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
            np.save("{}/{}".format(san_fn, i), getattr(mol, i).astype(np.float32))
        except AttributeError:
            pass


def save_gobasis_params(mol, san_fn):
    obasis = mol.obasis
    print san_fn + " = ",
    print(obasis.centers, obasis.shell_map, obasis.nprims, obasis.shell_types, obasis.alphas, obasis.con_coeffs)


def save_dms(mol, san_fn):
    dm = mol.get_dm_full().astype(np.float32)
    np.save(san_fn + "/dm", dm)


def save_exps(mol, san_fn):
    orba = mol.orb_alpha
    np.save(san_fn + "/orbs_a_coeffs", orba.coeffs)
    np.save(san_fn + "/orbs_a_occs", orba.occupations)
    np.save(san_fn + "/orbs_a_dms", orba.to_dm())
    orbb = mol.orb_beta
    np.save(san_fn + "/orbs_b_coeffs", orbb.coeffs)
    np.save(san_fn + "/orbs_b_occs", orbb.occupations)
    np.save(san_fn + "/orbs_b_dms", orbb.to_dm())


def save_moldata(mol, san_fn):
    d = {"coordinates": mol.coordinates,
         "numbers": mol.numbers,
         "pseudo_numbers": mol.pseudo_numbers}
    print san_fn + "=",
    print(d)


def save_quads(mol, san_fn):
    quad = mol.quadrupole_moment.astype(np.float32)
    np.save(san_fn + "/quads", quad)


def save_dipoles(mol, san_fn):
    dipole = mol.dipole_moment.astype(np.float32)
    np.save(san_fn + "/dipoles", dipole)


def save_json(json_fn, san_fn):
    with open(json_fn) as fh:
        arr = np.array(json.load(fh))
    np.save(san_fn + "/er", arr)


grep_file()
