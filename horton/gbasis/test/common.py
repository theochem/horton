import json

from os import path

import numpy as np

from .. import gobasis

import mol_data as mdata
import gobasis_data as gdata


def _compose_fn(subpath, fn, ext=".npy"):
    cur_pth = path.split(__file__)[0]
    pth = cur_pth + "/cached/{}/{}{}".format(fn, subpath, ext)
    return np.load(pth).astype(np.float64)


def load_json(fn):
    return _compose_fn("er", fn)
    # cur_pth = path.split(__file__)[0]
    # pth = cur_pth + "/cached/json/{}".format(fn)
    # with open(pth) as fh:
    #     a = np.array(json.load(fh))
    # return a


def load_quad(fn):
    return _compose_fn("quads", fn)


def load_dipole(fn):
    return _compose_fn("dipoles", fn)


def load_dm(fn):
    return _compose_fn("dm", fn)


def load_olp(fn):
    return _compose_fn("olp", fn)


def load_na(fn):
    return _compose_fn("na", fn)


def load_kin(fn):
    return _compose_fn("kin", fn)


def load_er(fn):
    return _compose_fn("er", fn)


def load_orbsa_coeffs(fn):
    return _compose_fn("orbs_a_coeffs", fn)


def load_orbsa_occs(fn):
    return _compose_fn("orbs_a_occs", fn)


def load_orbsa_dms(fn):
    return _compose_fn("orbs_a_dms", fn)


def load_orbsb_coeffs(fn):
    return _compose_fn("orbs_b_coeffs", fn)


def load_orbsb_occs(fn):
    return _compose_fn("orbs_b_occs", fn)


def load_orbsb_dms(fn):
    return _compose_fn("orbs_b_dms", fn)


def load_obasis(fn):
    return gobasis.GOBasis(*getattr(gdata, fn))


def load_mdata(fn):
    return getattr(mdata, fn)
