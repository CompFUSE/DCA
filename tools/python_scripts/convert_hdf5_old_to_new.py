#!/bin/bash/python3

# Partially converts in place an hdf5 file part of the test suite, or generated from main_dca, to
# the new hdf5 interface, so that it can be used for testing, or for the BSE solver. The changes are:
# - Vectors of vectors and DCA functions are stored as a single dataset.
# - Strings are human readable in hdfview.
# - Complex numbers are stored as such
#
# Usage: python convert_hdf5_old_to_new.py <file>
#
# Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)

import h5py
import numpy as np
from sys import argv
import copy
import os
import shutil


assert(len(argv) == 2)
filename = argv[1]

print("Converting file ", filename)

inp = h5py.File(filename, 'r')
outname = filename + ".hdf5tmp"
out = h5py.File(outname, 'w')

def isOldFunc(node : h5py.Group) :
    return ("name" in node) and ("domain-sizes" in node) and ("data" in node)

def isNewFunc(node : h5py.Dataset) :
    return "domain-sizes" in node.attrs

def isOldVecOfVec(node : h5py.Group) :
    return ("equal-size" in node) and ("size" in node) and ("data" in node)


def toString(set) :
    str = ""
    try :
        str = "".join([chr(a) for a in set[()]])
        return str
    except : return set


def iterator(name, node):
    path, dname = ["/", name]
    try : path, dname = name.rsplit("/", 1)
    except: pass
    if path == "" : path = "/"

    if isinstance(node, h5py.Dataset):  # node is a dataset

        if isNewFunc(node) :
            data = node[...]
            if len(node.attrs["domain-sizes"]) == len(node.shape) - 1 : # is complex
                data = node[..., 0] + 1j * node[..., 1]
            out.create_dataset(name, data = data)
            for key in node.attrs : out[name].attrs[key] = node.attrs[key]

        elif dname == "type" and path.find("four-point") >= 0 :
            out[path].create_dataset("channels", data = [toString(node).encode()])

        elif dname == "cluster_greens_function_G_k_w" and len(node.shape) == 2 :
            data = node[0, :] + 1j * node[1, :]
            out.create_dataset(name, data = data)

        else:
            inp.copy(node, out[path], name = dname)

    else:  # node is a group
        if dname == "channels" :
            values = []
            for a in inp[name]["data"] :
                val = (''.join(chr(i) for i in a)).encode()
                values.append(val)
            out.create_dataset(name, data = values)
            return

        elif isOldFunc(node) :
            data = node["data"][...]
            if len(node["domain-sizes"][:]) ==  len(node["data"].shape) - 1 : # is complex
                data = node["data"][..., 0] + 1j * node["data"][..., 1]

            out.create_dataset(name, data = data)
            out[name].attrs["name"] = toString(inp[name]["name"])
            out[name].attrs["domain-sizes"] = inp[name]["domain-sizes"]
            return

        elif isOldVecOfVec(node) :
            if node["equal-size"][0] : data = [line for line in node["data"][:]]
            else : data = [node["data"][key][:] for key in node["data"]]
            out.create_dataset(name, data = data)

        else :
            if name != "" : out.create_group(name)
            for key in node : iterator(name + "/" + key, node[key])

iterator("", inp)

inp.close()
out.close()

shutil.move(outname, filename)
