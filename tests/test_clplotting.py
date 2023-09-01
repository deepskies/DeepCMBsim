"""
tests clplotting.py
"""

import pytest
nmt = pytest.importorskip("pymaster") # skips tests if pymaster not installed
import numpy as np
from deepcmbsim.cl_plotting import flatmap


def test_flatmap_scalar():
    ells = int(1e4+1)
    t_only_dict = {'clTT': np.random.rand(ells)}
    pixels = 192
    degrees = 5
    out_shape = nmt.synfast_flat(pixels, pixels, degrees*np.pi/180, degrees*np.pi/180, [np.random.rand(ells)], [0]).shape
    assert flatmap(pixels, degrees, cl_dict=t_only_dict).flatmap('T').shape == out_shape


def test_flatmap_TEB():
    ells = int(1e4+1)
    teb_dict = {'clTT': np.random.rand(ells),
                'clEE': np.random.rand(ells),
                'clBB': np.random.rand(ells),
                'clTE': np.random.rand(ells),
                'clEB': np.zeros(ells),
                'clTB': np.zeros(ells)}
    pixels = 192
    degrees = 5
    out_shape = nmt.synfast_flat(pixels, pixels, degrees*np.pi/180, degrees*np.pi/180, np.random.rand(6, ells), [0, 0, 0]).shape
    assert flatmap(pixels, degrees, cl_dict=teb_dict).flatmap('TEB').shape == out_shape


def test_flatmap_TQU():
    ells = int(1e4+1)
    teb_dict = {'clTT': np.random.rand(ells),
                'clEE': np.random.rand(ells),
                'clBB': np.random.rand(ells),
                'clTE': np.random.rand(ells),
                'clEB': np.zeros(ells),
                'clTB': np.zeros(ells)}
    pixels = 192
    degrees = 5
    out_shape = nmt.synfast_flat(pixels, pixels, degrees*np.pi/180, degrees*np.pi/180, np.random.rand(6, ells), [0, 2]).shape
    assert flatmap(pixels, degrees, cl_dict=teb_dict).flatmap('TQU').shape == out_shape
