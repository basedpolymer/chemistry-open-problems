#!/usr/bin/env python3
"""Analytic cube/shear checks of tensor integration and affine certificates.

Python 3.12+, numpy/scipy/python-flint. Runtime < 1 second after imports.
These cases test off-diagonal moments, periodic lifts, scale and translation,
and compare the affine bound with the independently known cube minimum.
"""
from fractions import Fraction as Q
import unittest
import numpy as np
from check_affine_orbits import tensor_data, verify_data


def cube(transform=None, offset=None):
    vertices = np.array([[i,j,k] for i in (0,2) for j in (0,2) for k in (0,2)],dtype=np.int64)
    index = {tuple(v):i for i,v in enumerate(vertices)}
    faces = [((2,0,0),(2,2,0),(2,2,2),(2,0,2)),
             ((0,2,0),(0,2,2),(2,2,2),(2,2,0)),
             ((0,0,2),(2,0,2),(2,2,2),(0,2,2))]
    tri = [[index[f[0]],index[f[1]],index[f[2]]] for f in faces]
    tri += [[index[f[0]],index[f[2]],index[f[3]]] for f in faces]
    mat = np.eye(3,dtype=np.int64) if transform is None else np.array(transform,dtype=np.int64)
    off = np.zeros(3,dtype=np.int64) if offset is None else np.array(offset,dtype=np.int64)
    return {'L':2*mat,'X':vertices@mat.T+off,'C':np.array([[1,1,1]])@mat.T+off,
            'ncell':np.array(1),'K':np.array(1),'tri':np.array(tri),
            'tsh':np.zeros((6,3,3),dtype=np.int64),'tface':np.array([0,1,2,0,1,2]),
            'fi':np.zeros(3,dtype=np.int64),'fj':np.zeros(3,dtype=np.int64),
            'fsj':np.eye(3,dtype=np.int64)}


class AffineTests(unittest.TestCase):
    def test_cube_moments_and_global_affine_minimum(self):
        ns,s,v,scale = tensor_data(cube())
        self.assertEqual(v,1)
        self.assertEqual(s,[[Q(int(i==j),12) for j in range(3)] for i in range(3)])
        metric = [[Q(int(i==j)) for j in range(3)] for i in range(3)]
        lower,upper,_ = verify_data(ns,s,v,scale,metric,Q(12)+Q(1,10**10))
        # For a unit cube, A=3, J=1/4 and F^3=5/3.
        self.assertLess(lower**3,Q(5,3))
        self.assertGreater(upper**3,Q(5,3))
        self.assertLess(upper-lower,Q(1,10**8))

    def test_shear_and_translation_tensor_law(self):
        m = np.array([[1,1,0],[0,1,1],[1,0,2]],dtype=np.int64)
        _,s,v,_ = tensor_data(cube(m,offset=[17,-9,5]))
        self.assertEqual(v,3)
        self.assertEqual(s,[[Q(int(x),4) for x in row] for row in m@m.T])
        _,s2,v2,_ = tensor_data(cube(m))
        self.assertEqual((s,v),(s2,v2))

    def test_invalid_metric_and_false_dual_bound_rejected(self):
        ns,s,v,scale = tensor_data(cube())
        metric = [[Q(int(i==j)) for j in range(3)] for i in range(3)]
        with self.assertRaises(ValueError):
            verify_data(ns,s,v,scale,metric,Q(11))
        metric[2][2] = Q(-1)
        with self.assertRaises(ValueError):
            verify_data(ns,s,v,scale,metric,Q(13))


if __name__ == '__main__':
    unittest.main(verbosity=2)
