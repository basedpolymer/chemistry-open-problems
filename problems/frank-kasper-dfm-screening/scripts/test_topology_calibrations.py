"""Analytic cubic torus and deliberately false certificates; no optimizer."""
import copy
import itertools
import unittest
import numpy as np
from check_topology_calibrations import verify_field


class CalibrationTests(unittest.TestCase):
    def setUp(self):
        self.points = np.zeros((1, 3), dtype=np.int64)
        self.tets = []
        for order in itertools.permutations(range(3)):
            v = [0, 0, 0, 0]
            tet = [v.copy()]
            for axis in order:
                v[axis+1] += 1
                tet.append(v.copy())
            self.tets.append(tet)
        self.tets = np.array(self.tets, dtype=np.int64)
        # The six tetrahedra triangulate the cubic torus. W=I/2 satisfies
        # all edge constraints, including the body diagonal sqrt(3)/2.
        self.cert = dict(id='cube', W=np.eye(3, dtype=int).tolist(),
                         Z=[[0, 0, 0]], denominator=2, F_lower='0.7')

    def test_cube_and_gauge_translation(self):
        a = verify_field(self.cert, self.points, self.tets)
        shifted = copy.deepcopy(self.cert)
        shifted['Z'] = [[-127, 83, 17]]
        b = verify_field(shifted, self.points, self.tets)
        self.assertEqual(a, b)
        self.assertEqual(a['det_W'], '1/8')

    def test_interface_constraint_rejected(self):
        self.cert['denominator'] = 1
        with self.assertRaisesRegex(ValueError, 'norm constraint'):
            verify_field(self.cert, self.points, self.tets)

    def test_overclaim_rejected(self):
        self.cert['F_lower'] = '1.1'
        with self.assertRaisesRegex(ValueError, 'lower bound'):
            verify_field(self.cert, self.points, self.tets)

    def test_missing_tetrahedron_rejected(self):
        with self.assertRaisesRegex(ValueError, 'volume identity'):
            verify_field(self.cert, self.points, self.tets[:-1])


if __name__ == '__main__':
    unittest.main()
