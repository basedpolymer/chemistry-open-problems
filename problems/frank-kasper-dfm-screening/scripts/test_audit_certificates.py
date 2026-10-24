#!/usr/bin/env python3
"""Adversarial and invariance checks for the certificate checker (~10 s)."""
from copy import deepcopy
from pathlib import Path
import unittest
import numpy as np
from audit_certificates import canon, exact_energy, intended_box, surface, Q
from certify_catalogue import verify

ROOT = Path(__file__).resolve().parents[1]


class AuditTests(unittest.TestCase):
    def test_sphere_and_disk(self):
        sphere = [(1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)]
        surface(sphere)
        surface(sphere[1:], disk=True)
        with self.assertRaises(ValueError):
            surface(sphere[1:])
        with self.assertRaises(ValueError):
            surface(sphere+[sphere[0]])

    def test_pinched_link_is_rejected(self):
        first = [(1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)]
        second = [tuple(0 if v == 0 else v+3 for v in tri) for tri in first]
        with self.assertRaisesRegex(ValueError, 'vertex link'):
            surface(first+second)

    def test_periodic_canonicalization(self):
        v = [(0, 0, 0, 0), (1, 0, -1, 0), (0, 1, 0, 0)]
        shifted = [(i, x+7, y-3, z+1) for i, x, y, z in v]
        self.assertEqual(canon(v), canon(reversed(shifted)))

    def test_energy_invariance_and_orientation(self):
        with np.load(ROOT/'certificates'/'meshes'/'BCC_fixedbox_L4.npz') as z:
            c = {k: z[k] for k in z.files}
        lo, hi, _ = exact_energy(c)
        transformed = deepcopy(c)
        for k in ('X', 'L', 'C'):
            transformed[k] = transformed[k]*2
        transformed['X'] += np.array([19, -7, 13], dtype=np.int64)
        transformed['C'] += np.array([19, -7, 13], dtype=np.int64)
        lo2, hi2, _ = exact_energy(transformed)
        self.assertEqual((lo, hi), (lo2, hi2))
        bad = deepcopy(c)
        bad['tri'][0, [0, 1]] = bad['tri'][0, [1, 0]]
        bad['tsh'][0, [0, 1]] = bad['tsh'][0, [1, 0]]
        with self.assertRaisesRegex(ValueError, 'non-positive cone'):
            exact_energy(bad)

    def test_hexagonal_interval_transport(self):
        with np.load(ROOT/'certificates'/'meshes'/'mu_fixedbox_L2.npz') as z:
            c = {k: z[k] for k in z.files}
        bound, metadata = intended_box(c, 'mu', Q(1))
        self.assertLess(bound, Q(100000000001, 100000000000))
        bad = deepcopy(c)
        bad['L'][0, 0] += 10**9
        with self.assertRaises(ValueError):
            intended_box(bad, 'mu', Q(1))

    def test_catalogue_rejects_missing_tetrahedron(self):
        with np.load(ROOT/'certificates'/'catalogue_topology.npz') as z:
            P, T = z['I_Phys26_points'], z['I_Phys26_tets']
        verify(P, T)
        with self.assertRaises(ValueError):
            verify(P, T[:-1])


if __name__ == '__main__':
    unittest.main(verbosity=2)
