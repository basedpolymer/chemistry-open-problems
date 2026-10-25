"""Check the two-sided log barrier derivative and rejection at the boundary."""
import unittest
import numpy as np
from run_barrier import cone_barrier
from check_barrier_witnesses import centroid_check
from test_affine_orbits import cube


class BarrierTests(unittest.TestCase):
    def setUp(self):
        self.X = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
        self.C = np.array([[.2, .2, -1.], [.2, .2, 1.]])
        self.args = (np.eye(3), np.array([[0, 1, 2]]),
                     np.zeros((1, 3, 3), dtype=int), np.array([0]),
                     np.array([0]), np.array([1]), np.zeros((1, 3), dtype=int))

    def test_directional_gradient_and_translation(self):
        rng = np.random.default_rng(17)
        dx, dc = rng.normal(size=self.X.shape), rng.normal(size=self.C.shape)
        value, gx, gc, mind = cone_barrier(self.X, self.C, *self.args, .2)
        h = 1e-6
        plus = cone_barrier(self.X+h*dx, self.C+h*dc, *self.args, .2)[0]
        minus = cone_barrier(self.X-h*dx, self.C-h*dc, *self.args, .2)[0]
        self.assertAlmostEqual((plus-minus)/(2*h), np.sum(gx*dx)+np.sum(gc*dc), places=7)
        np.testing.assert_allclose(gx.sum(axis=0)+gc.sum(axis=0), 0, atol=1e-14)
        self.assertAlmostEqual(value, -np.log(.8), places=14)
        self.assertEqual(mind, 1.)

    def test_boundary_rejected(self):
        self.C[0, 2] = 0
        self.assertTrue(np.isinf(cone_barrier(self.X, self.C, *self.args)[0]))

    def test_true_centroid_independent_of_auxiliary_apex(self):
        c = cube()
        expected = centroid_check(c)
        for key in ('X', 'C', 'L'):
            c[key] *= 4
        c['C'][0] = [1, 2, 3]  # Interior, away from the true centroid [4,4,4].
        self.assertEqual(centroid_check(c), expected)
        self.assertEqual(expected['negative_centroid_cones'], 0)
        self.assertEqual(expected['minimum_centroid_cone_over_mean_cell_volume'], '1/12')


if __name__ == '__main__':
    unittest.main()
