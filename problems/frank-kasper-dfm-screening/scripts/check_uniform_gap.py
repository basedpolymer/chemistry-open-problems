#!/usr/bin/env python3
"""Exact constants for a uniform DFM gap above the spherical relaxation.

Python 3.10+, standard library only; runtime < 1 second.
The geometric proof is in paper.tex, section 'A uniform gap above 1'.
This checks its rational parameter conditions, not the proof in a proof
assistant. The stronger version invokes the established Kepler theorem;
it does not reverify Flyspeck. No numerical minimization enters the checks.
"""
from fractions import Fraction as Q
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def certificate(name, eta, lam, tau, target, packing_square_upper):
    """Prove the target lower bound using exact rational comparisons.

    rho is an upper bound for equal-ball packing density, with
    rho**2 <= packing_square_upper. Strict inequality in the last
    comparison ensures F > target, including when the packing bound
    itself is attained. All inputs must be exact Fractions.
    """
    if not all(isinstance(x, Q) for x in
               (eta, lam, tau, target, packing_square_upper)):
        raise TypeError('certificate inputs must be exact Fractions')
    if not (0 < eta < 1 and 0 < lam < 1 and 0 < tau < 1 and target > 1):
        raise ValueError('invalid parameter range')
    a, b = 1-eta, 1+eta
    # Lens volume of two radius-a balls at distance 2*a*(1-lam),
    # minus the maximum missing volume of two good cells, divided by 4*pi/3.
    overlap_margin = a**3*lam**2*(3-lam)/2 - 2*tau*b**3
    if overlap_margin <= 0:
        raise ValueError('missing-volume bound does not enforce separation')
    C = min(eta**2, Q(10, 9)*a*a*tau*tau)
    q = a*(1-lam)/b
    delta_target = target**3-1
    density_required = q**3*(1-delta_target/C)
    if density_required <= 0:
        raise ValueError('target allows all cells to be discarded')
    density_square_margin = density_required**2-packing_square_upper
    if density_square_margin <= 0:
        raise ValueError('packing comparison cannot prove the requested bound')
    return {
        'name': name, 'status': 'PASS',
        'eta': str(eta), 'lambda': str(lam), 'tau': str(tau),
        'strict_F_lower_bound': str(target),
        'packing_density_squared_upper': str(packing_square_upper),
        'overlap_margin_divided_by_4pi_over3': str(overlap_margin),
        'defect_coefficient_C': str(C), 'radius_ratio_q': str(q),
        'delta_if_F_equals_target': str(delta_target),
        'required_packing_density': str(density_required),
        'packing_density_squared_margin': str(density_square_margin),
    }


def main():
    rows = [
        certificate('elementary Blichfeldt bound, reproved in manuscript',
                    Q(3, 4000), Q(1, 32), Q(7, 10000),
                    Q('1.000000004'), Q(25, 32)),
        certificate('Kepler theorem and pi < 22/7',
                    Q(9, 2500), Q(7, 100), Q(7, 2000),
                    Q('1.00000025'), Q(242, 441)),
    ]
    for row in rows:
        print('PASS', row['name'])
        print('  F >', row['strict_F_lower_bound'],
              '; all overlap and packing margins strictly positive (exact)')

    # Guard against interpreting the small universal gap as a separator
    # near sigma. The same sufficient inequalities cannot certify 1.068.
    for args in [(Q(3, 4000), Q(1, 32), Q(7, 10000), Q(25, 32)),
                 (Q(9, 2500), Q(7, 100), Q(7, 2000), Q(242, 441))]:
        eta, lam, tau, rho2 = args
        try:
            certificate('overclaim', eta, lam, tau, Q('1.068'), rho2)
        except ValueError:
            pass
        else:
            raise AssertionError('unsupported near-sigma bound was accepted')
    print('PASS unsupported near-sigma targets rejected')

    report = {
        'status': 'PASS',
        'scope': ('every periodic partition in the manuscript class, uniformly '
                  'in topology, number of cells, positive cell volumes and box; '
                  'no exclusion of vanishing-volume sequences near sigma'),
        'method': 'exact rational checks of constants in the analytical proof',
        'external_theorem': 'Kepler sphere packing, DOI 10.1017/fmp.2017.1',
        'certificates': rows,
        'global_catalogue_enclosure': ['1.00000025', '1.068575909'],
        'topology_ordering': 'UNRESOLVED',
    }
    (ROOT/'certificates'/'uniform_gap_report.json').write_text(
        json.dumps(report, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
