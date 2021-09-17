import itertools
import logging

import numpy as np
import pytest
from mpmath import mp
from qecsim import paulitools as pt
from qecsim.models.generic import BiasedDepolarizingErrorModel, DepolarizingErrorModel

from qsdcolorxyz import Color666XYZCode, Color666XYZMPSDecoder


def _is_close(a, b, rtol=1e-05, atol=1e-08):
    # np.isclose for mp.mpf, i.e. absolute(a - b) <= (atol + rtol * absolute(b))
    try:
        return [mp.almosteq(le, ri, rel_eps=rtol, abs_eps=atol) for le, ri in itertools.zip_longest(a, b)]
    except TypeError:
        return mp.almosteq(a, b, rel_eps=rtol, abs_eps=atol)


def test_color666xyz_mps_decoder_properties():
    decoder = Color666XYZMPSDecoder(chi=8, tol=1e-8)
    assert isinstance(decoder.label, str)
    assert isinstance(repr(decoder), str)
    assert isinstance(str(decoder), str)


@pytest.mark.parametrize('chi, tol', [
    (None, None),
    (6, None),
    (None, 0.1),
    (None, 1),
])
def test_color666xyz_mps_decoder_new_valid_parameters(chi, tol):
    Color666XYZMPSDecoder(chi=chi, tol=tol)  # no error raised


@pytest.mark.parametrize('chi, tol', [
    (-1, None),  # invalid chi
    (0.1, None),  # invalid chi
    ('asdf', None),  # invalid chi
    (None, -1),  # invalid tol
    (None, 'asdf'),  # invalid tol
])
def test_color666xyz_mps_decoder_new_invalid_parameters(chi, tol):
    with pytest.raises((ValueError, TypeError), match=r"^Color666XYZMPSDecoder") as exc_info:
        Color666XYZMPSDecoder(chi=chi, tol=tol)
    print(exc_info)


@pytest.mark.parametrize('error_pauli', [
    Color666XYZCode(3).new_pauli().site('X', (2, 1)).site('Y', (3, 3)),
    Color666XYZCode(5).new_pauli().site('X', (3, 1)).site('Y', (2, 2)).site('Z', (6, 4)),
    Color666XYZCode(7).new_pauli().site('X', (4, 2)).site('Y', (4, 3)).site('Z', (8, 4), (8, 5)),
])
def test_color666xyz_mps_decoder_sample_recovery(error_pauli):
    print('\nerror:\n{}'.format(error_pauli))
    error = error_pauli.to_bsf()
    code = error_pauli.code
    syndrome = pt.bsp(error, code.stabilizers.T)
    print('\nsyndrome:\n{}'.format(code.ascii_art(syndrome=syndrome)))
    recovery_pauli = Color666XYZMPSDecoder.sample_recovery(code, syndrome)
    print('\nrecovery:\n{}'.format(recovery_pauli))
    recovery = recovery_pauli.to_bsf()
    assert np.array_equal(pt.bsp(recovery, code.stabilizers.T), syndrome), (
        'recovery {} does not give the same syndrome as the error {}'.format(recovery, error))
    assert np.all(pt.bsp(recovery ^ error, code.stabilizers.T) == 0), (
        'recovery ^ error ({} ^ {}) does not commute with stabilizers.'.format(recovery, error))


@pytest.mark.parametrize('chi, tol, rtol', [
    # varying chi
    (None, None, 1e-14),
    (2, None, 1e-14),
    (4, None, 1e-14),
    (8, None, 1e-14),
    (16, None, 1e-14),
    (32, None, 1e-14),
    # with tol
    (None, 1e-8, 1e-14),
    (2, 1e-8, 1e-14),
    (4, 1e-8, 1e-14),
    (8, 1e-8, 1e-14),
    (16, 1e-8, 1e-14),
    (32, 1e-8, 1e-14),
])
def test_color666xyz_mps_decoder_cosets_probability_inequality(chi, tol, rtol):
    code = Color666XYZCode(5)
    decoder = Color666XYZMPSDecoder(chi=chi, tol=tol)
    # probabilities
    prob_dist = DepolarizingErrorModel().probability_distribution(0.1)
    # coset probabilities for null Pauli
    coset_ps, _ = decoder._coset_probabilities(prob_dist, code.new_pauli())
    coset_i_p, coset_x_p, coset_y_p, coset_z_p = coset_ps
    # expect Pr(IG) > Pr(XG) ~= Pr(YG) ~= Pr(ZG)
    print()
    print('Pr(IG):{!r} > Pr(XG):{!r} ~= Pr(YG):{!r} ~= Pr(ZG):{!r}. rtol={}. rtol={}'.format(
        coset_i_p, coset_x_p, coset_y_p, coset_z_p,
        abs(coset_x_p - coset_y_p) / abs(coset_y_p),
        abs(coset_y_p - coset_z_p) / abs(coset_z_p)))
    print('types: Pr(IG):{}, Pr(XG):{}, Pr(YG):{}, Pr(ZG):{}'.format(
        type(coset_i_p), type(coset_x_p), type(coset_y_p), type(coset_z_p)))
    assert coset_i_p > coset_x_p, 'Coset probabilites do not satisfy Pr(IG) > Pr(XG)'
    assert coset_i_p > coset_y_p, 'Coset probabilites do not satisfy Pr(IG) > Pr(YG)'
    assert coset_i_p > coset_z_p, 'Coset probabilites do not satisfy Pr(IG) > Pr(ZG)'
    assert _is_close(coset_x_p, coset_y_p, rtol=rtol, atol=0), 'Coset probabilites do not satisfy Pr(XG) ~= Pr(YG)'
    assert _is_close(coset_y_p, coset_z_p, rtol=rtol, atol=0), 'Coset probabilites do not satisfy Pr(YG) ~= Pr(ZG)'


@pytest.mark.parametrize('chi, tol', [
    # varying chi
    # (Note: rtol is getting larger with larger chi)
    # (Note: approximate cosets values are getting closer to exact values with larger chi)
    (8, None),
    (16, None),
    (32, None),
    # with tol
    (8, 1e-8),
    (16, 1e-8),
    (32, 1e-8),
])
def test_color666xyz_mps_decoder_cosets_probability_inequality_biased_noise(chi, tol):
    code = Color666XYZCode(7)
    decoder = Color666XYZMPSDecoder(chi=chi, tol=tol)
    # probabilities
    prob_dist = BiasedDepolarizingErrorModel(bias=100).probability_distribution(0.1)
    # coset probabilities for null Pauli
    coset_ps, _ = decoder._coset_probabilities(prob_dist, code.new_pauli())
    coset_i_p, coset_x_p, coset_y_p, coset_z_p = coset_ps
    # expect Pr(IG) > Pr(YG) > Pr(XG) > Pr(ZG)
    # N.B. YG has 4 Y and 3 Z, XG has 3 Y and 4 X, ZG has 4 Z and 3 X
    print()
    print('Pr(IG):{!r} > Pr(YG):{!r} > Pr(XG):{!r} > Pr(ZG):{!r}'.format(
        coset_i_p, coset_y_p, coset_x_p, coset_z_p))
    print('types: Pr(IG):{}, Pr(YG):{}, Pr(XG):{}, Pr(ZG):{}'.format(
        type(coset_i_p), type(coset_y_p), type(coset_x_p), type(coset_z_p)))
    assert coset_i_p > coset_y_p, 'Coset probabilites do not satisfy Pr(IG) > Pr(YG)'
    assert coset_y_p > coset_x_p, 'Coset probabilites do not satisfy Pr(YG) > Pr(XG)'
    assert coset_x_p > coset_z_p, 'Coset probabilites do not satisfy Pr(XG) > Pr(ZG)'


@pytest.mark.parametrize('size, rtol', [
    # varying size
    # (Note: rtol is getting larger with larger size for fixed chi)
    (3, 1e-15),
    (5, 1e-10),
    (7, 1e-8),
])
def test_color666xyz_mps_decoder_cosets_probability_equality_biased_noise(size, rtol):
    # rtol = 1e-9
    code = Color666XYZCode(size)
    decoder = Color666XYZMPSDecoder(chi=16)
    # probabilities
    biasedx_prob_dist = BiasedDepolarizingErrorModel(bias=100, axis='X').probability_distribution(0.1)
    biasedy_prob_dist = BiasedDepolarizingErrorModel(bias=100, axis='Y').probability_distribution(0.1)
    biasedz_prob_dist = BiasedDepolarizingErrorModel(bias=100, axis='Z').probability_distribution(0.1)
    # coset probabilities for null Pauli
    biasedx_coset_ps, _ = decoder._coset_probabilities(biasedx_prob_dist, code.new_pauli())
    biasedy_coset_ps, _ = decoder._coset_probabilities(biasedy_prob_dist, code.new_pauli())
    biasedz_coset_ps, _ = decoder._coset_probabilities(biasedz_prob_dist, code.new_pauli())
    # expect Pr(IG)_x = Pr(IG)_y = Pr(IG)_z
    print('Pr(IG)_x:{!r} ~= Pr(IG)_y:{!r} ~= Pr(IG)_z:{!r}'.format(
        biasedx_coset_ps[0], biasedy_coset_ps[0], biasedz_coset_ps[0]
    ))
    assert (_is_close(biasedx_coset_ps[0], biasedy_coset_ps[0], rtol=rtol, atol=0))
    assert (_is_close(biasedy_coset_ps[0], biasedz_coset_ps[0], rtol=rtol, atol=0))
    # expect Pr(XG)_x = Pr(YG)_y = Pr(ZG)_z
    print('Pr(XG)_x:{!r} ~= Pr(YG)_y:{!r} ~= Pr(ZG)_z:{!r}'.format(
        biasedx_coset_ps[1], biasedy_coset_ps[2], biasedz_coset_ps[3]
    ))
    assert (_is_close(biasedx_coset_ps[1], biasedy_coset_ps[2], rtol=rtol, atol=0))
    assert (_is_close(biasedy_coset_ps[2], biasedz_coset_ps[3], rtol=rtol, atol=0))
    # expect Pr(YG)_x = Pr(ZG)_y = Pr(XG)_z
    print('Pr(YG)_x:{!r} ~= Pr(ZG)_y:{!r} ~= Pr(XG)_z:{!r}'.format(
        biasedx_coset_ps[2], biasedy_coset_ps[3], biasedz_coset_ps[1]
    ))
    assert (_is_close(biasedx_coset_ps[2], biasedy_coset_ps[3], rtol=rtol, atol=0))
    assert (_is_close(biasedy_coset_ps[3], biasedz_coset_ps[1], rtol=rtol, atol=0))
    # expect Pr(ZG)_x = Pr(XG)_y = Pr(YG)_z
    print('Pr(ZG)_x:{!r} ~= Pr(XG)_y:{!r} ~= Pr(YG)_z:{!r}'.format(
        biasedx_coset_ps[3], biasedy_coset_ps[1], biasedz_coset_ps[2]
    ))
    assert (_is_close(biasedx_coset_ps[3], biasedy_coset_ps[1], rtol=rtol, atol=0))
    assert (_is_close(biasedy_coset_ps[1], biasedz_coset_ps[2], rtol=rtol, atol=0))


def test_color666xyz_mps_decoder_cosets_probability_triplet_optimisation():
    code = Color666XYZCode(5)
    decoder = Color666XYZMPSDecoder()
    # probabilities
    prob_dist = BiasedDepolarizingErrorModel(bias=10).probability_distribution(probability=0.1)
    # coset probabilities for null Pauli
    coset_i_ps, _ = decoder._coset_probabilities(prob_dist, code.new_pauli())
    coset_x_ps, _ = decoder._coset_probabilities(prob_dist, code.new_pauli().logical_x())
    coset_y_ps, _ = decoder._coset_probabilities(prob_dist, code.new_pauli().logical_x().logical_z())
    coset_z_ps, _ = decoder._coset_probabilities(prob_dist, code.new_pauli().logical_z())
    # expect Pr(iIG) ~= Pr(xXG)
    assert _is_close(coset_i_ps[0], coset_x_ps[1], rtol=0, atol=0), (
        'Coset probabilites do not satisfy Pr(iIG) ~= Pr(xXG)')
    # expect Pr(iXG) ~= Pr(xIG)
    assert _is_close(coset_i_ps[1], coset_x_ps[0], rtol=0, atol=0), (
        'Coset probabilites do not satisfy Pr(iXG) ~= Pr(xIG)')
    # expect Pr(iIG) ~= Pr(yYG)
    assert _is_close(coset_i_ps[0], coset_y_ps[2], rtol=0, atol=0), (
        'Coset probabilites do not satisfy Pr(iIG) ~= Pr(yYG)')
    # expect Pr(iYG) ~= Pr(yIG)
    assert _is_close(coset_i_ps[2], coset_y_ps[0], rtol=0, atol=0), (
        'Coset probabilites do not satisfy Pr(iXG) ~= Pr(xIG)')
    # expect Pr(iIG) ~= Pr(zZG)
    assert _is_close(coset_i_ps[0], coset_z_ps[3], rtol=0, atol=0), (
        'Coset probabilites do not satisfy Pr(iIG) ~= Pr(zZG)')
    # expect Pr(iZG) ~= Pr(zIG)
    assert _is_close(coset_i_ps[3], coset_z_ps[0], rtol=0, atol=0), (
        'Coset probabilites do not satisfy Pr(iZG) ~= Pr(zIG)')


@pytest.mark.parametrize('sample_pauli_f, sample_pauli_g', [
    (Color666XYZCode(5).new_pauli(), Color666XYZCode(5).new_pauli()),
    (Color666XYZCode(5).new_pauli(), Color666XYZCode(5).new_pauli().plaquette('Z', (1, 1)).plaquette('X', (5, 3))),
    (Color666XYZCode(5).new_pauli().logical_x(),
     Color666XYZCode(5).new_pauli().logical_x().plaquette('X', (2, 0)).plaquette('X', (3, 2)).plaquette('X', (5, 3))),
    (Color666XYZCode(5).new_pauli().logical_z(),
     Color666XYZCode(5).new_pauli().logical_z().plaquette('Z', (2, 0)).plaquette('Z', (3, 2)).plaquette('Z', (5, 3))),
])
def test_color666xyz_mps_decoder_cosets_probability_equivalence(sample_pauli_f, sample_pauli_g):
    decoder = Color666XYZMPSDecoder()
    # probabilities
    prob_dist = BiasedDepolarizingErrorModel(bias=10).probability_distribution(probability=0.1)
    # coset probabilities
    coset_f_ps, _ = decoder._coset_probabilities(prob_dist, sample_pauli_f)
    coset_g_ps, _ = decoder._coset_probabilities(prob_dist, sample_pauli_g)
    print('#Pr(fG)=', coset_f_ps)
    print('#Pr(gG)=', coset_g_ps)
    assert all(_is_close(coset_f_ps, coset_g_ps, rtol=1e-14, atol=0)), (
        'Coset probabilites do not satisfy Pr(fG) ~= Pr(gG)')


@pytest.mark.parametrize('error_pauli, chi', [
    (Color666XYZCode(3).new_pauli().site('X', (2, 1)).site('Y', (3, 3)), None),
    (Color666XYZCode(3).new_pauli().site('X', (2, 1)).site('Y', (3, 3)), 8),
    (Color666XYZCode(5).new_pauli().site('X', (3, 1)).site('Y', (2, 2)).site('Z', (6, 4)), None),
    (Color666XYZCode(5).new_pauli().site('X', (3, 1)).site('Y', (2, 2)).site('Z', (6, 4)), 8),
    (Color666XYZCode(7).new_pauli().site('X', (4, 2)).site('Y', (4, 3)).site('Z', (8, 4), (8, 5)), 8),
])
def test_color666xyz_mps_decoder_decode(error_pauli, chi, caplog):
    with caplog.at_level(logging.WARN):
        error = error_pauli.to_bsf()
        code = error_pauli.code
        syndrome = pt.bsp(error, code.stabilizers.T)
        decoder = Color666XYZMPSDecoder(chi=chi)
        recovery = decoder.decode(code, syndrome)
        assert np.array_equal(pt.bsp(recovery, code.stabilizers.T), syndrome), (
            'recovery {} does not give the same syndrome as the error {}'.format(recovery, error))
        assert np.all(pt.bsp(recovery ^ error, code.stabilizers.T) == 0), (
            'recovery ^ error ({} ^ {}) does not commute with stabilizers.'.format(recovery, error))
        assert len(caplog.records) == 0, 'Unexpected log messages: {}'.format(caplog.text)


def test_color666xyz_mps_decoder_decode_small_codes_exact_approx():
    code = Color666XYZCode(5)
    exact_decoder = Color666XYZMPSDecoder()
    approx_decoder = Color666XYZMPSDecoder(chi=16)
    identity = code.new_pauli()
    # probabilities
    prob_dist = BiasedDepolarizingErrorModel(bias=10).probability_distribution(probability=0.1)
    # coset probabilities
    exact_coset_ps, _ = exact_decoder._coset_probabilities(prob_dist, identity)
    approx_coset_ps, _ = approx_decoder._coset_probabilities(prob_dist, identity)
    print('#exact Pr(G)=', exact_coset_ps)
    print('#approx Pr(G)=', approx_coset_ps)
    assert all(_is_close(exact_coset_ps, approx_coset_ps, rtol=1e-10, atol=0)), (
        'Coset probabilites do not satisfy exact Pr(G) ~= approx Pr(G)')


def test_color666xyz_mps_decoder_decode_value():
    # expected coset_ps (for regression)
    expected_coset_ps = (
        mp.mpf('1.0439958913009047e-11'),  # I
        mp.mpf('7.2763692299320466e-12'),  # X
        mp.mpf('1.9893345721051557e-9'),  # Y
        mp.mpf('2.6652981651379235e-12'),  # Z
    )
    code = Color666XYZCode(7)
    decoder = Color666XYZMPSDecoder(chi=16)
    # error
    error = code.new_pauli()
    error.site('X', (2, 1))
    error.site('Y', (7, 3))
    error.site('Z', (6, 0), (3, 3), (7, 6))
    # syndrome
    syndrome = pt.bsp(error.to_bsf(), code.stabilizers.T)
    # sample
    sample = decoder.sample_recovery(code, syndrome)
    print(sample)
    # probabilities
    prob_dist = (0.8, 0.04, 0.06, 0.1)
    # coset probabilities
    coset_ps, _ = decoder._coset_probabilities(prob_dist, sample)
    print('# expected Pr(G)=', expected_coset_ps)
    print('#   actual Pr(G)=', coset_ps)
    assert all(_is_close(expected_coset_ps, coset_ps, rtol=1e-10, atol=0)), (
        'Coset probabilites do not satisfy expected Pr(G) ~= Pr(G)')
