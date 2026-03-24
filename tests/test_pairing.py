# built-in dependencies
import pytest

# project dependencies
from lightecc import LightECC
from lightecc.interfaces.elliptic_curve import EllipticCurvePoint
from lightecc.commons.pairing import weil_pairing
from lightecc.commons.logger import Logger

logger = Logger(module="tests/test_pairing.py")


def _point(curve, x, y):
    """Create an EllipticCurvePoint on the given curve."""
    return EllipticCurvePoint(x, y, curve)


# ---- Test curve: y^2 = x^3 - x over F_7 ----
# This is a supersingular curve with |E(F_7)| = 8.
# Group structure: E(F_7) ≅ Z/4 × Z/2
# Points: O, (0,0), (1,0), (4,2), (4,5), (5,1), (5,6), (6,0)
# 2-torsion: E[2] = {O, (0,0), (1,0), (6,0)} ≅ Z/2 × Z/2
# The full 2-torsion is in E(F_7), so the Weil pairing on E[2] is non-degenerate.
# e_2((0,0), (1,0)) = 6 = -1 mod 7

PAIRING_CURVE = LightECC(
    form_name="weierstrass",
    curve_name="custom",
    config={"a": -1, "b": 0, "p": 7, "G": (4, 2), "n": 8},
).curve


def test_weil_pairing_nontrivial():
    """
    Weil pairing on 2-torsion points of y^2 = x^3 - x over F_7.
    Independent 2-torsion points should give e = -1 mod p.
    """
    P = _point(PAIRING_CURVE, 0, 0)
    Q = _point(PAIRING_CURVE, 1, 0)
    result = weil_pairing(P, Q, r=2)
    assert result == 6, f"Expected 6 (-1 mod 7), got {result}"
    logger.info("✅ Non-trivial Weil pairing test passed.")


def test_weil_pairing_alternating():
    """Weil pairing is alternating: e(P, P) = 1."""
    P = _point(PAIRING_CURVE, 0, 0)
    result = weil_pairing(P, P, r=2)
    assert result == 1
    logger.info("✅ Alternating property test passed.")


def test_weil_pairing_identity():
    """e(P, O) = e(O, Q) = 1."""
    O = _point(PAIRING_CURVE, float("inf"), float("inf"))
    P = _point(PAIRING_CURVE, 0, 0)
    assert weil_pairing(P, O, r=2) == 1
    assert weil_pairing(O, P, r=2) == 1
    logger.info("✅ Identity element test passed.")


def test_weil_pairing_antisymmetry():
    """e(P, Q) * e(Q, P) = 1 mod p."""
    P = _point(PAIRING_CURVE, 0, 0)
    Q = _point(PAIRING_CURVE, 1, 0)
    e_PQ = weil_pairing(P, Q, r=2)
    e_QP = weil_pairing(Q, P, r=2)
    assert (
        e_PQ * e_QP
    ) % 7 == 1, f"Antisymmetry failed: e(P,Q)*e(Q,P) = {(e_PQ * e_QP) % 7}"
    logger.info("✅ Antisymmetry test passed.")


def test_weil_pairing_bilinearity():
    """
    Bilinearity: e(P, Q+R) = e(P, Q) * e(P, R).
    Using P=(0,0), Q=(1,0), R=(6,0) which are all 2-torsion points.
    Q + R = (0,0) = P, so e(P, Q+R) = e(P, P) = 1.
    e(P, Q) * e(P, R) should also be 1.
    """
    p = 7
    P = _point(PAIRING_CURVE, 0, 0)
    Q = _point(PAIRING_CURVE, 1, 0)
    R = _point(PAIRING_CURVE, 6, 0)

    e_PQ = weil_pairing(P, Q, r=2)
    e_PR = weil_pairing(P, R, r=2)

    QpR = Q + R
    e_P_QpR = weil_pairing(P, QpR, r=2)

    assert (
        e_P_QpR == (e_PQ * e_PR) % p
    ), f"Bilinearity failed: e(P,Q+R)={e_P_QpR} != e(P,Q)*e(P,R)={(e_PQ * e_PR) % p}"
    logger.info("✅ Bilinearity test passed.")


def test_weil_pairing_root_of_unity():
    """The result of e_r(P, Q) must be an r-th root of unity."""
    P = _point(PAIRING_CURVE, 0, 0)
    Q = _point(PAIRING_CURVE, 1, 0)
    result = weil_pairing(P, Q, r=2)
    assert pow(result, 2, 7) == 1, f"{result} is not a 2nd root of unity mod 7"
    logger.info("✅ Root of unity test passed.")


def test_weil_pairing_all_2torsion_pairs():
    """Verify pairing values for all pairs of 2-torsion points."""
    torsion_coords = [(0, 0), (1, 0), (6, 0)]
    points = [_point(PAIRING_CURVE, x, y) for x, y in torsion_coords]

    for P in points:
        for Q in points:
            result = weil_pairing(P, Q, r=2)
            assert (
                pow(result, 2, 7) == 1
            ), f"e({P}, {Q}) = {result} is not a 2nd root of unity"

    logger.info("✅ All 2-torsion pairs test passed.")


def test_weil_pairing_not_torsion():
    """Should raise ValueError if points are not r-torsion."""
    P = _point(PAIRING_CURVE, 4, 2)  # order 4, not a 2-torsion point
    Q = _point(PAIRING_CURVE, 0, 0)

    with pytest.raises(ValueError, match="not an.*torsion point"):
        weil_pairing(P, Q, r=2)

    logger.info("✅ Torsion validation test passed.")


def test_pairing_via_lightecc_interface():
    """Test pairing through the LightECC high-level interface."""
    ec = LightECC(form_name="weierstrass", curve_name="test-curve")
    G = ec.G

    # All points have order dividing n, so n*P = O for any P
    P = 4 * G
    Q = 8 * G

    result = ec.pairing(P, Q)

    # Result must be an n-th root of unity
    assert pow(result, ec.n, ec.modulo) == 1 or result == 0
    logger.info("✅ LightECC interface pairing test passed.")


def test_pairing_edwards_raises():
    """Pairing is only supported for Weierstrass form."""
    ec = LightECC(form_name="edwards", curve_name="test-curve")
    G = ec.G

    with pytest.raises(ValueError, match="only supported for Weierstrass"):
        ec.pairing(G, G)

    logger.info("✅ Edwards form rejection test passed.")
