# 3rd party dependencies
import pytest

# project dependencies
from lightecc import LightECC
from lightecc.commons.errors import InvalidCurveOrder, PointNotOnCurve
from lightecc.interfaces.elliptic_curve import ProjectivePoint
from lightecc.commons.logger import Logger

logger = Logger(module="tests/test_projective.py")

# secp256k1 is y^2 = x^3 + 7 with prime order, so it has no point of order 2
SECP256K1 = LightECC(form_name="weierstrass", curve_name="secp256k1")
P_PRIME = SECP256K1.modulo
B = 7

# composite modulus n^2 with n = 11 * 17, both 2 mod 3
N = 11 * 17
N_SQUARE = N * N


def find_affine_point(b: int, p: int) -> tuple:
    for x in range(1, p):
        rhs = (x**3 + b) % p
        for y in range(1, p):
            if (y * y) % p == rhs:
                return x, y
    raise ValueError("no point found")


def count_points(b: int, p: int) -> int:
    squares = {}
    for y in range(p):
        squares[(y * y) % p] = squares.get((y * y) % p, 0) + 1
    return 1 + sum(squares.get((x**3 + b) % p, 0) for x in range(p))


def small_odd_order_curve() -> tuple:
    # y^2 = x^3 + b over a small prime p = 1 mod 3 with odd number of points
    p = 1009
    for b in range(1, p):
        if count_points(b, p) % 2 == 1:
            return p, b, find_affine_point(b, p)
    raise ValueError("no odd order curve found")


def build_secp256k1_projective() -> LightECC:
    return LightECC(
        form_name="projective",
        config={"modulus": P_PRIME, "b": B, "G": SECP256K1.G.get_point()},
    )


def test_projective_form_on_prime_field():
    x, y = SECP256K1.G.get_point()
    ec = build_secp256k1_projective()
    G, O = ec.G, ec.O

    assert ec.n is None
    assert G.to_affine() == (x, y)
    assert O.is_identity() and not G.is_identity()
    assert str(O) == "\U0001d4aa"

    # identity
    assert G + O == G and O + G == G
    assert G - G == O
    assert -O == O
    assert O * 7 == O and 7 * O == O
    assert G * 0 == O

    # arithmetic
    assert G + G == 2 * G
    assert 123 * G + 456 * G == 579 * G
    assert (G * 5) - (G * 2) == G * 3
    assert G * -5 == -(G * 5)
    assert SECP256K1.n * G == O

    # same point with a different Z
    assert ProjectivePoint(2 * x, 2 * y, 2, ec.curve) == G

    # matches affine arithmetic of weierstrass form
    for k in [2, 3, 17, 2**128 + 1, SECP256K1.n - 1]:
        assert (k * G).to_affine() == (k * SECP256K1.G).get_point()

    logger.info("✅ Projective form on prime field tested")


def test_precomputed_powers_match_double_and_add():
    ec = build_secp256k1_projective()
    curve = ec.curve
    G = ec.G.get_point()
    powers = curve.precompute_powers(G, 20)

    for k in [0, 1, 2, 3, 1000, 2**19 + 7]:
        expected = ProjectivePoint(*curve.double_and_add(G, k), curve=curve)
        actual = ProjectivePoint(*curve.multiply_with_powers(powers, k), curve=curve)
        assert actual == expected

    with pytest.raises(ValueError):
        curve.multiply_with_powers(powers, 2**20)

    with pytest.raises(ValueError):
        curve.multiply_with_powers(powers, -1)

    logger.info("✅ Precomputed powers tested")


def test_projective_form_with_order():
    p, b, (x, y) = small_odd_order_curve()
    ec = LightECC(form_name="projective", config={"modulus": p, "b": b, "G": (x, y)})
    order = ec.find_order()
    assert order * ec.G == ec.O
    assert count_points(b, p) % order == 0

    ec = LightECC(
        form_name="projective",
        config={"modulus": p, "b": b, "G": (x, y), "n": order},
    )
    assert ec.n == order

    with pytest.raises(InvalidCurveOrder):
        LightECC(
            form_name="projective",
            config={"modulus": p, "b": b, "G": (x, y), "n": order + 1},
        )

    logger.info("✅ Projective form with order tested")


def test_projective_form_over_composite_modulus():
    # E(Z_{n^2}) contains (kn : 1 : 0), reducing to infinity mod n but distinct mod n^2
    # and behave as k under addition, similar to (1 + n)^k = 1 + kn in Paillier
    ec = LightECC(form_name="projective", config={"modulus": N_SQUARE, "b": B, "G": (0, 1, 0)})

    def message_point(k):
        return ProjectivePoint(k * N, 1, 0, ec.curve)

    M1 = message_point(1)
    assert M1 != ec.O
    assert not M1.is_identity()
    assert M1 * N == ec.O
    assert M1 * (N + 3) == message_point(3)

    for k1, k2 in [(3, 5), (100, 150), (N - 1, 1)]:
        assert message_point(k1) + message_point(k2) == message_point((k1 + k2) % N)

    # those points cannot be written as affine points
    with pytest.raises(ValueError):
        M1.to_affine()

    logger.info("✅ Projective form over composite modulus tested")


def test_exceptional_case_with_point_of_order_2_is_rejected():
    # y^2 = x^3 + 7 over p = 2 mod 3 has p + 1 points, so a point of order 2
    p = 10007
    x, y = find_affine_point(B, p)
    ec = LightECC(form_name="projective", config={"modulus": p, "b": B, "G": (x, y)})
    T = 2502 * ec.G
    assert T.y == 0, "expected a point of order 2"

    # P - Q = T has order 2, so complete formulas return (0 : 0 : 0)
    assert ec.curve.add_points((T + ec.G).get_point(), ec.G.get_point()) == (0, 0, 0)
    with pytest.raises(PointNotOnCurve):
        (T + ec.G) + ec.G
    with pytest.raises(PointNotOnCurve):
        ProjectivePoint(0, 0, 0, ec.curve)

    logger.info("✅ Exceptional case of complete formulas tested")


def test_projective_form_validation():
    x, y = SECP256K1.G.get_point()

    with pytest.raises(ValueError):
        LightECC(form_name="projective")

    with pytest.raises(ValueError):
        LightECC(form_name="projective", config={"modulus": P_PRIME, "b": B})

    with pytest.raises(ValueError):
        LightECC(
            form_name="projective",
            curve_name="secp256k1",
            config={"modulus": P_PRIME, "b": B, "G": (x, y)},
        )

    with pytest.raises(PointNotOnCurve):
        LightECC(form_name="projective", config={"modulus": P_PRIME, "b": B, "G": (x, y + 1)})

    # points of different curves
    ec1 = build_secp256k1_projective()
    ec2 = LightECC(form_name="projective", config={"modulus": N_SQUARE, "b": B, "G": (0, 1, 0)})
    with pytest.raises(ValueError):
        ec1.G + ec2.G

    logger.info("✅ Projective form validation tested")
