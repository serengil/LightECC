# project dependencies
from lightecc import LightECC
from lightecc.commons.logger import Logger

logger = Logger(module="tests/test_arithmetic.py")

FORMS = ["weierstrass", "koblitz", "edwards"]


def test_adding_a_point_with_its_negative():
    for form in FORMS:
        cs = LightECC(form_name=form)
        G = cs.curve.G
        G_minus = cs.curve.negative_point(G)
        assert cs.curve.O == cs.curve.add_points(G, G_minus)

        logger.info(f"✅ Adding a point with its negative tested for {form}")


def test_zero_times_base_point():
    for form in FORMS:
        cs = LightECC(form_name=form)
        G = cs.curve.G
        assert cs.curve.double_and_add(G, 0) == cs.curve.O
        logger.info(f"✅ Test 0 x G = O done for {form}")


def test_double_and_add_with_negative_input():
    for form in FORMS:
        cs = LightECC(form_name=form)
        G = cs.curve.G
        assert cs.curve.double_and_add(G, -10) == cs.curve.negative_point(
            cs.curve.double_and_add(G, 10)
        )
        logger.info(f"✅ Test (-10) x G = -(10 x G) done for {form}")


def test_elliptic_curve_cyclic_group_on_real_curves():
    for form in FORMS:
        ec = LightECC(form_name=form)

        # base point
        G = ec.G

        # order of the elliptic curve
        n = ec.n

        # neutral point
        nG = n * G

        _17G = 17 * G

        assert _17G == _17G + nG, f"17G + nG != 17G for {form}"
        assert (n + 1) * G == G, f"(n + 1)G != G for {form}"
        assert (n + 2) * G == 2 * G, f"(n + 2)G != 2G for {form}"

        logger.info(f"✅ Test elliptic curve cyclic group on real {form} curve done.")


def test_elliptic_curve_cyclic_group_on_test_curve():

    curves = ["weierstrass", "koblitz", "edwards"]

    for form in curves:
        logger.debug(f"ℹ️ Testing elliptic curve cyclic group on {form} test curve")
        cs = LightECC(form_name=form, curve_name="test-curve")

        for k in range(0, 2 * cs.curve.n + 1):
            P = cs.curve.double_and_add(cs.curve.G, k)
            logger.debug(f"{k} x G = {P}")

            if k in [0, cs.curve.n]:
                assert P == cs.curve.O

        logger.info(f"✅ Test elliptic curve cyclic group on test {form} curve done.")


def test_weierstrass_point_addition_returning_point_at_infinity():
    cs = LightECC(form_name="weierstrass", curve_name="test-curve")

    # we know that 20G + 8G = 28G = point at infinity
    P = cs.curve.add_points(
        cs.curve.double_and_add(cs.curve.G, 20),
        cs.curve.double_and_add(cs.curve.G, 8),
    )
    assert P == cs.curve.O

    _14G = cs.curve.double_and_add(cs.curve.G, 14)
    Q = cs.curve.add_points(_14G, _14G)
    assert Q == cs.curve.O

    logger.info("✅ Test weierstras point addition returning point at infinity done.")


def test_koblitz_point_addition_returning_point_at_infinity():
    cs = LightECC(form_name="koblitz", curve_name="test-curve")

    # we know that 12G + 4G = 16G = point at infinity
    P = cs.curve.add_points(
        cs.curve.double_and_add(cs.curve.G, 12),
        cs.curve.double_and_add(cs.curve.G, 4),
    )
    assert P == cs.curve.O

    _8G = cs.curve.double_and_add(cs.curve.G, 8)
    Q = cs.curve.add_points(_8G, _8G)
    assert Q == cs.curve.O

    logger.info("✅ Test koblitz point addition returning point at infinity done.")


def test_edwards_point_addition_returning_point_at_infinity():
    cs = LightECC(form_name="edwards", curve_name="test-curve")

    # we know that 6G + 2G = 8G = point at infinity
    P = cs.curve.add_points(
        cs.curve.double_and_add(cs.curve.G, 6),
        cs.curve.double_and_add(cs.curve.G, 2),
    )
    assert P == cs.curve.O

    _4G = cs.curve.double_and_add(cs.curve.G, 4)
    Q = cs.curve.add_points(_4G, _4G)
    assert Q == cs.curve.O

    logger.info("✅ Test edwards point addition returning point at infinity done.")


def test_double_and_add_for_k_close_to_n():
    for form in FORMS:
        cs = LightECC(form_name=form)

        _ = cs.curve.double_and_add(cs.curve.G, cs.curve.n - 1)
        assert cs.curve.double_and_add(cs.curve.G, cs.curve.n) == cs.curve.O
        assert cs.curve.double_and_add(cs.curve.G, cs.curve.n + 1) == cs.curve.G

        logger.info(
            f"✅ Double and add for k being close to order test done for {form}"
        )


def test_add_neutral_point():
    for form in FORMS:
        cs = LightECC(form_name=form)

        _7G = cs.curve.double_and_add(cs.curve.G, 7)

        assert cs.curve.add_points(_7G, cs.curve.O) == _7G
        assert cs.curve.add_points(cs.curve.O, _7G) == _7G
        assert cs.curve.add_points(cs.curve.O, cs.curve.O) == cs.curve.O

        logger.info(f"✅ Adding neutral point test done for {form}")
