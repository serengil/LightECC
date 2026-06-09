# built-in dependencies
import pickle

# project dependencies
from lightecc import LightECC
from lightecc.interfaces.elliptic_curve import EllipticCurvePoint
from lightecc.commons.logger import Logger

logger = Logger(module="tests/test_curve_equality.py")

FORMS = ["weierstrass", "edwards", "koblitz"]


def test_curve_equality_is_value_based():
    # two freshly constructed instances of the same curve are distinct
    # objects but must compare equal (value based, not object identity)
    for form in FORMS:
        ec1 = LightECC(form_name=form, curve_name="test-curve")
        ec2 = LightECC(form_name=form, curve_name="test-curve")

        assert ec1.curve is not ec2.curve
        assert ec1.curve == ec2.curve
        assert not (ec1.curve != ec2.curve)

        logger.info(f"✅ Value based curve equality tested for {form} form.")


def test_curves_with_different_parameters_are_not_equal():
    # different forms describe different curves
    ws = LightECC(form_name="weierstrass", curve_name="test-curve")
    ed = LightECC(form_name="edwards", curve_name="test-curve")
    assert ws.curve != ed.curve

    # same form but different named curve must differ as well
    ws_test = LightECC(form_name="weierstrass", curve_name="test-curve")
    ws_real = LightECC(form_name="weierstrass", curve_name="secp256k1")
    assert ws_test.curve != ws_real.curve

    logger.info("✅ Curves with different parameters are not equal.")


def test_equal_curves_share_hash_and_are_usable_in_sets():
    for form in FORMS:
        ec1 = LightECC(form_name=form, curve_name="test-curve")
        ec2 = LightECC(form_name=form, curve_name="test-curve")

        # equal objects must have equal hashes
        assert hash(ec1.curve) == hash(ec2.curve)

        # and must collapse to a single entry in a set
        assert len({ec1.curve, ec2.curve}) == 1

        logger.info(f"✅ Curve hashing tested for {form} form.")


def test_curve_survives_pickle_round_trip():
    for form in FORMS:
        ec = LightECC(form_name=form, curve_name="test-curve")
        restored = pickle.loads(pickle.dumps(ec.curve))

        assert restored is not ec.curve
        assert restored == ec.curve
        assert hash(restored) == hash(ec.curve)

        logger.info(f"✅ Curve pickle round-trip tested for {form} form.")


def test_point_equality_survives_pickle_round_trip():
    # regression test: an EllipticCurvePoint carries its curve, and pickling
    # rebuilds that curve as a fresh instance. Without value based curve
    # equality the restored point would never compare equal to a point
    # recomputed on the original curve, breaking any consumer that relies on
    # `==` after a point crosses a serialization boundary (e.g. multiprocessing).
    for form in FORMS:
        ec = LightECC(form_name=form, curve_name="test-curve")

        P = ec.G * 5
        restored = pickle.loads(pickle.dumps(P))

        assert restored is not P
        assert (restored.x, restored.y) == (P.x, P.y)
        assert restored == P

        # a point recomputed on the original curve must also match the
        # restored point - this mirrors the decryption comparison loop
        recomputed = ec.G * 5
        assert recomputed == restored

        logger.info(f"✅ Point equality across pickle tested for {form} form.")


def test_point_equality_across_independent_curve_instances():
    # same coordinates declared against two independent (but equal) curve
    # instances must compare equal
    for form in FORMS:
        ec1 = LightECC(form_name=form, curve_name="test-curve")
        ec2 = LightECC(form_name=form, curve_name="test-curve")

        P = ec1.G * 3
        Q = EllipticCurvePoint(x=P.x, y=P.y, curve=ec2.curve)

        assert P == Q

        logger.info(
            f"✅ Point equality across independent curve instances tested for {form} form."
        )
