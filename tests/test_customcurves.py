# 3rd party dependencies
import pytest

# project dependencies
from lightecc import LightECC
from lightecc.curves import inventory
from lightecc.curves.weierstrass import Secp256k1, SS_Weierstrass_307
from lightecc.curves.edwards import Ed25519
from lightecc.curves.koblitz import K163
from lightecc.commons.errors import InvalidCurveOrder, PointNotOnCurve
from lightecc.commons.logger import Logger

logger = Logger(module="tests/test_customcurves.py")

FORMS = ["weierstrass", "edwards", "koblitz"]


def test_build_curves():
    for form in FORMS:
        curves = inventory.list_curves(form)
        for curve in curves:
            if curve in [
                "customweierstrasscurve",
                "customedwardscurve",
                "customkoblitzcurve",
            ]:
                continue

            ec = LightECC(form_name=form, curve_name=curve)

            # base point
            G = ec.G

            # additions
            _2G = G + G
            _3G = G + _2G
            _5G = _2G + _3G
            _10G = _5G + _5G
            _20G = _10G + _10G

            # subtractions
            _18G = _20G - _2G
            _17G = _18G - G

            assert _17G == 17 * G
            assert _17G == G * 17

            # multiplications
            _34G = 34 * G

            # division
            assert _34G / _17G == 2

            logger.info(f"✅ {form}-{curve} tests done")


def test_custom_weierstrass_curve_with_wrong_generator():
    curve = Secp256k1()
    config = {
        "p": curve.p,
        "a": curve.a,
        "b": curve.b,
        "G": (curve.G[0], curve.G[1] + 1),  # wrong generator
        "n": curve.n,  # wrong order
    }
    with pytest.raises(PointNotOnCurve):
        LightECC(form_name="weierstrass", curve_name="custom", config=config)

    logger.info("✅ custom weierstrass curve with wrong generator test done")


def test_custom_weierstrass_curve_with_wrong_order():
    curve = Secp256k1()
    config = {
        "p": curve.p,
        "a": curve.a,
        "b": curve.b,
        "G": curve.G,
        "n": curve.n - 3,  # wrong order
    }
    with pytest.raises(
        InvalidCurveOrder,
        match="Invalid curve order: n x G does not equal to the point at infinity / identity element.",
    ):
        LightECC(form_name="weierstrass", curve_name="custom", config=config)

    logger.info("✅ custom weierstrass curve with wrong order test done")


def test_custom_edwards_curve_with_wrong_order():
    curve = Ed25519()
    config = {
        "p": curve.p,
        "a": curve.a,
        "d": curve.d,
        "G": curve.G,
        "n": curve.n - 3,  # wrong order
    }
    with pytest.raises(
        InvalidCurveOrder,
        match="Invalid curve order: n x G does not equal to the neutral element / identity element.",
    ):
        LightECC(form_name="edwards", curve_name="custom", config=config)

    logger.info("✅ custom edwards curve with wrong order test done")


def test_custom_koblitz_curve_with_wrong_order():
    curve = K163()
    config = {
        "m": curve.m,
        "coefficients": curve.coefficients,
        "a": curve.a,
        "b": curve.b,
        "G": curve.G,
        "n": curve.n - 3,  # wrong order
    }
    with pytest.raises(
        InvalidCurveOrder,
        match="Invalid curve order: n x G does not equal to the point at infinity / identity element.",
    ):
        LightECC(form_name="koblitz", curve_name="custom", config=config)

    logger.info("✅ custom koblitz curve with wrong order test done")


def test_custom_curve_with_none_order():
    # when n is None, order validation should be skipped
    config = {
        "p": 97,
        "a": 2,
        "b": 3,
        "G": (3, 6),
        "n": None,
    }
    ec = LightECC(form_name="weierstrass", curve_name="custom", config=config)
    assert ec.n is None

    logger.info("✅ custom curve with None order test done")


def test_find_order():
    curve = SS_Weierstrass_307()
    config = {
        "p": curve.p,
        "a": curve.a,
        "b": curve.b,
        "G": curve.G,
        "n": None,  # order will be found automatically
    }
    ec = LightECC(form_name="weierstrass", curve_name="custom", config=config)
    n = ec.find_order()
    assert n is not None
    assert n * ec.G == ec.O
    assert n == curve.n

    logger.info("✅ find order test done")
