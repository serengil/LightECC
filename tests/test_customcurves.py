# project dependencies
from lightecc import LightECC
from lightecc.curves import inventory
from lightecc.commons.logger import Logger

logger = Logger(module="tests/test_customcurves.py")

FORMS = ["weierstrass", "edwards", "koblitz"]


def test_build_curves():
    for form in FORMS:
        curves = inventory.list_curves(form)
        for curve in curves:
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
