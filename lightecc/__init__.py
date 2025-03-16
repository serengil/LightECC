# built-in dependencies
from typing import Optional

# project dependencies
from lightecc.forms.weierstrass import Weierstrass
from lightecc.forms.edwards import TwistedEdwards
from lightecc.forms.koblitz import Koblitz
from lightecc.interfaces.elliptic_curve import EllipticCurvePoint
from lightecc.commons.logger import Logger

logger = Logger(module="lightecc/__init__.py")

VERSION = "0.0.2"


# pylint: disable=too-few-public-methods
class LightECC:
    __version__ = VERSION

    def __init__(
        self, form_name: Optional[str] = None, curve_name: Optional[str] = None
    ):
        """
        Construct an Elliptic Curve over a finite field (prime or binary)
        Args:
            form_name (str): specifies the form of the elliptic curve.
                Options: 'weierstrass' (default), 'edwards', 'koblitz'.
            curve_name (str): specifies the elliptic curve to use.
                Options:
                 - e.g. ed25519, ed448 for edwards form
                 - e.g. secp256k1 for weierstrass form
                 - e.g. k-409 for koblitz form
                List of all available curves:
                    github.com/serengil/LightECC
        """
        if form_name is None or form_name == "weierstrass":
            self.curve = Weierstrass(curve=curve_name)
        elif form_name in "edwards":
            self.curve = TwistedEdwards(curve=curve_name)
        elif form_name in "koblitz":
            self.curve = Koblitz(curve=curve_name)
        else:
            raise ValueError(f"unimplemented curve form - {form_name}")

        # base point
        self.G = EllipticCurvePoint(self.curve.G[0], self.curve.G[1], self.curve)

        # order of the curve
        self.n = self.curve.n

        # point at infinity or neutral / identity element
        self.O = EllipticCurvePoint(self.curve.O[0], self.curve.O[1], self.curve)

        # modulo
        self.modulo = self.curve.modulo
