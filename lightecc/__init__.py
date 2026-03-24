# built-in dependencies
from typing import Optional, Tuple, Union

# project dependencies
from lightecc.forms.weierstrass import Weierstrass
from lightecc.forms.edwards import TwistedEdwards
from lightecc.forms.koblitz import Koblitz
from lightecc.interfaces.elliptic_curve import EllipticCurvePoint
from lightecc.commons.pairing import weil_pairing
from lightecc.commons.logger import Logger

logger = Logger(module="lightecc/__init__.py")

VERSION = "0.0.5"


# pylint: disable=too-few-public-methods
class LightECC:
    __version__ = VERSION

    def __init__(
        self,
        form_name: Optional[str] = None,
        curve_name: Optional[str] = None,
        config: Optional[dict] = None,
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
            config (dict): optional configuration parameters for the curve.
                 This can include custom parameters for the curve, such as:
                 - 'a', 'b' for weierstrass curves
                 - 'd' for edwards curves
                 - 'k' for koblitz curves
                 If curve name is custom, then these parameters are required to define the curve.
        """
        if form_name is None or form_name == "weierstrass":
            self.curve = Weierstrass(curve=curve_name, config=config)
        elif form_name == "edwards":
            self.curve = TwistedEdwards(curve=curve_name, config=config)
        elif form_name == "koblitz":
            self.curve = Koblitz(curve=curve_name, config=config)
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

    def pairing(
        self, P: EllipticCurvePoint, Q: EllipticCurvePoint, n: Optional[int] = None
    ) -> Union[Tuple[int, int], int]:
        """
        Compute the pairing e_n(P, Q).

        For standard curves, computes the Weil pairing over F_p.
        For supersingular curves with embedding degree 2, computes the
        modified Tate pairing over F_{p^2} using a distortion map.

        Args:
            P: first n-torsion point on the curve
            Q: second n-torsion point on the curve
            n: torsion order (default: curve order).
                P and Q must satisfy n*P = O and n*Q = O.

        Returns:
            int: an n-th root of unity in F_p (standard curves)
            tuple[int, int]: an element (a, b) representing a + b*i
                in F_{p^2} (supersingular curves with embedding degree 2)
        """
        if not isinstance(self.curve, Weierstrass):
            raise ValueError("Pairing is only supported for Weierstrass form curves.")

        if n is None:
            n = self.n

        result = weil_pairing(P, Q, n)

        return result
