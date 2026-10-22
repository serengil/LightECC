# built-in dependencies
from typing import Tuple
from abc import ABC, abstractmethod

# project dependencies
from lightecc.commons.errors import PointNotOnCurve
from lightecc.commons.logger import Logger

logger = Logger(module="lightecc/interfaces/elliptic_curve.py")

# Signature for elliptic curve


class EllipticCurve(ABC):
    # point at infinity or neutral / identity element
    O: Tuple[int, int]

    # base point G
    G: Tuple[int, int]

    # modulo (prime p or polynomial fx)
    modulo: int

    # order of the curve
    n: int

    # coefficients
    a: int
    b: int = 0  # for weierstrass & koblitz form
    d: int = 0  # for edwards form

    @abstractmethod
    def add_points(
        self,
        P: Tuple[int, int],
        Q: Tuple[int, int],
    ) -> Tuple[int, int]:
        pass

    @abstractmethod
    def double_point(self, P: Tuple[int, int]) -> Tuple[int, int]:
        pass

    @abstractmethod
    def is_on_curve(self, P: Tuple[int, int]) -> bool:
        pass

    @abstractmethod
    def negative_point(self, P: Tuple[int, int]) -> Tuple[int, int]:
        pass

    def double_and_add(self, G: Tuple[int, int], k: int) -> Tuple[int, int]:
        """
        Perform scalar multiplication over elliptic curve
        Args:
            G (Tuple[int, int]): a point on an elliptic curve
            k (int): scalar value
        Returns
            kxG (Tuple[int, int]): a point on an elliptic curve
        """
        target_point = G

        if self.n is not None and k >= self.n:
            k = k % self.n
        if k == 0:
            return self.O
        if k < 0:
            return self.negative_point(self.double_and_add(G, abs(k)))

        k_binary = bin(k)[2:]

        for i in range(1, len(k_binary)):
            current_bit = k_binary[i : i + 1]

            # doubling - always
            target_point = self.double_point(target_point)

            if current_bit == "1":
                target_point = self.add_points(target_point, G)

        assert (
            self.is_on_curve(target_point) is True
        ), f"{target_point} is not on the curve!"

        return target_point

    def __eq__(self, other: object) -> bool:
        """
        Two curves are equal when they share the same form and defining
        parameters (equation coefficients, modulo and order). This is a
        value-based comparison rather than object identity so that a curve
        keeps comparing equal to itself after being pickled - e.g. when an
        EllipticCurvePoint is sent across a multiprocessing boundary.
        """
        if not isinstance(other, EllipticCurve):
            return NotImplemented
        return (
            type(self) is type(other)
            and self.modulo == other.modulo
            and self.n == other.n
            and self.a == other.a
            and self.b == other.b
            and self.d == other.d
        )

    def __hash__(self) -> int:
        return hash((type(self).__name__, self.modulo, self.n, self.a, self.b, self.d))


class EllipticCurvePoint:
    """
    Define a point on an elliptic curve
    """

    def __init__(self, x: int, y: int, curve: EllipticCurve):
        self.x = x
        self.y = y
        self.curve = curve

        if not self.curve.is_on_curve((x, y)):
            raise PointNotOnCurve(f"({x}, {y}) is not on the curve!")

    def get_point(self) -> Tuple[int, int]:
        return (self.x, self.y)

    def is_identity(self) -> bool:
        """
        Check if the point is the neutral / identity element O of its curve
        """
        return self.x == self.curve.O[0] and self.y == self.curve.O[1]

    def __repr__(self):

        if self.is_identity() and self.x == float("inf") and self.y == float("inf"):
            # because Edwards has neutral / identity element instead of point at infinity
            return "\U0001d4aa"  # unicode for "𝒪" (circle O)
        return f"({self.x}, {self.y})"

    def __str__(self):
        return self.__repr__()

    def __add__(self, other: "EllipticCurvePoint") -> "EllipticCurvePoint":
        """
        Calculate P + Q for two given points P and Q
        """
        if not isinstance(other, EllipticCurvePoint):
            raise ValueError("Addition is only defined for 2 points")

        if self.curve != other.curve:
            raise ValueError("Points are not on the same curve")

        x, y = self.curve.add_points((self.x, self.y), (other.x, other.y))
        return EllipticCurvePoint(x, y, self.curve)

    def __sub__(self, other: "EllipticCurvePoint") -> "EllipticCurvePoint":
        """
        Calculate P - Q for two given points P and Q
        """
        return self.__add__(other=other.__neg__())

    def __mul__(self, k: int) -> "EllipticCurvePoint":
        """
        Calculate k*P for a given k and P
        """
        if not isinstance(k, int):
            raise ValueError("Multiplication is only defined for an integer")

        # k x O = O for any k
        if self.is_identity():
            return self

        x, y = self.curve.double_and_add((self.x, self.y), k)
        return EllipticCurvePoint(x, y, self.curve)

    def __rmul__(self, k: int) -> "EllipticCurvePoint":
        """
        Calculate k*P for a given k and P
        """
        if not isinstance(k, int):
            raise ValueError("Multiplication is only defined for an integer")
        return self.__mul__(k)

    def __neg__(self) -> "EllipticCurvePoint":
        """
        Calculate -P for a given P
        """
        # -O = O
        if self.is_identity():
            return self

        x, y = self.curve.negative_point((self.x, self.y))
        return EllipticCurvePoint(x, y, self.curve)

    def __eq__(self, other: "EllipticCurvePoint") -> bool:
        """
        Check if two points are equal
        """
        return self.x == other.x and self.y == other.y and self.curve == other.curve

    def __truediv__(self, other: "EllipticCurvePoint") -> int:
        """
        Resolve ECDLP - this is a hard problem!
        """
        # TODO: you may consider to use baby-step-giant-step instead of brute force

        logger.debug(f"Find k from ({self.x}, {self.y}) = k x ({other.x}, {other.y})")

        # O = 0 x G
        if self.is_identity():
            return 0
        if self.x == other.x and self.y == other.y:
            return 1

        # base point
        gx, gy = other.get_point()

        k = 2
        while True:
            kG = self.curve.double_and_add((gx, gy), k)

            if kG[0] == self.x and kG[1] == self.y:
                return k

            k = k + 1

            if k > self.curve.n:
                raise ValueError(
                    f"Cannot restore scalar from ({self.x}, {self.y}) = k x ({other.x}, {other.y})"
                )


class ProjectivePoint:
    """
    Define a point (X : Y : Z) on an elliptic curve in projective coordinates
    """

    def __init__(self, x: int, y: int, z: int, curve: EllipticCurve):
        self.x = x % curve.modulo
        self.y = y % curve.modulo
        self.z = z % curve.modulo
        self.curve = curve

        if not self.curve.is_on_curve((self.x, self.y, self.z)):
            raise PointNotOnCurve(f"({x} : {y} : {z}) is not on the curve!")

    def get_point(self) -> Tuple[int, int, int]:
        return (self.x, self.y, self.z)

    def is_identity(self) -> bool:
        """
        Check if the point is the point at infinity (0 : 1 : 0)
        """
        return self.x == 0 and self.z == 0

    def to_affine(self) -> Tuple[int, int]:
        """
        Convert the point to affine (x, y). Raises ValueError if Z is not invertible.
        """
        x, y, z = self.get_point()
        try:
            z_inv = pow(z, -1, self.curve.modulo)
        except ValueError as err:
            raise ValueError(f"{self} cannot be written as an affine point") from err
        return (x * z_inv) % self.curve.modulo, (y * z_inv) % self.curve.modulo

    def __repr__(self):
        if self.is_identity():
            return "\U0001d4aa"  # unicode for "𝒪" (circle O)
        return f"({self.x} : {self.y} : {self.z})"

    def __str__(self):
        return self.__repr__()

    def __add__(self, other: "ProjectivePoint") -> "ProjectivePoint":
        """
        Calculate P + Q for two given points P and Q
        """
        if not isinstance(other, ProjectivePoint):
            raise ValueError("Addition is only defined for 2 points")

        if self.curve != other.curve:
            raise ValueError("Points are not on the same curve")

        x, y, z = self.curve.add_points(self.get_point(), other.get_point())
        return ProjectivePoint(x, y, z, self.curve)

    def __sub__(self, other: "ProjectivePoint") -> "ProjectivePoint":
        """
        Calculate P - Q for two given points P and Q
        """
        return self.__add__(other=other.__neg__())

    def __mul__(self, k: int) -> "ProjectivePoint":
        """
        Calculate k*P for a given k and P
        """
        if not isinstance(k, int):
            raise ValueError("Multiplication is only defined for an integer")

        x, y, z = self.curve.double_and_add(self.get_point(), k)
        return ProjectivePoint(x, y, z, self.curve)

    def __rmul__(self, k: int) -> "ProjectivePoint":
        """
        Calculate k*P for a given k and P
        """
        if not isinstance(k, int):
            raise ValueError("Multiplication is only defined for an integer")
        return self.__mul__(k)

    def __neg__(self) -> "ProjectivePoint":
        """
        Calculate -P for a given P
        """
        x, y, z = self.curve.negative_point(self.get_point())
        return ProjectivePoint(x, y, z, self.curve)

    def __eq__(self, other: object) -> bool:
        """
        Check if two points are equal. (X1 : Y1 : Z1) and (X2 : Y2 : Z2) are the same
        point when one is a scalar multiple of the other, i.e. all 2x2 minors vanish.
        Checking only X1 Z2 = X2 Z1 is not enough over Z_{n^2}, e.g. for (0 : 1 : 0)
        and (n : 1 : 0).
        """
        if not isinstance(other, ProjectivePoint):
            return NotImplemented
        if self.curve != other.curve:
            return False

        m = self.curve.modulo
        x1, y1, z1 = self.get_point()
        x2, y2, z2 = other.get_point()
        return (
            (x1 * y2 - x2 * y1) % m == 0
            and (x1 * z2 - x2 * z1) % m == 0
            and (y1 * z2 - y2 * z1) % m == 0
        )
