import math
from typing import List, Optional, Tuple

from lightecc.interfaces.elliptic_curve import EllipticCurve


class Projective(EllipticCurve):
    def __init__(self, config: Optional[dict] = None):
        """
        Create Elliptic Curve satisfying y^2 = x^3 + b over Z_m in projective coordinates

        Points are (X : Y : Z) triples with x = X / Z and y = Y / Z, and the point at
        infinity is (0 : 1 : 0). Additions use complete formulas, which work for all
        points including doubling and the point at infinity, without any division.

        So, the modulus m does not have to be prime. This suits rings such as Z_{n^2} with
        n = pq, where division is not always possible and many points reduce to the point
        at infinity modulo n, as in Elliptic Curve Paillier and KMOV cryptosystems.
        Affine arithmetic of weierstrass form fails for those points.

        Note: complete formulas assume the subgroup has no point of order 2, e.g. odd order.
        Otherwise, P + Q returns the invalid triple (0 : 0 : 0) when P - Q has order 2, and
        is_on_curve rejects it.

        Ref: Renes, J., Costello, C., Batina, L. (2016). Complete Addition Formulas
            for Prime Order Elliptic Curves. EUROCRYPT 2016.
        Args:
            config (dict): curve parameters
                - modulus (int): modulus m, not necessarily prime
                - b (int): curve coefficient
                - G (tuple): base point as (x, y) or (X, Y, Z)
                - n (int): optional order of G. Leave it unset when it is secret,
                    e.g. it is the private key in Elliptic Curve Paillier.
        """
        if config is None:
            raise ValueError("Projective curve requires a configuration dictionary")

        if not all(k in config for k in ("modulus", "b", "G")):
            raise ValueError(
                "Projective curve configuration must include 'modulus', 'b', and 'G'"
            )

        # equation parameters
        self.a = 0
        self.b = config["b"] % config["modulus"]

        # modulos
        self.modulo = config["modulus"]

        # base point G
        G = tuple(config["G"])
        if len(G) == 2:
            G = (G[0], G[1], 1)
        if len(G) != 3:
            raise ValueError(f"Base point must be (x, y) or (X, Y, Z) but {G} given")
        self.G = G

        # order of G - optional
        self.n = config.get("n")

        # point at infinity
        self.O = (0, 1, 0)

    def add_points(
        self, P: Tuple[int, int, int], Q: Tuple[int, int, int]
    ) -> Tuple[int, int, int]:
        """
        Find P + Q for two given projective points with complete formulas
            (Algorithm 7 of Renes-Costello-Batina for a = 0)
        Args:
            P (Tuple[int, int, int]): (X1 : Y1 : Z1)
            Q (Tuple[int, int, int]): (X2 : Y2 : Z2)
        Returns:
            P+Q (Tuple[int, int, int]): (X3 : Y3 : Z3)
        """
        mod = self.modulo
        b3 = 3 * self.b
        x1, y1, z1 = P
        x2, y2, z2 = Q

        t0 = (x1 * x2) % mod
        t1 = (y1 * y2) % mod
        t2 = (z1 * z2) % mod
        t3 = ((x1 + y1) * (x2 + y2) - t0 - t1) % mod
        t4 = ((y1 + z1) * (y2 + z2) - t1 - t2) % mod
        y3 = ((x1 + z1) * (x2 + z2) - t0 - t2) % mod
        t0 = (3 * t0) % mod
        t2 = (b3 * t2) % mod
        z3 = (t1 + t2) % mod
        t1 = (t1 - t2) % mod
        y3 = (b3 * y3) % mod
        x3 = (t3 * t1 - t4 * y3) % mod
        y3 = (t1 * z3 + y3 * t0) % mod
        z3 = (z3 * t4 + t0 * t3) % mod
        return x3, y3, z3

    def double_point(self, P: Tuple[int, int, int]) -> Tuple[int, int, int]:
        """
        Find 2P for a given projective point with complete formulas
        """
        return self.add_points(P, P)

    def negative_point(self, P: Tuple[int, int, int]) -> Tuple[int, int, int]:
        x, y, z = P
        return x % self.modulo, (-y) % self.modulo, z % self.modulo

    def is_on_curve(self, P: Tuple[int, int, int]) -> bool:
        """
        Check Y^2 Z = X^3 + b Z^3 for a projective point. Coordinates must not share a
        common factor with the modulus, so (0 : 0 : 0) is not a point.
        Args:
            P (Tuple[int, int, int]): (X : Y : Z)
        Returns
            is_on_curve (boolean): returns True if point is on the curve
        """
        x, y, z = P
        if math.gcd(math.gcd(x, y), math.gcd(z, self.modulo)) != 1:
            return False
        return (y * y * z - x**3 - self.b * z**3) % self.modulo == 0

    def double_and_add(self, G: Tuple[int, int, int], k: int) -> Tuple[int, int, int]:
        """
        Perform scalar multiplication. Unlike other forms, k is not reduced modulo n,
        so it works when the order is unknown or G is not the generator of the group.
        Args:
            G (Tuple[int, int, int]): (X : Y : Z)
            k (int): scalar value, negative values give k x G = |k| x (-G)
        Returns
            kxG (Tuple[int, int, int]): a point on the elliptic curve
        """
        if k < 0:
            return self.double_and_add(self.negative_point(G), -k)

        result = self.O
        for bit in bin(k)[2:]:
            result = self.double_point(result)
            if bit == "1":
                result = self.add_points(result, G)
        return result

    def precompute_powers(
        self, P: Tuple[int, int, int], bits: int
    ) -> List[Tuple[int, int, int]]:
        """
        Precompute 2^i x P for i in [0, bits) to multiply a fixed point faster
        Args:
            P (Tuple[int, int, int]): (X : Y : Z)
            bits (int): number of powers, i.e. max bit length of scalars
        Returns:
            powers (list): [P, 2P, 4P, ..., 2^(bits-1) P]
        """
        powers = [P]
        for _ in range(bits - 1):
            powers.append(self.double_point(powers[-1]))
        return powers

    def multiply_with_powers(
        self, powers: List[Tuple[int, int, int]], k: int
    ) -> Tuple[int, int, int]:
        """
        Multiply a fixed point with a scalar using its precomputed powers, which avoids doublings
        Args:
            powers (list): output of precompute_powers for the point
            k (int): non-negative scalar with bit length at most len(powers)
        Returns:
            kxP (Tuple[int, int, int])
        """
        if k < 0:
            raise ValueError("Scalar must be non-negative for precomputed multiplication")
        if k.bit_length() > len(powers):
            raise ValueError(
                f"Scalar has {k.bit_length()} bits but {len(powers)} powers are given"
            )

        result = self.O
        for i in range(k.bit_length()):
            if (k >> i) & 1:
                result = self.add_points(result, powers[i])
        return result
