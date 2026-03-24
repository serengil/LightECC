"""
Pairing implementations using Miller's algorithm for Weierstrass curves.

Supports:
    - Weil pairing over F_p (for curves where full r-torsion is in E(F_p))
    - Modified Tate pairing with distortion map over F_{p^2}
      (for supersingular curves y^2 = x^3 + ax with embedding degree 2)

The pairing e_r(P, Q) maps two r-torsion points P, Q on an elliptic curve
to an r-th root of unity.

Properties:
    - Bilinear: e(aP, bQ) = e(P, Q)^(ab)
    - Alternating: e(P, P) = 1
    - Non-degenerate: if e(P, Q) = 1 for all Q, then P = O

Reference: Washington, "Elliptic Curves: Number Theory and Cryptography", Chapter 11
"""

from typing import Tuple, Union

from lightecc.interfaces.elliptic_curve import EllipticCurvePoint
from lightecc.commons.logger import Logger

logger = Logger(module="lightecc/pairing.py")


# ===========================================================================
# F_{p^2} arithmetic: elements are (a, b) representing a + b*i where i^2 = -1
# Used for supersingular curves with embedding degree 2 and p ≡ 3 (mod 4)
# ===========================================================================


def _fp2_add(x, y, p):
    """Add two F_{p^2} elements."""
    return ((x[0] + y[0]) % p, (x[1] + y[1]) % p)


def _fp2_sub(x, y, p):
    """Subtract two F_{p^2} elements."""
    return ((x[0] - y[0]) % p, (x[1] - y[1]) % p)


def _fp2_mul(x, y, p):
    """Multiply two F_{p^2} elements. (a+bi)(c+di) = (ac-bd) + (ad+bc)i"""
    a, b = x
    c, d = y
    return ((a * c - b * d) % p, (a * d + b * c) % p)


def _fp2_inv(x, p):
    """Invert an F_{p^2} element. (a+bi)^{-1} = (a-bi)/(a^2+b^2)"""
    a, b = x
    norm = (a * a + b * b) % p
    if norm == 0:
        raise ZeroDivisionError("Cannot invert zero in F_{p^2}")
    norm_inv = pow(norm, -1, p)
    return ((a * norm_inv) % p, ((-b) * norm_inv) % p)


def _fp2_pow(x, n, p):
    """Exponentiate an F_{p^2} element by a non-negative integer n."""
    if n == 0:
        return (1, 0)
    if n < 0:
        x = _fp2_inv(x, p)
        n = -n
    result = (1, 0)
    base = x
    while n > 0:
        if n & 1:
            result = _fp2_mul(result, base, p)
        base = _fp2_mul(base, base, p)
        n >>= 1
    return result


# ===========================================================================
# Supersingular curve detection
# ===========================================================================


def _is_supersingular_embedding2(curve) -> bool:
    """
    Check if a Weierstrass curve y^2 = x^3 + ax + b is supersingular
    with embedding degree 2.

    For y^2 = x^3 + ax (b=0, a!=0) with p ≡ 3 (mod 4), the curve is
    supersingular with embedding degree 2, meaning p^2 ≡ 1 (mod r)
    but p ≢ 1 (mod r) for the relevant torsion orders.
    """
    return curve.b == 0 and curve.a != 0 and curve.modulo % 4 == 3


# ===========================================================================
# Standard Weil pairing (F_p) — for curves with full torsion in E(F_p)
# ===========================================================================


def _is_identity(P: EllipticCurvePoint) -> bool:
    """Check if a point is the identity element (point at infinity)."""
    return P.x == P.curve.O[0] and P.y == P.curve.O[1]


def _line_eval(
    T: EllipticCurvePoint, P: EllipticCurvePoint, R: EllipticCurvePoint
) -> int:
    """
    Evaluate the line function g_{T,P}(R) used in Miller's algorithm.

    g_{T,P} = l_{T,P} / v_{T+P} where:
        l_{T,P} is the line through T and P (tangent if T == P)
        v_{T+P} is the vertical line through the sum T + P

    Args:
        T: first point on the line
        P: second point on the line (same as T for tangent/doubling)
        R: evaluation point

    Returns:
        value of g_{T,P}(R) in F_p
    """
    if _is_identity(T) or _is_identity(P) or _is_identity(R):
        return 1

    a = T.curve.a
    p = T.curve.modulo

    xT, yT = T.x, T.y
    xP, yP = P.x, P.y
    xR, yR = R.x, R.y

    # compute T + P using the library's point arithmetic
    sum_point = T + P

    if _is_identity(sum_point):
        # T + P = O means vertical line (tangent at 2-torsion or T = -P)
        return (xR - xT) % p

    # slope of line through T and P (tangent if T == P, secant otherwise)
    if T == P:
        lam = (3 * xT * xT + a) * pow(2 * yT, -1, p) % p
    else:
        lam = (yP - yT) * pow(xP - xT, -1, p) % p

    # Line l evaluated at R: yR - yT - lam * (xR - xT)
    num = (yR - yT - lam * (xR - xT)) % p

    # Vertical line v through T+P evaluated at R: xR - x_{T+P}
    den = (xR - sum_point.x) % p

    if den == 0:
        return num if num != 0 else 1

    return num * pow(den, -1, p) % p


def _miller(P: EllipticCurvePoint, R: EllipticCurvePoint, r: int) -> int:
    """
    Miller's algorithm to compute f_{r,P}(R) over F_p.

    Args:
        P: base point (r-torsion)
        R: evaluation point
        r: torsion order

    Returns:
        f_{r,P}(R) in F_p
    """
    p = P.curve.modulo

    if _is_identity(P) or _is_identity(R):
        return 1

    f = 1
    T = P
    bits = bin(r)[2:]

    for i in range(1, len(bits)):
        # Doubling step
        if _is_identity(T):
            f = (f * f) % p
        else:
            g = _line_eval(T, T, R)
            f = (f * f * g) % p
            T = T + T

        if bits[i] == "1":
            # Addition step
            if _is_identity(T):
                T = P
            else:
                g = _line_eval(T, P, R)
                f = (f * g) % p
                T = T + P

    return f % p


def _find_auxiliary_point(
    P: EllipticCurvePoint, Q: EllipticCurvePoint
) -> EllipticCurvePoint:
    """
    Find an auxiliary point S on the curve for the Weil pairing computation.

    S must avoid degeneracies: S should not be in {O, P, -P, Q, -Q, P-Q, Q-P, P+Q}
    and the derived points Q+S, -S, P-S must also not be O.

    Args:
        P: first pairing input point
        Q: second pairing input point

    Returns:
        a safe auxiliary point S
    """
    curve = P.curve
    G = EllipticCurvePoint(curve.G[0], curve.G[1], curve)

    bad_points = [P, -P, Q, -Q, P - Q, Q - P, P + Q]

    for k in range(1, min(curve.n, 1000)):
        S = k * G

        if _is_identity(S):
            continue

        if any(S == bp for bp in bad_points):
            continue

        QpS = Q + S
        PmS = P - S

        if _is_identity(QpS) or _is_identity(PmS):
            continue

        return S

    raise ValueError(
        "Could not find a suitable auxiliary point for the Weil pairing. "
        "The curve may have too few points."
    )


# ===========================================================================
# Modified Tate pairing with distortion map (F_{p^2})
# For supersingular curves y^2 = x^3 + ax with embedding degree 2
# ===========================================================================


def _line_eval_fp2(T: EllipticCurvePoint, P: EllipticCurvePoint, xR, yR, p, a):
    """
    Evaluate the line function g_{T,P} at a point R = (xR, yR) in F_{p^2}.

    T and P are points on E(F_p). R is a point in E(F_{p^2}) given by
    its coordinates as F_{p^2} elements (tuples).

    Returns an F_{p^2} element (tuple).
    """
    if _is_identity(T) or _is_identity(P):
        return (1, 0)

    xT, yT = T.x, T.y
    xP, yP = P.x, P.y

    sum_point = T + P

    if _is_identity(sum_point):
        # Vertical line: xR - xT
        return _fp2_sub(xR, (xT, 0), p)

    # Slope lambda (in F_p, embedded as (lam, 0))
    if T == P:
        lam = (3 * xT * xT + a) * pow(2 * yT, -1, p) % p
    else:
        lam = (yP - yT) * pow(xP - xT, -1, p) % p

    # Numerator: yR - yT - lam * (xR - xT)
    num = _fp2_sub(
        _fp2_sub(yR, (yT, 0), p),
        _fp2_mul((lam, 0), _fp2_sub(xR, (xT, 0), p), p),
        p,
    )

    # Denominator: xR - x_{T+P}
    den = _fp2_sub(xR, (sum_point.x, 0), p)

    if den == (0, 0):
        return num if num != (0, 0) else (1, 0)

    return _fp2_mul(num, _fp2_inv(den, p), p)


def _miller_fp2(P: EllipticCurvePoint, xR, yR, r: int, p: int, a: int):
    """
    Miller's algorithm computing f_{r,P}(R) where R is in E(F_{p^2}).

    P is a point on E(F_p). R = (xR, yR) with coordinates in F_{p^2}.

    Returns an F_{p^2} element.
    """
    if _is_identity(P):
        return (1, 0)

    f = (1, 0)
    T = P
    bits = bin(r)[2:]

    for i in range(1, len(bits)):
        # Doubling step
        if _is_identity(T):
            f = _fp2_mul(f, f, p)
        else:
            g = _line_eval_fp2(T, T, xR, yR, p, a)
            f = _fp2_mul(_fp2_mul(f, f, p), g, p)
            T = T + T

        if bits[i] == "1":
            # Addition step
            if _is_identity(T):
                T = P
            else:
                g = _line_eval_fp2(T, P, xR, yR, p, a)
                f = _fp2_mul(f, g, p)
                T = T + P

    return f


def _tate_pairing_supersingular(
    P: EllipticCurvePoint, Q: EllipticCurvePoint, r: int
):
    """
    Modified Tate pairing for supersingular curves y^2 = x^3 + ax
    with embedding degree 2, using the distortion map phi(x, y) = (-x, i*y).

    Computes e_r(P, Q) = f_{r,P}(phi(Q))^{(p^2-1)/r} in F_{p^2}.

    Args:
        P: first r-torsion point on E(F_p)
        Q: second r-torsion point on E(F_p)
        r: torsion order

    Returns:
        An F_{p^2} element (a, b) representing a + b*i, which is an r-th root of unity.
    """
    p = P.curve.modulo
    a = P.curve.a

    # Apply distortion map: phi(Q) = (-x_Q, i * y_Q)
    # In F_{p^2}: x_phi = (-x_Q mod p, 0), y_phi = (0, y_Q)
    x_phi = ((-Q.x) % p, 0)
    y_phi = (0, Q.y % p)

    # Run Miller's algorithm over F_{p^2}
    f = _miller_fp2(P, x_phi, y_phi, r, p, a)

    # Final exponentiation: f^{(p^2 - 1) / r}
    exp = (p * p - 1) // r
    result = _fp2_pow(f, exp, p)

    return result


# ===========================================================================
# Public API
# ===========================================================================


def weil_pairing(
    P: EllipticCurvePoint, Q: EllipticCurvePoint, r: int
) -> Union[int, Tuple[int, int]]:
    """
    Compute the pairing e_r(P, Q).

    For supersingular curves with embedding degree 2 (y^2 = x^3 + ax, p ≡ 3 mod 4),
    automatically uses the modified Tate pairing with distortion map over F_{p^2}.

    For other curves, uses the standard Weil pairing over F_p.

    Args:
        P: first r-torsion point
        Q: second r-torsion point
        r: torsion order (r*P = O and r*Q = O must hold)

    Returns:
        int: an r-th root of unity in F_p (standard curves)
        tuple[int, int]: an element (a, b) representing a + b*i
            in F_{p^2} (supersingular curves with embedding degree 2)
    """
    p = P.curve.modulo

    if _is_identity(P) or _is_identity(Q):
        return 1

    # Use modified Tate pairing for supersingular curves with embedding degree 2,
    # but only when the full r-torsion is NOT in E(F_p).
    # If r | (p-1), the full r-torsion is in E(F_p) and the standard Weil pairing works.
    # Note: P == Q is allowed here because the distortion map phi(Q) != Q,
    # so e(P, phi(P)) is non-trivial (the alternating property does not apply).
    if _is_supersingular_embedding2(P.curve) and (p - 1) % r != 0:
        # Verify P and Q are r-torsion points
        rP = r * P
        if not _is_identity(rP):
            raise ValueError(
                f"P = ({P.x}, {P.y}) is not an {r}-torsion point"
                f" (r*P = ({rP.x}, {rP.y}))"
            )
        if P != Q:
            rQ = r * Q
            if not _is_identity(rQ):
                raise ValueError(
                    f"Q = ({Q.x}, {Q.y}) is not an {r}-torsion point"
                    f" (r*Q = ({rQ.x}, {rQ.y}))"
                )
        return _tate_pairing_supersingular(P, Q, r)

    if P == Q:
        return 1

    # Verify P and Q are r-torsion points
    rP = r * P
    rQ = r * Q
    if not _is_identity(rP):
        raise ValueError(
            f"P = ({P.x}, {P.y}) is not an {r}-torsion point"
            f" (r*P = ({rP.x}, {rP.y}))"
        )
    if not _is_identity(rQ):
        raise ValueError(
            f"Q = ({Q.x}, {Q.y}) is not an {r}-torsion point"
            f" (r*Q = ({rQ.x}, {rQ.y}))"
        )

    # Standard Weil pairing for curves with full torsion in E(F_p)
    S = _find_auxiliary_point(P, Q)
    nS = -S

    QpS = Q + S
    PmS = P - S

    fP_QpS = _miller(P, QpS, r)
    fP_S = _miller(P, S, r)
    fQ_PmS = _miller(Q, PmS, r)
    fQ_nS = _miller(Q, nS, r)

    if fP_S == 0 or fQ_PmS == 0:
        raise ValueError(
            "Degenerate case in Weil pairing computation. Try different input points."
        )

    num = (fP_QpS * fQ_nS) % p
    den = (fP_S * fQ_PmS) % p

    if den == 0:
        raise ValueError("Degenerate case: denominator is zero in Weil pairing.")

    return num * pow(den, -1, p) % p
