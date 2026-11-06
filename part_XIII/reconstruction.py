"""Bounded rational reconstruction by the extended Euclidean algorithm.

Returns a Fraction a/b with |a|, |b| <= floor(sqrt(m/2)) and a == x*b
modulo m, or None. Callers still verify recovered network weights exactly.
"""
from fractions import Fraction
from math import isqrt


def rational_reconstruction(x, modulus):
    if modulus <= 1:
        raise ValueError('modulus must exceed one')
    bound = isqrt(modulus // 2)
    previous, remainder = modulus, x % modulus
    previous_denominator, denominator = 0, 1
    while remainder > bound:
        quotient = previous // remainder
        previous, remainder = remainder, previous - quotient * remainder
        previous_denominator, denominator = denominator, previous_denominator - quotient * denominator
    if not denominator or abs(denominator) > bound:
        return None
    result = Fraction(remainder, denominator)
    if (result.numerator - x * result.denominator) % modulus:
        return None
    return result
