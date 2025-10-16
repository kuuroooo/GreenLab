"""
Project Euler Problem 234: https://projecteuler.net/problem=234

For any integer n, consider the three functions

f1,n(x,y,z) = x^(n+1) + y^(n+1) - z^(n+1)
f2,n(x,y,z) = (xy + yz + zx)*(x^(n-1) + y^(n-1) - z^(n-1))
f3,n(x,y,z) = xyz*(xn-2 + yn-2 - zn-2)

and their combination

fn(x,y,z) = f1,n(x,y,z) + f2,n(x,y,z) - f3,n(x,y,z)

We call (x,y,z) a golden triple of order k if x, y, and z are all rational numbers
of the form a / b with 0 < a < b ≤ k and there is (at least) one integer n,
so that fn(x,y,z) = 0.

Let s(x,y,z) = x + y + z.
Let t = u / v be the sum of all distinct s(x,y,z) for all golden triples
(x,y,z) of order 35.
All the s(x,y,z) and t must be in reduced form.

Find u + v.


Solution:

By expanding the brackets it is easy to show that
fn(x, y, z) = (x + y + z) * (x^n + y^n - z^n).

Since x,y,z are positive, the requirement fn(x, y, z) = 0 is fulfilled if and
only if x^n + y^n = z^n.

By Fermat's Last Theorem, this means that the absolute value of n can not
exceed 2, i.e. n is in {-2, -1, 0, 1, 2}. We can eliminate n = 0 since then the
equation would reduce to 1 + 1 = 1, for which there are no solutions.

So all we have to do is iterate through the possible numerators and denominators
of x and y, calculate the corresponding z, and check if the corresponding numerator and
denominator are integer and satisfy 0 < z_num < z_den <= 0. We use a set "uniquq_s"
to make sure there are no duplicates, and the fractions.Fraction class to make sure
we get the right numerator and denominator.

Reference:
https://en.wikipedia.org/wiki/Fermat%27s_Last_Theorem
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt
from functools import lru_cache


def add_three(
    x_num: int, x_den: int, y_num: int, y_den: int, z_num: int, z_den: int
) -> tuple[int, int]:
    """
    Given the numerators and denominators of three fractions, return the
    numerator and denominator of their sum in lowest form.
    >>> add_three(1, 3, 1, 3, 1, 3)
    (1, 1)
    >>> add_three(2, 5, 4, 11, 12, 3)
    (262, 55)
    """
    top: int = x_num * y_den * z_den + y_num * x_den * z_den + z_num * x_den * y_den
    bottom: int = x_den * y_den * z_den
    hcf: int = gcd(top, bottom)
    top //= hcf
    bottom //= hcf
    return top, bottom


@lru_cache(maxsize=1000)
def cached_gcd(a: int, b: int) -> int:
    """
    G14: Cache frequently used GCD calculations.
    """
    return gcd(a, b)


@lru_cache(maxsize=1000)
def cached_isqrt(n: int) -> int:
    """
    G14: Cache frequently used square root calculations.
    """
    return isqrt(n)


def solution(order: int = 35) -> int:
    """
    Find the sum of the numerator and denominator of the sum of all s(x,y,z) for
    golden triples (x,y,z) of the given order.

    >>> solution(5)
    296
    >>> solution(10)
    12519
    >>> solution(20)
    19408891927
    """
    # G14: Cache all valid fractions and their squares
    fractions = [(a, b) for b in range(2, order + 1) for a in range(1, b)]
    # G14: Cache squares for all fractions to avoid repeated computation
    squares = [(a * a, b * b) for a, b in fractions]
    
    # G14: Cache frequently used values
    order_cache = order
    
    unique_s: set = set()
    total: Fraction = Fraction(0)

    for i, (x_num, x_den) in enumerate(fractions):
        xn2, xd2 = squares[i]  # G14: Use cached squares
        for j, (y_num, y_den) in enumerate(fractions):
            yn2, yd2 = squares[j]  # G14: Use cached squares
            
            # n=1
            z_num = x_num * y_den + x_den * y_num
            z_den = x_den * y_den
            # G14: Use cached GCD calculation
            g = cached_gcd(z_num, z_den)
            z_num //= g
            z_den //= g
            if 0 < z_num < z_den <= order_cache:
                unique_s.add(add_three(
                    x_num, x_den, y_num, y_den, z_num, z_den
                ))

            # n=2
            z_num = xn2 * yd2 + xd2 * yn2  # G14: Use cached squares
            z_den = xd2 * yd2  # G14: Use cached squares
            if z_num > 0 and z_den > 0:
                # G14: Use cached square root calculation
                z_num_sqrt = cached_isqrt(z_num)
                z_den_sqrt = cached_isqrt(z_den)
                if z_num_sqrt * z_num_sqrt == z_num and z_den_sqrt * z_den_sqrt == z_den:
                    z_num = z_num_sqrt
                    z_den = z_den_sqrt
                    g = cached_gcd(z_num, z_den)
                    z_num //= g
                    z_den //= g
                    if 0 < z_num < z_den <= order_cache:
                        unique_s.add(add_three(
                            x_num, x_den, y_num, y_den, z_num, z_den
                        ))

            # n=-1
            z_num = x_num * y_num
            z_den = x_den * y_num + x_num * y_den
            g = cached_gcd(z_num, z_den)
            z_num //= g
            z_den //= g
            if 0 < z_num < z_den <= order_cache:
                unique_s.add(add_three(
                    x_num, x_den, y_num, y_den, z_num, z_den
                ))

            # n=-2
            z_num = xn2 * yn2  # G14: Use cached squares
            z_den = xd2 * yn2 + xn2 * yd2  # G14: Use cached squares
            if z_num > 0 and z_den > 0:
                # G14: Use cached square root calculation
                z_num_sqrt = cached_isqrt(z_num)
                z_den_sqrt = cached_isqrt(z_den)
                if z_num_sqrt * z_num_sqrt == z_num and z_den_sqrt * z_den_sqrt == z_den:
                    z_num = z_num_sqrt
                    z_den = z_den_sqrt
                    g = cached_gcd(z_num, z_den)
                    z_num //= g
                    z_den //= g
                    if 0 < z_num < z_den <= order_cache:
                        unique_s.add(add_three(
                            x_num, x_den, y_num, y_den, z_num, z_den
                        ))

    for num, den in unique_s:
        total += Fraction(num, den)

    return total.denominator + total.numerator


if __name__ == "__main__":
    print(f"{solution() = }")





