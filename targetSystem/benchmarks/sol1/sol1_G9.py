"""
https://projecteuler.net/problem=234

For an integer n ≥ 4, we define the lower prime square root of n, denoted by
lps(n), as the largest prime ≤ √n and the upper prime square root of n, ups(n),
as the smallest prime ≥ √n.

So, for example, lps(4) = 2 = ups(4), lps(1000) = 31, ups(1000) = 37. Let us
call an integer n ≥ 4 semidivisible, if one of lps(n) and ups(n) divides n,
but not both.

The sum of the semidivisible numbers not exceeding 15 is 30, the numbers are 8,
10 and 12. 15 is not semidivisible because it is a multiple of both lps(15) = 3
and ups(15) = 5. As a further example, the sum of the 92 semidivisible numbers
up to 1000 is 34825.

What is the sum of all semidivisible numbers not exceeding 999966663333 ?

G9 OPTIMIZATION: When possible, reduce data dependencies between iterations 
of the same code section to make the code parallelizable.
"""

import math


def prime_sieve(n: int) -> list:
    """
    Sieve of Erotosthenes
    Function to return all the prime numbers up to a certain number
    https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
    >>> prime_sieve(3)
    [2]
    >>> prime_sieve(50)
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    """
    is_prime = [True] * n
    is_prime[0] = False
    is_prime[1] = False
    is_prime[2] = True

    for i in range(3, int(n**0.5 + 1), 2):
        index = i * 2
        while index < n:
            is_prime[index] = False
            index = index + i

    primes = [2]

    for i in range(3, n, 2):
        if is_prime[i]:
            primes.append(i)

    return primes


def compute_interval_contribution(p: int, q: int, limit: int) -> int:
    """
    G9: Compute contribution for a single prime interval (p^2, q^2)
    This function is independent and can be parallelized
    """
    p_squared = p * p
    q_squared = q * q
    
    if p_squared > limit:
        return 0
    
    hi = min(limit, q_squared - 1)
    if hi <= p_squared:
        return 0
    
    contribution = 0
    
    # Get numbers divisible by lps(current)
    current = p_squared + p
    while current <= hi:
        contribution += current
        current += p

    # Add the numbers divisible by ups(current)
    current = hi
    while current > p_squared and current % q == 0:
        current -= q
    while current > p_squared:
        contribution += current
        current -= q

    # Remove the numbers divisible by both ups and lps
    current = p * q
    while current <= hi:
        if current > p_squared:
            contribution -= 2 * current
        current += p * q

    return contribution


def solution(limit: int = 999_966_663_333) -> int:
    """
    Computes the solution to the problem up to the specified limit
    >>> solution(1000)
    34825

    >>> solution(10_000)
    1134942

    >>> solution(100_000)
    36393008
    """
    primes_upper_bound = math.floor(math.sqrt(limit)) + 100
    primes = prime_sieve(primes_upper_bound)

    matches_sum = 0
    
    # G9: Each prime interval is independent - can be parallelized
    for i in range(len(primes) - 1):
        p = primes[i]
        q = primes[i + 1]
        
        # G9: Each interval computation is independent
        contribution = compute_interval_contribution(p, q, limit)
        matches_sum += contribution

    return matches_sum


if __name__ == "__main__":
    print(solution())
