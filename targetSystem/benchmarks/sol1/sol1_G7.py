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

G7 OPTIMIZATION: Use bulk operations to reduce overhead from handling individual 
tasks, e.g., generating a batch of integers at once rather than one by one.
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


def sum_ap_multiples_in_range(k: int, L: int, R: int) -> int:
    """
    G7: Bulk operation - sum of multiples of k in [L, R], inclusive
    Uses arithmetic progression formula instead of individual iterations
    """
    if L > R:
        return 0
    first = ((L + k - 1) // k) * k
    if first > R:
        return 0
    last = (R // k) * k
    n = (last - first) // k + 1
    return n * (first + last) // 2


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
    prime_index = 0
    last_prime = primes[prime_index]

    while (last_prime**2) <= limit:
        next_prime = primes[prime_index + 1]

        lower_bound = last_prime**2
        upper_bound = next_prime**2

        # G7: Use bulk operations instead of individual loops
        # Get numbers divisible by lps(current) - use arithmetic progression
        hi = min(limit, upper_bound - 1)
        if hi >= lower_bound + last_prime:
            sum_p = sum_ap_multiples_in_range(last_prime, lower_bound + 1, hi)
            matches_sum += sum_p

        # Add the numbers divisible by ups(current) - use arithmetic progression
        if hi >= lower_bound + next_prime:
            sum_q = sum_ap_multiples_in_range(next_prime, lower_bound + 1, hi)
            matches_sum += sum_q

        # Remove the numbers divisible by both ups and lps - use arithmetic progression
        if hi >= lower_bound + last_prime * next_prime:
            sum_pq = sum_ap_multiples_in_range(last_prime * next_prime, lower_bound + 1, hi)
            matches_sum -= 2 * sum_pq

        # Setup for next pair
        last_prime = next_prime
        prime_index += 1

    return matches_sum


if __name__ == "__main__":
    print(solution())
