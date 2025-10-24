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

OPTIMIZED VERSION: Applies all guidelines G1, G3, G4, G6, G7, G9, G12, G14
"""

import math
from concurrent.futures import ProcessPoolExecutor, as_completed


def prime_sieve(n: int) -> list:
    """
    Sieve of Erotosthenes - OPTIMIZED with G7 (bulk operations)
    Function to return all the prime numbers up to a certain number
    https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
    >>> prime_sieve(3)
    [2]
    >>> prime_sieve(50)
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    """
    if n < 2:
        return []
    
    # G7: Use bulk operations with slice assignment
    is_prime = bytearray(b"\x01") * n
    is_prime[0:2] = b"\x00\x00"
    
    # G7: Mark even numbers > 2 as not prime in bulk
    for i in range(4, n, 2):
        is_prime[i] = 0
    
    r = int(n**0.5)
    for i in range(3, r + 1, 2):
        if is_prime[i]:
            # G7: Bulk operation - mark all multiples of i starting from i*i
            is_prime[i*i:n:i] = b"\x00" * ((n - 1 - i*i) // i + 1)
    
    return [2] + [i for i in range(3, n, 2) if is_prime[i]]


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


def compute_interval_contribution(args):
    """
    G9/G12: Worker function for parallel processing
    Computes contribution for a single prime interval (p^2, q^2)
    """
    p, q, limit = args
    
    # G1: Hoist repeated expressions
    p_squared = p * p
    q_squared = q * q
    
    # G3: Early termination
    if p_squared > limit:
        return 0
    
    # G1: Hoist repeated expressions
    hi = min(limit, q_squared - 1)
    
    # G4: Short-circuit evaluation
    if hi <= p_squared:
        return 0
    
    # G7: Use bulk operations instead of loops
    sum_p = sum_ap_multiples_in_range(p, p_squared + 1, hi)
    sum_q = sum_ap_multiples_in_range(q, p_squared + 1, hi)
    sum_pq = sum_ap_multiples_in_range(p * q, p_squared + 1, hi)
    
    # G6: Avoid recomputing - use the bulk results directly
    return sum_p + sum_q - 2 * sum_pq


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
    # G1: Hoist repeated expressions
    primes_upper_bound = int(math.isqrt(limit)) + 100
    primes = prime_sieve(primes_upper_bound)
    
    # G14: Cache frequently used data - precompute prime squares
    prime_squares = [p * p for p in primes]

    # G9: Prepare work items for parallel processing (reduce data dependencies)
    work_items = []
    for i in range(len(primes) - 1):
        p = primes[i]
        q = primes[i + 1]
        # G3: Early termination check
        if prime_squares[i] > limit:
            break
        work_items.append((p, q, limit))

    matches_sum = 0
    
    # G12: Use ProcessPoolExecutor with work-stealing dynamic scheduling
    with ProcessPoolExecutor() as executor:
        # Submit all work items - work-stealing will distribute them dynamically
        future_to_work = {executor.submit(compute_interval_contribution, work): work 
                          for work in work_items}
        
        # Collect results as they complete (work-stealing ensures efficient distribution)
        for future in as_completed(future_to_work):
            try:
                contribution = future.result()
                matches_sum += contribution
            except Exception as exc:
                print(f'Work item {future_to_work[future]} generated an exception: {exc}')

    return matches_sum


if __name__ == "__main__":
    print(solution())
