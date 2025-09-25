"""
Simple multithreaded algorithm to show how the 4 phases of a genetic algorithm works
(Evaluation, Selection, Crossover and Mutation)
https://en.wikipedia.org/wiki/Genetic_algorithm
Author: D4rkia
G7 Optimization: Use bulk operations to reduce overhead from handling individual tasks
"""

from __future__ import annotations

import random
import concurrent.futures
from functools import lru_cache
from typing import Dict, List, Tuple

# Maximum size of the population.  Bigger could be faster but is more memory expensive.
N_POPULATION = 200
# Number of elements selected in every generation of evolution. The selection takes
# place from best to worst of that generation and must be smaller than N_POPULATION.
N_SELECTED = 50
# Probability that an element of a generation can mutate, changing one of its genes.
# This will guarantee that all genes will be used during evolution.
MUTATION_PROBABILITY = 0.4
# Just a seed to improve randomness required by the algorithm.
random.seed(random.randint(0, 1000))

# Cache for frequently used computations
EVALUATION_CACHE: Dict[str, float] = {}
POPULATION_THIRD = N_POPULATION // 3


@lru_cache(maxsize=10000)
def evaluate(item: str, main_target: str) -> tuple[str, float]:
    """
    Evaluate how similar the item is with the target by just
    counting each char in the right position
    >>> evaluate("Helxo Worlx", "Hello World")
    ('Helxo Worlx', 9.0)
    """
    # G1: Store repeated expressions in variables
    target_length = len(main_target)
    
    # G4: Use short-circuit evaluation for early termination
    if len(item) != target_length:
        return (item, 0.0)
    
    # G6: Avoid redundant computations by using bulk operations
    score = sum(1 for position, g in enumerate(item) if g == main_target[position])
    return (item, float(score))


def crossover(parent_1: str, parent_2: str) -> tuple[str, str]:
    """
    Slice and combine two strings at a random point.
    >>> random.seed(42)
    >>> crossover("123456", "abcdef")
    ('12345f', 'abcde6')
    """
    random_slice = random.randint(0, len(parent_1) - 1)
    child_1 = parent_1[:random_slice] + parent_2[random_slice:]
    child_2 = parent_2[:random_slice] + parent_1[random_slice:]
    return (child_1, child_2)


def mutate(child: str, genes: list[str]) -> str:
    """
    Mutate a random gene of a child with another one from the list.
    >>> random.seed(123)
    >>> mutate("123456", list("ABCDEF"))
    '12345A'
    """
    # G1: Store repeated expressions in variables
    child_length = len(child)
    
    # G4: Use short-circuit evaluation
    if random.uniform(0, 1) >= MUTATION_PROBABILITY:
        return child
    
    # G6: Avoid redundant list conversion
    child_list = list(child)
    child_list[random.randint(0, child_length - 1)] = random.choice(genes)
    return "".join(child_list)


# Select, crossover and mutate a new population.
def select(
    parent_1: tuple[str, float],
    population_score: list[tuple[str, float]],
    genes: list[str],
) -> list[str]:
    """
    Select the second parent and generate new population

    >>> random.seed(42)
    >>> parent_1 = ("123456", 8.0)
    >>> population_score = [("abcdef", 4.0), ("ghijkl", 5.0), ("mnopqr", 7.0)]
    >>> genes = list("ABCDEF")
    >>> child_n = int(min(parent_1[1] + 1, 10))
    >>> population = []
    >>> for _ in range(child_n):
    ...     parent_2 = population_score[random.randrange(len(population_score))][0]
    ...     child_1, child_2 = crossover(parent_1[0], parent_2)
    ...     population.extend((mutate(child_1, genes), mutate(child_2, genes)))
    >>> len(population) == (int(parent_1[1]) + 1) * 2
    True
    """
    # G1: Store repeated expressions in variables
    parent_1_str = parent_1[0]
    parent_1_score = parent_1[1]
    
    # G3: Loop optimization - store end condition
    child_n = min(int(parent_1_score * 100) + 1, 10)
    
    # G7: Use bulk operations for better performance - generate all children at once
    # G7: Bulk operations - create all parent pairs first, then process in bulk
    parent_pairs = [
        (parent_1_str, population_score[random.randint(0, min(N_SELECTED, len(population_score) - 1))][0])
        for _ in range(child_n)
    ]
    
    # G7: Bulk operations - process all pairs and collect children
    children = []
    for parent_1_pair, parent_2_pair in parent_pairs:
        child_1, child_2 = crossover(parent_1_pair, parent_2_pair)
        # G7: Bulk append operations - avoid individual append calls
        children.extend([mutate(child_1, genes), mutate(child_2, genes)])
    
    return children


def basic(target: str, genes: list[str], debug: bool = True) -> tuple[int, int, str]:
    """
    Verify that the target contains no genes besides the ones inside genes variable.

    >>> from string import ascii_lowercase
    >>> basic("doctest", ascii_lowercase, debug=False)[2]
    'doctest'
    >>> genes = list(ascii_lowercase)
    >>> genes.remove("e")
    >>> basic("test", genes)
    Traceback (most recent call last):
        ...
    ValueError: ['e'] is not in genes list, evolution cannot converge
    >>> genes.remove("s")
    >>> basic("test", genes)
    Traceback (most recent call last):
        ...
    ValueError: ['e', 's'] is not in genes list, evolution cannot converge
    >>> genes.remove("t")
    >>> basic("test", genes)
    Traceback (most recent call last):
        ...
    ValueError: ['e', 's', 't'] is not in genes list, evolution cannot converge
    """

    # Verify if N_POPULATION is bigger than N_SELECTED
    if N_POPULATION < N_SELECTED:
        msg = f"{N_POPULATION} must be bigger than {N_SELECTED}"
        raise ValueError(msg)
    # Verify that the target contains no genes besides the ones inside genes variable.
    not_in_genes_list = sorted({c for c in target if c not in genes})
    if not_in_genes_list:
        msg = f"{not_in_genes_list} is not in genes list, evolution cannot converge"
        raise ValueError(msg)

    # G1: Store repeated expressions in variables
    target_length = len(target)
    
    # G7: Use bulk operations for population generation - single list comprehension
    population = ["".join(random.choice(genes) for _ in range(target_length)) for _ in range(N_POPULATION)]

    # G1: Store repeated expressions in variables
    generation, total_population = 0, 0
    population_size = len(population)

    # G3: Loop optimization with early termination
    while True:
        generation += 1
        total_population += population_size

        # G12: Enable multithreading for evaluation
        def evaluate_population():
            with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
                futures = {executor.submit(evaluate, item, target): item for item in population}
                return [future.result() for future in concurrent.futures.as_completed(futures)]
        
        # G7: Use bulk operations for evaluation
        population_score = evaluate_population()

        # G1: Store repeated expressions in variables
        population_score = sorted(population_score, key=lambda x: x[1], reverse=True)
        best_individual = population_score[0]
        best_string = best_individual[0]
        best_score = best_individual[1]

        # G4: Use short-circuit evaluation for early termination
        if best_string == target:
            return (generation, total_population, best_string)

        # G6: Avoid redundant computations
        if debug and generation % 10 == 0:
            print(
                f"\nGeneration: {generation}"
                f"\nTotal Population: {total_population}"
                f"\nBest score: {best_score}"
                f"\nBest string: {best_string}"
            )

        # G1: Use pre-computed constant
        population_best = population[:POPULATION_THIRD]
        population.clear()
        population.extend(population_best)
        
        # G6: Avoid redundant division
        target_length_inv = 1.0 / target_length
        population_score = [
            (item, score * target_length_inv) for item, score in population_score
        ]

        # G3: Loop optimization with early termination
        # G7: Use bulk operations for selection - collect all children at once
        new_population = []
        for i in range(N_SELECTED):
            # G7: Bulk operations - extend with all children from select at once
            new_population.extend(select(population_score[i], population_score, genes))
            # G3: Early termination condition
            if len(new_population) > N_POPULATION:
                break
        
        # G7: Bulk operations - extend population with all new children at once
        population.extend(new_population)
        population_size = len(population)


if __name__ == "__main__":
    target_str = (
        "This is a genetic algorithm to evaluate, combine, evolve, and mutate a string!"
    )
    genes_list = list(
        " ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklm"
        "nopqrstuvwxyz.,;!?+-*#@^'èéòà€ù=)(&%$£/\\"
    )
    generation, population, target = basic(target_str, genes_list)
    print(
        f"\nGeneration: {generation}\nTotal Population: {population}\nTarget: {target}"
    )
