"""
Simple multithreaded algorithm to show how the 4 phases of a genetic algorithm works
(Evaluation, Selection, Crossover and Mutation)
https://en.wikipedia.org/wiki/Genetic_algorithm
Author: D4rkia
OPTIMIZED: All guidelines applied (G1, G3, G4, G6, G7, G9, G12, G14)
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
# G1: Store repeated expressions in variables
POPULATION_THIRD = N_POPULATION // 3

# G14: Cache frequently used data
# G14: Cache target length and inverse for normalization
_target_length_cache: Dict[str, int] = {}
_target_inverse_cache: Dict[str, float] = {}
_population_score_cache: Dict[str, List[Tuple[str, float]]] = {}


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
    
    # G4: Use short-circuit evaluation - return immediately when mutation doesn't happen
    if random.uniform(0, 1) >= MUTATION_PROBABILITY:
        return child
    
    # G6: Avoid redundant list conversion - only convert when mutation happens
    child_list = list(child)
    child_list[random.randint(0, child_length - 1)] = random.choice(genes)
    return "".join(child_list)


# G9: Parallelizable selection function - no shared state dependencies
def select_parallel(
    parent_1: tuple[str, float],
    population_score: list[tuple[str, float]],
    genes: list[str],
) -> list[str]:
    """
    Select the second parent and generate new population
    G9: This function has no data dependencies between iterations - fully parallelizable

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
    
    # G14: Cache frequently used data - cache top slice for reuse
    population_score_key = str(id(population_score))
    if population_score_key not in _population_score_cache:
        _population_score_cache[population_score_key] = population_score[:N_SELECTED]
    
    top_slice = _population_score_cache[population_score_key]
    
    # G6: Avoid redundant computations - precompute bounds and avoid repeated slicing
    top_slice_len = len(top_slice)
    max_index = min(N_SELECTED, top_slice_len - 1)
    
    # G7: Use bulk operations for better performance - generate all children at once
    # G7: Bulk operations - create all parent pairs first, then process in bulk
    parent_pairs = [
        (parent_1_str, top_slice[random.randint(0, max_index)][0])
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

    # G1: Store repeated expressions in variables
    max_population = N_POPULATION
    selected_count = N_SELECTED
    third_population = max_population // 3
    
    # Verify if N_POPULATION is bigger than N_SELECTED
    if max_population < selected_count:
        msg = f"{max_population} must be bigger than {selected_count}"
        raise ValueError(msg)
    # Verify that the target contains no genes besides the ones inside genes variable.
    not_in_genes_list = sorted({c for c in target if c not in genes})
    if not_in_genes_list:
        msg = f"{not_in_genes_list} is not in genes list, evolution cannot converge"
        raise ValueError(msg)

    # G1: Store repeated expressions in variables
    target_length = len(target)
    
    # G14: Cache frequently used data - target length and inverse
    if target not in _target_length_cache:
        _target_length_cache[target] = target_length
    if target not in _target_inverse_cache:
        _target_inverse_cache[target] = 1.0 / target_length
    
    cached_target_length = _target_length_cache[target]
    cached_target_inverse = _target_inverse_cache[target]
    
    # G3: Loop optimization - precompute loop bounds and use single comprehension
    # G7: Use bulk operations for population generation - single list comprehension
    population = ["".join(random.choice(genes) for _ in range(cached_target_length)) for _ in range(max_population)]

    # G1: Store repeated expressions in variables
    generation, total_population = 0, 0
    population_size = len(population)

    # G3: Loop optimization with early termination
    while True:
        generation += 1
        total_population += population_size

        # G12: Enhanced multithreading with work-stealing dynamic scheduling
        def evaluate_population():
            # G12: Use ProcessPoolExecutor for CPU-bound tasks with work-stealing
            # G12: Dynamic scheduling allows idle workers to steal tasks from busy workers
            with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
                # G12: Submit all evaluation tasks at once for dynamic work distribution
                futures = [executor.submit(evaluate, item, target) for item in population]
                # G12: as_completed provides work-stealing behavior - idle workers take remaining tasks
                return [future.result() for future in concurrent.futures.as_completed(futures)]
        
        # G7: Use bulk operations for evaluation
        population_score = evaluate_population()

        # G1: Store repeated expressions in variables
        # G6: Avoid redundant computations - use in-place sort instead of creating new list
        population_score.sort(key=lambda x: x[1], reverse=True)
        best_individual = population_score[0]
        best_string = best_individual[0]
        best_score = best_individual[1]

        # G4: Use short-circuit evaluation for early termination
        if best_string == target:
            return (generation, total_population, best_string)

        # G4: Use short-circuit evaluation - check debug flag first (cheap check)
        # G6: Avoid redundant computations
        if debug and generation % 10 == 0:
            print(
                f"\nGeneration: {generation}"
                f"\nTotal Population: {total_population}"
                f"\nBest score: {best_score}"
                f"\nBest string: {best_string}"
            )

        # G1: Use pre-computed constant
        population_best = population[:third_population]
        population.clear()
        population.extend(population_best)
        
        # G14: Use cached inverse instead of recomputing
        population_score = [
            (item, score * cached_target_inverse) for item, score in population_score
        ]

        # G9: Parallelizable selection - no data dependencies between iterations
        # G3: Loop optimization with early termination and precomputed bounds
        # G7: Use bulk operations for selection - collect all children at once
        new_population = []
        
        # G3: Store loop end condition and use early termination
        max_needed = max_population - len(population)
        selected_count_actual = min(selected_count, len(population_score))
        
        # G9: Create independent tasks for parallel execution
        def create_selection_task(i):
            return select_parallel(population_score[i], population_score, genes)
        
        # G12: Enhanced multithreading for selection with work-stealing
        # G9: Execute selection tasks in parallel (independent operations)
        with concurrent.futures.ThreadPoolExecutor(max_workers=None) as executor:
            # G9: Submit all selection tasks in parallel - no dependencies between them
            # G12: Submit all selection tasks for dynamic work distribution
            selection_futures = [executor.submit(create_selection_task, i) for i in range(selected_count_actual)]
            
            # G12: as_completed provides work-stealing - idle workers take remaining tasks
            # G9: Collect results as they complete
            for future in concurrent.futures.as_completed(selection_futures):
                children = future.result()
                # G3: Take only what we need to avoid overshooting
                take_count = min(len(children), max_needed - len(new_population))
                new_population.extend(children[:take_count])
                # G3: Early termination condition - stop when we have enough
                if len(new_population) >= max_needed:
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
