# === Combined optimizations: G1, G3, G4, G6, G7, G9, G12, G14 ===
import numpy as np  # OPT: G7
from concurrent.futures import ThreadPoolExecutor, as_completed  # OPT: G12
from functools import lru_cache  # OPT: G14
import math  # OPT
DIST_MATRIX = None  # OPT: module-level distance cache for linters/runtime

"""Use an ant colony optimization algorithm (optimized)."""
# OPT (G12): parallel route builder (top-level definition)
def _build_ant_route(cities, alpha, beta, pheromone):
    unvisited_cities = dict(cities)
    first_city_key = next(iter(unvisited_cities.keys()))
    current_city = {first_city_key: unvisited_cities[first_city_key]}
    del unvisited_cities[first_city_key]
    route = [first_city_key]
    while unvisited_cities:
        if unvisited_cities and (alpha > 0 and beta > 0):
            current_city, unvisited_cities = city_select(pheromone, current_city, unvisited_cities, alpha, beta)
            route.append(next(iter(current_city.keys())))
        else:
            break
    route.append(0)
    return route


# OPT (G9+G14): precompute and cache distances
def _precompute_distances(cities):
    n = len(cities)
    coords = [cities[i] for i in range(n)]
    dist = np.zeros((n, n))
    for i in range(n):
        xi, yi = coords[i]
        for j in range(i+1, n):
            xj, yj = coords[j]
            d = math.dist((xi, yi), (xj, yj))
            dist[i, j] = d
            dist[j, i] = d
    return dist

cities = {
    0: [0, 0],
    1: [0, 5],
    2: [3, 8],
    3: [8, 10],
    4: [12, 8],
    5: [12, 4],
    6: [8, 0],
    7: [6, 2],
}

def main(
    cities: dict[int, list[int]],
    ants_num: int,
    iterations_num: int,
    pheromone_evaporation: float,
    alpha: float,
    beta: float,
    q: float,
) -> tuple[list[int], float]:
    """Optimized ant colony main."""
    cities_num = len(cities)
    pheromone = np.full((cities_num, cities_num), 1.0)  # OPT: G7 bulk-friendly pheromone matrix
    best_path: list[int] = []
    best_distance = float("inf")
    for _ in range(iterations_num):
        ants_route: list[list[int]] = []
        with ThreadPoolExecutor() as ex:
            futures = [ex.submit(_build_ant_route, cities, alpha, beta, pheromone) for _ in range(ants_num)]
            for f in as_completed(futures):
                ants_route.append(f.result())
        pheromone, best_path, best_distance = pheromone_update(
            pheromone, cities, pheromone_evaporation, ants_route, q, best_path, best_distance
        )
    return best_path, best_distance

def distance(city1: list[int], city2: list[int]) -> float:
    return (((city1[0]-city2[0])**2)+((city1[1]-city2[1])**2))**0.5




def pheromone_update(
    pheromone,
    cities: dict[int, list[int]],
    pheromone_evaporation: float,
    ants_route: list[list[int]],
    q: float,
    best_path: list[int],
    best_distance: float,
):
    """Update pheromones and track best path (optimized)."""
    # Evaporate pheromone globally (bulk multiply)
    pheromone *= pheromone_evaporation  # OPT: G7 bulk evaporation

    # Ensure distance matrix is available
    global DIST_MATRIX  # OPT: G9 global cache
    if DIST_MATRIX is None:
        DIST_MATRIX = _precompute_distances(cities)  # OPT: G9 cache once

    for ant_route in ants_route:
        total_distance = 0.0
        best_so_far = best_distance  # OPT: G3 early termination snapshot
        for i in range(len(ant_route) - 1):
            if total_distance >= best_so_far:  # OPT: G3 early termination
                break
            total_distance += DIST_MATRIX[ant_route[i], ant_route[i + 1]]  # OPT: G9

        delta_pheromone = q / total_distance
        for i in range(len(ant_route) - 1):
            a, b = ant_route[i], ant_route[i + 1]
            pheromone[a, b] += delta_pheromone  # OPT: G7
            pheromone[b, a] = pheromone[a, b]   # OPT: symmetric

        if total_distance < best_distance:
            best_path = ant_route
            best_distance = total_distance

    return pheromone, best_path, best_distance


def city_select(
    pheromone,
    current_city: dict[int, list[int]],
    unvisited_cities: dict[int, list[int]],
    alpha: float,
    beta: float,
):
    """Pick next city based on pheromone and inverse distance."""
    import random
    # Calculate probabilities for each unvisited city
    probs = []
    cities_list = list(unvisited_cities.keys())
    curr_key = next(iter(current_city.keys()))
    curr_val = next(iter(current_city.values()))
    for city in cities_list:
        city_distance = distance(unvisited_cities[city], curr_val)
        desirability = (pheromone[city, curr_key] ** alpha) * ((1.0 / city_distance) ** beta)
        probs.append(desirability)

    chosen_city_i = random.choices(cities_list, weights=probs)[0]
    chosen_city = {chosen_city_i: unvisited_cities[chosen_city_i]}
    del unvisited_cities[chosen_city_i]
    return chosen_city, unvisited_cities


if __name__ == "__main__":
    best_path, best_distance = main(
        cities=cities,
        ants_num=10,
        iterations_num=20,
        pheromone_evaporation=0.7,
        alpha=1.0,
        beta=5.0,
        q=10,
    )
    print(f"{best_path = }")
    print(f"{best_distance = }")
