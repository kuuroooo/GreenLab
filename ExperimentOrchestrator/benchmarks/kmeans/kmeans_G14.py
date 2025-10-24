import warnings
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.metrics import pairwise_distances

warnings.filterwarnings("ignore")

def get_initial_centroids(data, k, seed=None):
    rng = np.random.default_rng(seed)
    n = data.shape[0]
    rand_indices = rng.integers(0, n, k)
    centroids = data[rand_indices, :]
    return centroids

# G14: cache distances to avoid recomputation
def assign_clusters(data, centroids, cached_distances=None):
    if cached_distances is None:
        cached_distances = pairwise_distances(data, centroids, metric="euclidean")  # [G14 applied here]
    cluster_assignment = np.argmin(cached_distances, axis=1)
    return cluster_assignment, cached_distances

def revise_centroids(data, k, cluster_assignment):
    new_centroids = []
    for i in range(k):
        member_data_points = data[cluster_assignment == i]
        centroid = member_data_points.mean(axis=0)
        new_centroids.append(centroid)
    new_centroids = np.array(new_centroids)
    return new_centroids

def compute_heterogeneity(data, k, centroids, cluster_assignment, cached_distances=None):
    heterogeneity = 0.0
    for i in range(k):
        member_data_points = data[cluster_assignment == i, :]
        if member_data_points.shape[0] > 0:
            # reuse cached distances if available
            if cached_distances is not None:
                distances = cached_distances[cluster_assignment == i, i][:, None]
            else:
                distances = pairwise_distances(member_data_points, [centroids[i]], metric="euclidean")
            squared_distances = distances**2
            heterogeneity += np.sum(squared_distances)
    return heterogeneity

def kmeans(data, k, initial_centroids, maxiter=500, record_heterogeneity=None, verbose=False):
    centroids = initial_centroids[:]
    prev_cluster_assignment = None
    for itr in range(maxiter):
        cluster_assignment, cached_distances = assign_clusters(data, centroids)  # [G14 applied here]
        centroids = revise_centroids(data, k, cluster_assignment)
        if prev_cluster_assignment is not None and (prev_cluster_assignment == cluster_assignment).all():
            break
        if record_heterogeneity is not None:
            score = compute_heterogeneity(data, k, centroids, cluster_assignment, cached_distances)
            record_heterogeneity.append(score)
        prev_cluster_assignment = cluster_assignment[:]
    return centroids, cluster_assignment

if __name__ == "__main__":
    from sklearn import datasets
    X = datasets.load_iris().data[:, :3]  # first 3 features
    k = 3
    initial_centroids = get_initial_centroids(X, k, seed=0)
    heterogeneity = []
    centroids, cluster_assignment = kmeans(
        X, k, initial_centroids, maxiter=20, record_heterogeneity=heterogeneity, verbose=True
    )
    print("Final centroids:\n", centroids)
    print("Heterogeneity history:", heterogeneity)
