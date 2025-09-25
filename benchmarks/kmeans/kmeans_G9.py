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

def assign_clusters(data, centroids):
    distances_from_centroids = pairwise_distances(data, centroids, metric="euclidean")
    cluster_assignment = np.argmin(distances_from_centroids, axis=1)
    return cluster_assignment

# G9: reduce dependencies between iterations - centroids updated independently
def revise_centroids(data, k, cluster_assignment):
    # [G9 applied here: each cluster's centroid can be computed independently]
    new_centroids = []
    for i in range(k):
        member_data_points = data[cluster_assignment == i]
        centroid = member_data_points.mean(axis=0)
        new_centroids.append(centroid)
    new_centroids = np.array(new_centroids)
    return new_centroids

def compute_heterogeneity(data, k, centroids, cluster_assignment):
    heterogeneity = 0.0
    for i in range(k):
        member_data_points = data[cluster_assignment == i, :]
        if member_data_points.shape[0] > 0:
            distances = pairwise_distances(member_data_points, [centroids[i]], metric="euclidean")
            squared_distances = distances**2
            heterogeneity += np.sum(squared_distances)
    return heterogeneity

def kmeans(data, k, initial_centroids, maxiter=500, record_heterogeneity=None, verbose=False):
    centroids = initial_centroids[:]
    prev_cluster_assignment = None
    for itr in range(maxiter):
        cluster_assignment = assign_clusters(data, centroids)
        # centroids updated independently (conceptually parallelizable)
        centroids = revise_centroids(data, k, cluster_assignment)  # [G9 applied here]
        if prev_cluster_assignment is not None and (prev_cluster_assignment == cluster_assignment).all():
            break
        if record_heterogeneity is not None:
            score = compute_heterogeneity(data, k, centroids, cluster_assignment)
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
