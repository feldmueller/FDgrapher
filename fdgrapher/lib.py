from gensim.models import Word2Vec
import pickle
from datetime import datetime
from sklearn import cluster
import pandas as pd
import numpy as np
from scipy.spatial import distance
from collections import Counter


def load_WEMs(meta_path_wem):
    meta_model = {}
    for meta, path in meta_path_wem:
        meta_model[meta] = Word2Vec.load(path).wv
    return meta_model


def load_kMeans(meta_path_km):
    meta_km = {}
    for meta, path in meta_path_km:
        with open(path, "rb") as f:
            meta_km[meta] = pickle.load(f)
    return meta_km


def train_kMeans(meta_model, meta_path_wem, k):
    meta_km = {}
    for meta, model in meta_model.items():
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print(f"Starting to cluster model {meta}: {current_time}")

        try:
            embeddings = model[model.key_to_index]
        except AttributeError:
            embeddings = model[model.index2word]

        n_clusters = int(len(embeddings) * k) if type(k) == float else k

        kmeans = cluster.KMeans(n_clusters=n_clusters)
        kmeans.fit(embeddings)
        meta_km[meta] = kmeans

        with open(meta_path_wem[meta] + ".kmeans.pkl", "wb") as f:
            pickle.dump(kmeans, f)

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print(f"Finished clustering model {meta}: {current_time}")

    return meta_km


def sort_dict(dic, reverse=True):
    return dict(sorted(dic.items(), key=lambda item: item[1], reverse=reverse))


def get_most_similar_clusters(wv, cluster_centroid):
    cluster_cosim = {}
    for cluster, centroid in cluster_centroid.items():
        cluster_cosim[cluster] = 1 - distance.cosine(centroid, wv)

    return sort_dict(cluster_cosim)


def get_t_centroid_clusters(
    meta,
    meta_model,
    meta_word_cluster,
    meta_cluster_words,
    meta_cluster_centroid,
    meta_cluster_label,
    seed_words,
    t_sim,
    neighbors,
    cutoff,
    shorthand=None,
    write_to_disc=True,
):
    def calculate_t_centroid(seedwords, meta, neighbors, cutoff):

        t_clusters = []

        for seedword in seedwords:
            if seedword not in meta_model[meta]:
                print(f"The seed word {seedword} is not contained in model {meta}!")
                continue
            try:
                topn = len(meta_model[meta].key_to_index)
            except AttributeError:
                topn = len(meta_model[meta].index2word)
            t_clusters.extend(
                [
                    meta_word_cluster[meta][word_dist[0]]
                    for word_dist in meta_model[meta].most_similar(
                        seedword,
                        topn=int(topn / neighbors),
                    )
                ]
            )

        cluster_freq = sort_dict(Counter(t_clusters))
        for c, freq in cluster_freq.items():
            cluster_freq[c] = freq / len(meta_cluster_words[meta][c])

        t_clusters = {c for c, freq in cluster_freq.items() if freq >= cutoff}

        return sort_dict(cluster_freq), t_clusters

    t_clusters = calculate_t_centroid(seed_words, meta, neighbors, cutoff)[1]
    t_centroid = np.mean([meta_cluster_centroid[meta][c] for c in t_clusters], axis=0)
    t_clusters_expanded = {
        c
        for c, sim in get_most_similar_clusters(
            t_centroid, meta_cluster_centroid[meta]
        ).items()
        if sim > t_sim
    } | t_clusters
    t_centroid = np.mean(
        [meta_cluster_centroid[meta][c] for c in t_clusters_expanded], axis=0
    )

    cluster_sim_t_centroid = pd.DataFrame(
        {
            "label": [
                (
                    shorthand[meta] + "_" + c
                    if shorthand and meta in shorthand
                    else meta + "_" + c
                )
                for c in meta_cluster_label[meta].values()
            ],
            "cluster": meta_cluster_label[meta].keys(),
            "n_words": [len(words) for words in meta_cluster_words[meta].values()],
        }
    )

    cluster_sim_t_centroid["thematic_similarity"] = cluster_sim_t_centroid.apply(
        lambda row: 1
        - distance.cosine(meta_cluster_centroid[meta][row.cluster], t_centroid),
        axis=1,
    )
    cluster_sim_t_centroid = cluster_sim_t_centroid.sort_values(
        "thematic_similarity", ascending=False
    )
    if write_to_disc:
        cluster_sim_t_centroid.to_csv(
            f"output/t_clusters_{meta}.csv",
            index=False,
        )

    return (t_clusters_expanded, t_centroid, cluster_sim_t_centroid)
