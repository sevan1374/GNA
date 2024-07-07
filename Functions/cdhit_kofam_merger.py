import os
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from collections import defaultdict

def parse_cdhit_clusters(file_path):
    clusters = defaultdict(list)
    current_cluster = None
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>Cluster'):
                current_cluster = line.strip().split()[1]
            else:
                sequence_id = line.split('>')[1].split('...')[0]
                clusters[current_cluster].append(sequence_id)
    return clusters

def parse_kofam(file_path):
    kofam_kos = defaultdict(list)
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            sequence_id = parts[0]
            kos = parts[1]
            kofam_kos[sequence_id].append(kos)
    return kofam_kos

def cluster_vectors(clusters, kofam_kos):
    cluster_vec = {}
    all_kos = set()
    for seqs in clusters.values():
        for seq in seqs:
            all_kos.update(kofam_kos[seq])
    all_kos = list(all_kos)
    kos_index = {ko: idx for idx, ko in enumerate(all_kos)}

    for cluster_id, seqs in clusters.items():
        vec = [0] * len(all_kos)
        for seq in seqs:
            for ko in kofam_kos[seq]:
                vec[kos_index[ko]] += 1
        cluster_vec[cluster_id] = vec
    
    return cluster_vec, all_kos

def merge_clusters(cluster_vec, clusters):
    vecs = np.array(list(cluster_vec.values()))
    sim = cosine_similarity(vecs)
    np.fill_diagonal(sim, 0)

    threshold = 0.5
    merged = {}
    for i in range(len(sim)):
        for j in range(i + 1, len(sim)):
            if sim[i][j] > threshold:
                merged[i] = merged.get(i, set()) | {j}
                merged[j] = merged.get(j, set()) | {i}

    new_clusters = {}
    visited = set()
    for i in range(len(sim)):
        if i not in visited:
            all_indices = {i} | merged.get(i, set())
            combined_seqs = []
            for idx in all_indices:
                combined_seqs.extend(clusters[str(idx)])
                visited.add(idx)
            new_clusters[i] = combined_seqs
    return new_clusters

def write_clusters(clusters, output_file):
    with open(output_file, 'w') as file:
        for i, seqs in clusters.items():
            file.write(f'>Cluster {i}\n')
            for seq in seqs:
                file.write(f'{seq}\n')

def run_cdhit_kofam_analysis(cdhit_path, kofam_path, base_output_name):
    network_files_folder = os.path.join(os.getcwd(), "network_files")
    os.makedirs(network_files_folder, exist_ok=True)

    output_path = os.path.join(network_files_folder, f'combined_functional_cluster_{base_output_name}.txt')

    clusters = parse_cdhit_clusters(cdhit_path)
    kofam_kos = parse_kofam(kofam_path)
    cluster_vec, all_kos = cluster_vectors(clusters, kofam_kos)
    new_clusters = merge_clusters(cluster_vec, clusters)
    write_clusters(new_clusters, output_path)
    
    return output_path
