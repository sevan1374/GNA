import os
from collections import defaultdict

def read_cluster_file(cluster_file_path):
    clusters = {}
    current_cluster = None
    with open(cluster_file_path, 'r') as file:
        for line in file:
            if line.startswith('>Cluster'):
                current_cluster = line.strip().split()[1]
                clusters[current_cluster] = []
            else:
                sequence_id = line
                clusters[current_cluster].append(sequence_id.strip('\n'))
    return clusters

def read_kofam_file(kofam_file_path):
    kofam_annotations = {}
    with open(kofam_file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            sequence_id = parts[0]
            annotation = parts[1]
            kofam_annotations[sequence_id] = annotation
    return kofam_annotations

def map_kofam_to_clusters(clusters, kofam_annotations):
    mapped_clusters = {}
    missing = 0
    for cluster_id, sequences in clusters.items():
        mapped_clusters[cluster_id] = []
        for seq in sequences:
            if seq in kofam_annotations:
                annotation = kofam_annotations[seq]
            else:
                annotation = 'Na'
                missing += 1
            mapped_clusters[cluster_id].append((seq, annotation))
    print(f"Missing Kofam annotations for {missing} sequences.")
    return mapped_clusters

def write_output(mapped_clusters, output_file_path):
    with open(output_file_path, 'w') as file:
        for cluster_id, seq_annotations in mapped_clusters.items():
            file.write(f'>Cluster {cluster_id}\n')
            for seq, annotation in seq_annotations:
                file.write(f'{seq}\t{annotation}\n')

def run_cluster_representation_mapper(cluster_file_path, kofam_file_path, base_output_name):
    network_files_folder = os.path.join(os.getcwd(), "network_files")
    os.makedirs(network_files_folder, exist_ok=True)

    output_file_path = os.path.join(network_files_folder, f"Mapped_combined_functional_cluster_{base_output_name}.txt")

    clusters = read_cluster_file(cluster_file_path)
    kofam_annotations = read_kofam_file(kofam_file_path)
    mapped_clusters = map_kofam_to_clusters(clusters, kofam_annotations)
    write_output(mapped_clusters, output_file_path)

    print(f"Output has been written to {output_file_path}")
    return output_file_path
