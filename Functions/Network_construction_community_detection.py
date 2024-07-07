import os
import networkx as nx
from networkx.algorithms import bipartite
import community as community_louvain
import leidenalg
import igraph as ig

def sequence_file_reader(input_fasta):
    list_fasta = []
    with open(input_fasta, 'r') as i1:
        for line in i1:
            if '>' in line:
                line = line.strip('\n').split(' ')[0].strip('>')
                list_fasta.append(line)
    return list_fasta

def read_indices_file(indices_file):
    indices_dict = {}
    with open(indices_file, 'r') as file:
        for line in file:
            line = line.strip()
            if '>' in line:
                line = line.strip('\n').split(' ')[0].strip('>')
                sequence_id = line.rsplit('_', 1)[0]
                index = int(line.rsplit('_', 1)[1])
                if sequence_id not in indices_dict:
                    indices_dict[sequence_id] = []
                indices_dict[sequence_id].append(index)
    return indices_dict

def cluster_file_reader(cluster_file, list_fasta, indices_dict):
    cluster_dict = {}
    temp_dict = {}
    with open(cluster_file, 'r') as file:
        current_cluster = None
        current_cluster_members = []
        for line in file:
            line = line.strip()
            if line.startswith(">Cluster"):
                if current_cluster and len(current_cluster_members) > 9:
                    for member in current_cluster_members:
                        if member not in list_fasta:
                            sequence_id = member.split()[0]
                            cluster_id = sequence_id.rsplit('_', 1)[0]
                            index = int(sequence_id.rsplit('_', 1)[1])
                            if cluster_id in indices_dict:
                                for i in indices_dict[cluster_id]:
                                    if index - 10 <= i <= index + 10:
                                        if f"{cluster_id}_{i}" not in temp_dict:
                                            temp_dict[f"{cluster_id}_{i}"] = []
                                        temp_dict[f"{cluster_id}_{i}"].append(current_cluster)
                current_cluster = line
                current_cluster_members = []
            else:
                current_cluster_members.append(line)
        if current_cluster and len(current_cluster_members) > 9:
            for member in current_cluster_members:
                if member not in list_fasta:
                    sequence_id = member.split()[0]
                    cluster_id = sequence_id.rsplit('_', 1)[0]
                    index = int(sequence_id.rsplit('_', 1)[1])
                    if cluster_id in indices_dict:
                        for i in indices_dict[cluster_id]:
                            if index - 10 <= i <= index + 10:
                                if f"{cluster_id}_{i}" not in temp_dict:
                                    temp_dict[f"{cluster_id}_{i}"] = []
                                temp_dict[f"{cluster_id}_{i}"].append(current_cluster)
    for key, value in temp_dict.items():
        cluster_dict[key] = list(set(value))
    return cluster_dict

def create_bipartite_network(cluster_dict):
    B = nx.Graph()
    for k, v in cluster_dict.items():
        if len(v) > 2:
            B.add_node(k, bipartite=0)  
            for node in v:
                B.add_node(node, bipartite=1) 
                B.add_edge(k, node)
    return B

def write_edgelist(B, filename):
    with open(filename, 'w') as f:
        for u, v, data in B.edges(data=True):
            if B.nodes[u]['bipartite'] == 0:  
                f.write(f"{u} {v} {data}\n")
            else:  
                f.write(f"{v} {u} {data}\n")

def write_communities(partition, filename):
    with open(filename, 'w') as f:
        for com in set(partition.values()):
            list_nodes = [nodes for nodes in partition.keys() if partition[nodes] == com]
            f.write(f"community {com+1}:\n")
            for node in list_nodes:
                f.write(f"{node}\n")
            f.write("\n")

def run_network_analysis(anchor_file, final_output_path, base_output_name):
    network_files_folder = os.path.join(os.getcwd(), "network_files")
    os.makedirs(network_files_folder, exist_ok=True)

    cluster_representation_file = os.path.join(network_files_folder, f"Cluster_representation_{base_output_name}.txt")
    louvain_communities_file = os.path.join(network_files_folder, f"louvain_communities_{base_output_name}.txt")
    leiden_communities_file = os.path.join(network_files_folder, f"leiden_communities_{base_output_name}.txt")
    network_combined_file = os.path.join(network_files_folder, f"network_{base_output_name}.txt")

    indices_dict = read_indices_file(anchor_file)
    list_fasta = sequence_file_reader(anchor_file)
    cluster_dict = cluster_file_reader(final_output_path, list_fasta, indices_dict)

    with open(cluster_representation_file, 'w') as o1:
        for k in cluster_dict:
            o1.write(f"{k}\t{cluster_dict[k]}\n")

    B = create_bipartite_network(cluster_dict)

    assert bipartite.is_bipartite(B)

    # Louvain community detection
    partition_louvain = community_louvain.best_partition(B)
    write_communities(partition_louvain, louvain_communities_file)

    # Leiden community detection
    g = ig.Graph.TupleList(B.edges(), directed=False)
    partition_leiden = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition)
    with open(leiden_communities_file, 'w') as f:
        for com in range(len(partition_leiden)):
            f.write(f"community {com+1}:\n")
            for node in partition_leiden[com]:
                f.write(f"{g.vs[node]['name']}\n")
            f.write("\n")

    write_edgelist(B, network_combined_file)

    return cluster_representation_file, louvain_communities_file, leiden_communities_file, network_combined_file
