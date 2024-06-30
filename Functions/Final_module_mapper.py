import os
import re
from collections import defaultdict

# Function to parse community data
def parse_community_data(community_file):
    communities = defaultdict(lambda: {'bacteria': [], 'clusters': []})
    current_community = None
    with open(community_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("community"):
                current_community = line.split()[1]
            elif line.startswith(">Cluster"):
                communities[current_community]['clusters'].append(line.split()[1])
            elif line:
                communities[current_community]['bacteria'].append(line)
    return communities

# Function to parse cluster data
def parse_cluster_data(cluster_file):
    clusters = defaultdict(list)
    current_cluster = None
    with open(cluster_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_cluster = line.split()[1]
            elif line:
                parts = line.split('\t')
                genome_id = parts[0].rsplit('_', 1)[0]
                kegg_info = re.search(r"\('(\w+)',", parts[1]) if len(parts) > 1 else None
                kegg_id = kegg_info.group(1) if kegg_info else 'NA'
                clusters[current_cluster].append((genome_id, kegg_id))
    return clusters

# Function to map bacteria to KEGG IDs based on community clusters
def map_bacteria_to_kegg(communities, clusters):
    community_kegg_map = defaultdict(list)
    for community, data in communities.items():
        bacteria_ids = data['bacteria']
        community_clusters = data['clusters']
        
        for bacteria_id in bacteria_ids:
            genome_id = bacteria_id.rsplit('_', 1)[0]
            kegg_ids = set()
            
            for cluster in community_clusters:
                cluster_entries = clusters.get(cluster, [])
                for entry in cluster_entries:
                    if entry[0] == genome_id:
                        kegg_ids.add(entry[1] if entry[1] != 'NA' else f"Cluster {cluster}")
            
            kegg_list = ','.join(kegg_ids)
            community_kegg_map[community].append((bacteria_id, kegg_list))
    
    return community_kegg_map

# Function to find similar KEGG modules and group bacteria
def find_similar_kegg_groups(community_kegg_map, similarity_threshold=0.50):
    output_data = []
    for community, bacteria_kegg_list in community_kegg_map.items():
        bacteria_ids = [b[0] for b in bacteria_kegg_list]
        kegg_lists = [set(b[1].split(',')) for b in bacteria_kegg_list]
        
        groups = []
        used = set()
        
        for i in range(len(bacteria_ids)):
            if i in used:
                continue
            group = [bacteria_ids[i]]
            for j in range(i+1, len(bacteria_ids)):
                if j in used:
                    continue
                intersection = kegg_lists[i].intersection(kegg_lists[j])
                union = kegg_lists[i].union(kegg_lists[j])
                similarity = len(intersection) / len(union)
                if similarity >= similarity_threshold:
                    group.append(bacteria_ids[j])
                    used.add(j)
            used.add(i)
            groups.append(group)
        
        for group in groups:
            if len(group) > 1:
                common_keggs = set.intersection(*[kegg_lists[bacteria_ids.index(b)] for b in group])
                common_keggs_str = ','.join(common_keggs) if common_keggs else "No common KEGG modules"
                output_data.append(f"{','.join(group)}\t{len(group)}\t{common_keggs_str}\tcommunity {community}")
    
    return output_data

# Function to write output to file
def write_output(output_file, output_data):
    with open(output_file, 'w') as file:
        file.write("Bacteria_IDs\tcount\tKEGG_module\tcommunity_name\n")
        for line in output_data:
            file.write(line + "\n")

# Main function to process files
def run_final_module_mapper(community_file, cluster_file, base_output_name):
    visualization_folder = os.path.join(os.getcwd(), "visualization")
    os.makedirs(visualization_folder, exist_ok=True)

    output_file = os.path.join(visualization_folder, f"leiden_communities_{base_output_name}_summary_20.txt")

    communities = parse_community_data(community_file)
    clusters = parse_cluster_data(cluster_file)
    community_kegg_map = map_bacteria_to_kegg(communities, clusters)
    output_data = find_similar_kegg_groups(community_kegg_map)
    write_output(output_file, output_data)
    return output_file
