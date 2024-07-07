import os
from collections import defaultdict
import pandas as pd

def process_communities(community_file, combined_file, base_output_name):
    detailed_communities_folder = os.path.join(os.getcwd(), "detailed_communities")
    os.makedirs(detailed_communities_folder, exist_ok=True)

    with open(community_file, 'r') as f:
        communities = f.read().split('\n\n')

    for community in communities:
        lines = community.split('\n')
        community_name = lines[0].strip(':')
        ids = [line for line in lines[1:] if not line.startswith('>Cluster')]

        id_dict = {}
        for id in ids:
            id_base, id_index = id.rsplit('_', 1)
            id_dict[id_base] = int(id_index)

        lines_dict = {id: [] for id in id_dict}

        with open(combined_file, 'r') as f:
            for line in f:
                if line.startswith('>Cluster'):
                    continue
                id_full, _ = line.split('\t', 1)
                id_base, id_index = id_full.rsplit('_', 1)
                id_index = int(id_index)
                if id_base in id_dict and abs(id_dict[id_base] - id_index) <= 20:
                    lines_dict[id_base].append((id_index, line))

        output_file_path = os.path.join(detailed_communities_folder, f'{community_name}_output_leiden_{base_output_name}.txt')
        with open(output_file_path, 'w') as out:
            for id in ids:
                id_base = id.rsplit('_', 1)[0]
                lines_dict[id_base].sort()
                for _, line in lines_dict[id_base]:
                    out.write(line)

        if os.stat(output_file_path).st_size == 0:
            os.remove(output_file_path)

    writer = pd.ExcelWriter(os.path.join(detailed_communities_folder, f'detailed_communities_{base_output_name}.xlsx'), engine='xlsxwriter')

    for community in communities:
        lines = community.split('\n')
        community_name = lines[0].strip(':')
        ids = [line for line in lines[1:] if not line.startswith('>Cluster')]

        id_dict = {}
        for id in ids:
            id_base, id_index = id.rsplit('_', 1)
            id_dict[id_base] = int(id_index)

        lines_dict = {id: [] for id in id_dict}

        with open(combined_file, 'r') as f:
            for line in f:
                if line.startswith('>Cluster'):
                    continue
                id_full, _ = line.split('\t', 1)
                id_base, id_index = id_full.rsplit('_', 1)
                id_index = int(id_index)
                if id_base in id_dict and abs(id_dict[id_base] - id_index) <= 20:
                    lines_dict[id_base].append((id_index, line))

        ko_dict = defaultdict(list)

        for id, lines in lines_dict.items():
            for _, line in lines:
                id_full, ko_info = line.split('\t', 1)
                if ko_info.startswith('Na'):
                    continue
                ko = ko_info.strip('\n')
                ko_dict[ko].append(id_full)

        sorted_kos = sorted(ko_dict.items(), key=lambda item: len(item[1]), reverse=True)
        output_data = []
        for ko, ids in sorted_kos:
            ko1 = ko.split(',')[0].strip("'(").strip("'").strip("['")
            ko2 = ko.split(',')[1].strip("'").strip(")").strip("]'")
            output_data.append({'KO1': ko1, 'KO2': ko2, 'Count': len(ids), 'IDs': ",".join(ids)})
        
        output_df = pd.DataFrame(output_data)
        output_df.to_excel(writer, sheet_name=community_name[:31], index=False)  

    writer.close()
