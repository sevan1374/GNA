import os
import re

def parse_kofamscan_results(file_path):
    results = {}
    else_dic = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('*'):
                parts = re.split(r'\s+', line.strip())
                gene_name = parts[1]
                ko = parts[2]
                ko_definition = ' '.join(parts[6:])
                results[gene_name] = (ko, ko_definition)
            else:
                parts = re.split(r'\s+', line.strip())
                gene_name = parts[0]
                ko = parts[2]
                ko_definition = ' '.join(parts[6:])
                else_dic[gene_name] = (ko, ko_definition)
    return results, else_dic

def check_for_missing(raw_input, kofamscan_results):
    missed_ids = []
    results, else_dic = kofamscan_results
    with open(raw_input, 'r') as file_1:
        for line_1 in file_1:
            if '>' in line_1:
                line_1 = line_1.strip('\n').strip('>')
                if line_1 not in results.keys() and line_1 not in else_dic.keys():
                    missed_ids.append(line_1)
    return results, missed_ids

def parse_for_missed_genes(file_path, results):
    list_positive = set(results.keys())
    final_results = {}

    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#') and not line.startswith('*'):
                parts = re.split(r'\s+', line.strip())
                gene_name = parts[0]
                ko = parts[1]
                threshold = parts[2]
                score = parts[3]
                ko_definition = ' '.join(parts[5:])
                if gene_name not in list_positive and str(threshold) != '-':
                    pass_threshold = float((float(threshold) - float(score)) / float(threshold))
                    if pass_threshold <= 0.10:
                        if gene_name not in final_results:
                            final_results[gene_name] = (pass_threshold, ko, ko_definition)
                        else:
                            if pass_threshold < final_results[gene_name][0]:
                                final_results[gene_name] = (pass_threshold, ko, ko_definition)

    return final_results

def write_results_to_files(base_output_name, results, missed_ids, final_results):
    parsed_kofam_folder = os.path.join(os.getcwd(), "parsed_kofam")
    os.makedirs(parsed_kofam_folder, exist_ok=True)

    missed_ids_file = os.path.join(parsed_kofam_folder, f"Missed_ids_{base_output_name}.txt")
    results_file = os.path.join(parsed_kofam_folder, f"results_{base_output_name}.txt")
    missed_genes_file = os.path.join(parsed_kofam_folder, f"Missed_{base_output_name}.txt")
    combined_results_file = os.path.join(parsed_kofam_folder, f"results_{base_output_name}_10.txt")

    # Write results to files
    with open(results_file, 'w') as o2:
        for gene_name, (ko, ko_definition) in results.items():
            o2.write(f"{gene_name}\t{ko}\t{ko_definition}\n")

    # Write missed genes to a separate file
    with open(missed_genes_file, 'w') as o3:
        for gene_name, (pass_threshold, ko, ko_definition) in final_results.items():
            o3.write(f"{gene_name}\t{ko}\t{ko_definition}\n")

    # Combine results and missed genes into one file
    with open(combined_results_file, 'w') as combined_output:
        for gene_name in results:
            combined_output.write(f"{gene_name}\t{results[gene_name]}\n")

        for gene_name in final_results:
            if gene_name not in results:  # Ensure no duplicates
                combined_output.write(f"{gene_name}\t{final_results[gene_name][1:]}\n")

    return results_file, missed_genes_file, combined_results_file
