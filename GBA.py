import argparse
import os
from collections import defaultdict
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import networkx as nx
from networkx.algorithms import bipartite
import community as community_louvain
import leidenalg
import igraph as ig
from Functions.Kofam_parser import parse_kofamscan_results, check_for_missing, parse_for_missed_genes, write_results_to_files
from Functions.cdhit_kofam_merger import run_cdhit_kofam_analysis
from Functions.Network_construction_community_detection import run_network_analysis
from Functions.Cluster_representation_mapper import run_cluster_representation_mapper
from Functions.Final_module_mapper import run_final_module_mapper
from Functions.Final_module_mapper_vizulaize import visualize_summary
from Functions.General_community_mapper import process_communities

"""
Arguments:
-a: Anchor fasta sequences  containing the collected upstream and downstream genes of the our genes of interest.
-o: file name for gene of interst.
--kofamscan: kofamscam annotation file for the collected collected upstream and downstream fasta file.
--cdhit: cluster .clstr file for the collected collected upstream and downstream fasta file.
"""
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze genes')

    parser.add_argument('-a', '--anchor', type=str, required=True, help='The anchor sequence file')
    parser.add_argument('-o', '--output', type=str, required=True, help='The output file')
    parser.add_argument('--kofamscan', type=str, required=True, help='The KofamScan results file')
    parser.add_argument('--cdhit', type=str, required=True, help='The CD-HIT cluster file')

    args = parser.parse_args()
    print(f"Parsed arguments: {args}\n")
    open(args.output, 'w').close()

    kofamscan_file_path = args.kofamscan
    input_file_path = args.output  

    # Parse kofamscan results
    results, else_dic = parse_kofamscan_results(kofamscan_file_path)

    # Check for missing genes
    results, missed_ids = check_for_missing(input_file_path, (results, else_dic))

    # Parse for missed genes
    final_results = parse_for_missed_genes(kofamscan_file_path, results)

    # Generate output filenames based on the user-specified output name
    base_output_name = os.path.splitext(args.output)[0]
    
    # Write Kofam scan results to files
    results_file, missed_genes_file, combined_results_file = write_results_to_files(base_output_name, results, missed_ids, final_results)

    # Process the CD-HIT clusters and combine clusters with similar KOs.
    cdhit_path = args.cdhit
    kofam_path = combined_results_file
    final_output_path = run_cdhit_kofam_analysis(cdhit_path, kofam_path, base_output_name)

    print(f"Final clusters written to {final_output_path}")

    # Network construction and community detection
    cluster_representation_file, louvain_communities_file, leiden_communities_file, network_combined_file = run_network_analysis(
        args.anchor, final_output_path, base_output_name)

    print(f"Network and community files written: {network_combined_file}, {louvain_communities_file}, {leiden_communities_file}, {cluster_representation_file}")

    # Map protein cluster files
    mapped_cluster_file = run_cluster_representation_mapper(final_output_path, combined_results_file, base_output_name)
    print(f"Mapped clusters written to {mapped_cluster_file}")

    # Module finder in all identified communities
    summary_file = run_final_module_mapper(leiden_communities_file, mapped_cluster_file, base_output_name)
    print(f"Summary written to {summary_file}")

    # visualization of major modules
    visualize_summary(summary_file, base_output_name)

    # Gene content summary for identified communities
    process_communities(leiden_communities_file, mapped_cluster_file, base_output_name)
    print(f"Detailed community files written in 'detailed_communities' subfolder")

    if os.path.exists(args.output):
        os.remove(args.output)