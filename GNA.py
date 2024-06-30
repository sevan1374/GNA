import argparse
import os
from Functions.Gene_neighborhood_collector import GeneAnalyzer

if __name__ == '__main__':
    # Create the parser
    parser = argparse.ArgumentParser(description='Analyze genes')

    # Add the arguments
    parser.add_argument('-i', '--input', type=str, help='The input file')
    parser.add_argument('-a', '--anchor', type=str, required=True, help='The anchor sequence file')
    parser.add_argument('-o', '--output', type=str, required=True, help='The output file')
    parser.add_argument('-d', '--directory', type=str, help='The directory of input files')

    # Parse the arguments
    args = parser.parse_args()
    print(f"Parsed arguments: {args}")

    # Clear the output file before processing
    open(args.output, 'w').close()

    if args.directory:
        print(f"Processing directory: {args.directory}")
        # If a directory is provided, process all files in the directory
        for filename in os.listdir(args.directory):
            if filename.endswith('.txt'):
                fasta_file_path = os.path.join(args.directory, filename)
                print(f"Processing file: {fasta_file_path}")
                gene_analyzer = GeneAnalyzer(args.anchor, fasta_file_path)
                gene_dict = gene_analyzer.anchor_gene(args.anchor)
                gene_analyzer.parse_fasta()
                analyzed_genes = gene_analyzer.analyze_genes(gene_dict)
                gene_analyzer.write_to_fasta(analyzed_genes, args.output)
    elif args.input:
        print(f"Processing single input file: {args.input}")
        # If an input file is provided, process the single input file
        gene_analyzer = GeneAnalyzer(args.anchor, args.input)
        gene_dict = gene_analyzer.anchor_gene(args.anchor)
        gene_analyzer.parse_fasta()
        analyzed_genes = gene_analyzer.analyze_genes(gene_dict)
        gene_analyzer.write_to_fasta(analyzed_genes, args.output)
    else:
        print("Please provide either an input file (-i) or a directory of input files (-d).")
