from collections import OrderedDict

class GeneAnalyzer:
    def __init__(self, anchor_seq_file, fasta_file):
        self.fasta_file = fasta_file
        self.anchor_seq_file = anchor_seq_file
        self.gene_data = OrderedDict()

    def anchor_gene(self, anchor_seq_file):
        gene_dict = OrderedDict()
        with open(anchor_seq_file, 'r') as file:
            combined_ids = None
            start_end = None

            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    start_end_raw = line.split(' # ')
                    start_end = start_end_raw[1] + ':' + start_end_raw[2]
                    combined_ids = line.split(' ')[0].strip('>')
                    gene_dict[combined_ids] = start_end

        return gene_dict

    def parse_fasta(self):
        with open(self.fasta_file, 'r') as file:
            current_id = None
            sequence = ""

            for line in file:
                line = line.strip()

                if line.startswith('>'):
                    if current_id and len(sequence) > 99:  
                        self.gene_data[current_id] = [raw_id, sequence]
                    start_end_raw = line.split(' # ')
                    current_id = start_end_raw[1] + ':' + start_end_raw[2]
                    raw_id = line.split(' ')[0].strip('>')
                    sequence = ""
                else:
                    sequence += line.replace('\n', '')
            if current_id and len(sequence) > 99:
                self.gene_data[current_id] = [raw_id, sequence]

    def analyze_genes(self, gene_dict, max_upstream=20, max_downstream=20):
        merged_ordered_dict = OrderedDict()

        gene_data_values_set = set(value[0] for value in self.gene_data.values())

        for ids, startend in gene_dict.items():
            primary_id = ids.rsplit('_',1)[0]
            if ids in gene_data_values_set:
                sequence = self.gene_data.get(startend, "")
                if startend in self.gene_data.keys():
                    index = list(self.gene_data.keys()).index(startend)

                    if sequence:
                        items = list(self.gene_data.items())
                        filtered_items = [item for item in items[max(0, index - max_upstream):index + max_downstream + 1]
                                          if item[1][0].rsplit('_',1)[0].startswith(primary_id)]
                        merged_ordered_dict.update(OrderedDict(filtered_items))

        return merged_ordered_dict

    def write_to_fasta(self, merged_ordered_dict, output_file):
        with open(output_file, 'a') as file:
            for gene_id, gene_info in merged_ordered_dict.items():
                file.write(f">{gene_info[0]}\n{gene_info[1]}\n")

if __name__ == '__main__':
    fasta_file_path = 'GCA_019232735.1_genomic_result.txt'
    anchor_seq_file_path = 'anchor_seq_2.fasta'

    gene_analyzer = GeneAnalyzer(anchor_seq_file_path, fasta_file_path)
    gene_dict = gene_analyzer.anchor_gene(anchor_seq_file_path)
    gene_analyzer.parse_fasta()

    analyzed_genes = gene_analyzer.analyze_genes(gene_dict)

    output_file_path = 'analyzed_genes.txt'
    gene_analyzer.write_to_fasta(analyzed_genes, output_file_path)
