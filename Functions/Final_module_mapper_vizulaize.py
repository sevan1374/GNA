import pandas as pd
import plotly.graph_objects as go
import math
import os

def visualize_summary(file_path, base_output_name):
    # Load the data
    data = pd.read_csv(file_path, sep='\t')

    # Filter data to only include rows with count greater than 2
    filtered_data = data[data['count'] > 2]

    # Function to create gene block entries for each group
    def create_gene_blocks(data):
        gene_blocks = []
        y_position = 0
        for idx, row in data.iterrows():
            bacteria_ids = row['Bacteria_IDs'].split(',')
            kegg_modules = row['KEGG_module'].split(',')
            community = row['community_name']
            for i, module in enumerate(kegg_modules):
                gene_blocks.append({
                    'x_start': i,
                    'x_end': i + 1,
                    'y_start': y_position,
                    'y_end': y_position + 1,
                    'module': module,
                    'community': community,
                    'bacteria_ids': ','.join(bacteria_ids)
                })
            y_position += 1  # Move to the next row for the next group
        return gene_blocks

    # Define a fixed set of colors
    fixed_colors = [
        'rgba(31, 119, 180, 0.8)',  # Blue
        'rgba(255, 127, 14, 0.8)',  # Orange
        'rgba(44, 160, 44, 0.8)',   # Green
        'rgba(214, 39, 40, 0.8)',   # Red
        'rgba(148, 103, 189, 0.8)'  # Purple
    ]

    # Create the subfolder if it does not exist
    subfolder = os.path.join(os.getcwd(), "visualization")
    os.makedirs(subfolder, exist_ok=True)

    # Function to create and save the plot to HTML
    def create_plot(data_chunk, filename):
        gene_blocks = create_gene_blocks(data_chunk)
        unique_communities = sorted(set(block['community'] for block in gene_blocks))
        colors = {community: fixed_colors[i % len(fixed_colors)] for i, community in enumerate(unique_communities)}

        fig = go.Figure()
        
        for block in gene_blocks:
            fillcolor = colors[block['community']]
            
            fig.add_trace(go.Scatter(
                x=[block['x_start'], block['x_end'], block['x_end'], block['x_start'], block['x_start']],
                y=[block['y_start'], block['y_start'], block['y_end'], block['y_end'], block['y_start']],
                fill='toself',
                fillcolor=fillcolor,
                line=dict(color='rgb(0,0,0)'),
                text=f"Module: {block['module']}<br>Community: {block['community']}<br>Bacteria: {block['bacteria_ids']}",
                hoverinfo='text',
                mode='lines'
            ))
            fig.add_trace(go.Scatter(
                x=[(block['x_start'] + block['x_end']) / 2],
                y=[(block['y_start'] + block['y_end']) / 2],
                text=block['module'],
                mode='text',
                showlegend=False,
                textfont=dict(color='black')  # Use black color for text
            ))

        fig.update_layout(
            title="KEGG KO and Cluster Visualization",
            xaxis_title="KEGG Modules/Clusters",
            yaxis=dict(
                tickmode='array',
                tickvals=[block['y_start'] + 0.5 for block in gene_blocks if block['x_start'] == 0],
                ticktext=[block['community'] for block in gene_blocks if block['x_start'] == 0]
            ),
            showlegend=False
        )
        
        # Save the plot to an HTML file in the visualization subfolder
        file_path = os.path.join(subfolder, filename)
        fig.write_html(file_path)
        print(f"Plot saved to {file_path}")

    # Split the filtered data into chunks of 20 lines and create plots
    chunk_size = 20
    num_chunks = math.ceil(len(filtered_data) / chunk_size)

    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(filtered_data))
        data_chunk = filtered_data.iloc[start_idx:end_idx]
        create_plot(data_chunk, f"genome_plot_part_{base_output_name}_{i + 1}.html")
