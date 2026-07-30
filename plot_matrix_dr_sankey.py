import argparse
import pandas as pd
import re
import plotly.graph_objects as go
import os

# Strict Columns: DRB1 -> DRB6 -> DRB7 -> DRB5 -> DRB8 -> DRB4 -> DRB2 -> DRB3 -> DRB9
PATHS = {
    '01': ['DRB1', 'DRB6', 'DRB9'],
    '10': ['DRB1', 'DRB6', 'DRB9'],
    '15': ['DRB1', 'DRB6', 'DRB5', 'DRB9'],
    '16': ['DRB1', 'DRB6', 'DRB5', 'DRB9'],
    '03': ['DRB1', 'DRB2', 'DRB3', 'DRB9'],
    '11': ['DRB1', 'DRB2', 'DRB3', 'DRB9'],
    '12': ['DRB1', 'DRB2', 'DRB3', 'DRB9'],
    '13': ['DRB1', 'DRB2', 'DRB3', 'DRB9'],
    '14': ['DRB1', 'DRB2', 'DRB3', 'DRB9'],
    '04': ['DRB1', 'DRB7', 'DRB8', 'DRB4', 'DRB9'],
    '07': ['DRB1', 'DRB7', 'DRB8', 'DRB4', 'DRB9'],
    '09': ['DRB1', 'DRB7', 'DRB8', 'DRB4', 'DRB9'],
    '08': ['DRB1', 'DRB9']
}

# Hex palette for coloring
COLORS = {
    '01': '#1f77b4', '10': '#aec7e8', 
    '15': '#ff7f0e', '16': '#ffbb78', 
    '03': '#2ca02c', '11': '#98df8a', '12': '#8c564b', '13': '#c49c94', '14': '#e377c2',
    '04': '#d62728', '07': '#ff9896', '09': '#bcbd22',
    '08': '#9467bd'
}

def hex_to_rgba(hex_color, alpha=0.4):
    hex_color = hex_color.lstrip('#')
    if len(hex_color) == 6:
        r, g, b = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
        return f'rgba({r},{g},{b},{alpha})'
    return f'rgba(200,200,200,{alpha})'

def extract_dr_groups(cell_value):
    val = str(cell_value).strip()
    if val.lower() in ['nan', 'none', 'na', '']: return []
    matches = re.findall(r'\*(\d{2})', val)
    if not matches:
        matches = re.findall(r'(?:^|[, ])(\d{2}):', val)
    return matches

def main():
    parser = argparse.ArgumentParser(description="Generate a Straight-Pipe Segmented Sankey.")
    parser.add_argument("--input", required=True, help="Path to cds_summary_matrix (.tsv or .csv)")
    parser.add_argument("--col", required=True, help="Name of the DRB column")
    parser.add_argument("--outdir", default="Sankey_Output", help="Directory to save the images")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print(f"Loading data from {args.input}...")
    
    sep = '\t' if args.input.endswith('.tsv') else ','
    df = pd.read_csv(args.input, sep=sep)
    
    if args.col not in df.columns:
        print(f"Error: Column '{args.col}' not found.")
        return

    all_haplotypes = []
    for val in df[args.col]:
        all_haplotypes.extend(extract_dr_groups(val))
        
    valid_haplotypes = [g for g in all_haplotypes if g in PATHS]
    group_counts = pd.Series(valid_haplotypes).value_counts()
    total_assemblies = sum(group_counts.values)

    print(f"Extracted {total_assemblies} valid DRB1 assembled haplotypes.")

    # 1. SETUP MANUAL COLUMNS
    COLUMN_ORDER = ['Total', 'DRB1', 'DRB6', 'DRB7', 'DRB5', 'DRB8', 'DRB4', 'DRB2', 'DRB3', 'DRB9']
    col_dict = {col: [] for col in COLUMN_ORDER}

    for grp, count in group_counts.items():
        color = COLORS.get(grp, '#888888')
        for gene in PATHS[grp]:
            node_id = f"{grp}_{gene}"
            col_dict[gene].append( (node_id, count, grp, color) )

    # 2. CALCULATE STRICT HORIZONTAL LANES FOR EACH DR LINEAGE
  
    sorted_groups = sorted(group_counts.keys())
    lane_centers = {}
    
    
    current_y = 0.01 
    usable_height = 0.98
    
    # Calculate where on the Y-axis every lineage belongs based on its size
    for grp in sorted_groups:
        count = group_counts[grp]
        proportion = count / total_assemblies
        lane_height = proportion * usable_height
        
        # The mathematical center of this specific DR lineage's horizontal lane
        lane_centers[grp] = current_y + (lane_height / 2.0)
        
        # 0.5% pad between groups so the pipes don't bleed together visually
        current_y += lane_height + 0.005 

    # 3. ASSIGN X AND Y COORDINATES TO ALL NODES
    node_x, node_y, node_labels, node_colors, node_indices = [], [], [], [], {}
    current_index = 0
    
    for col_idx, col_name in enumerate(COLUMN_ORDER):
        x_pos = 0.01 + (col_idx / (len(COLUMN_ORDER) - 1)) * 0.98
        
        if col_name == 'Total':
            node_indices['Total Assemblies'] = current_index
            node_x.append(x_pos)
            node_y.append(0.5) # Master node stays centered
            node_labels.append('')
            node_colors.append('#333333')
            current_index += 1
            continue
            
        for node_id, val, grp, color in col_dict[col_name]:
            node_indices[node_id] = current_index
            node_x.append(x_pos)
            
            # Lock the Y coordinate to its specific lineage lane
            node_y.append(lane_centers[grp])
            
            label = f"DRB1*{grp}" if col_name == 'DRB1' else col_name
            node_labels.append(label)
            node_colors.append(color)
            current_index += 1

    # 4. BUILD THE FLOW LINKS
    source_indices, target_indices, link_values, link_colors = [], [], [], []
    
    for grp, count in group_counts.items():
        source_indices.append(node_indices['Total Assemblies'])
        target_indices.append(node_indices[f"{grp}_DRB1"])
        link_values.append(count)
        link_colors.append(hex_to_rgba(COLORS.get(grp, '#888888'), 0.5))

        path = PATHS[grp]
        for i in range(len(path) - 1):
            source_indices.append(node_indices[f"{grp}_{path[i]}"])
            target_indices.append(node_indices[f"{grp}_{path[i+1]}"])
            link_values.append(count)
            link_colors.append(hex_to_rgba(COLORS.get(grp, '#888888'), 0.4))

    print("making straight-pipe diagram")
    
    fig = go.Figure(data=[go.Sankey(
        arrangement = "fixed", 
        node = dict(
          pad = 0, 
          thickness = 25,
          line = dict(color = "black", width = 0.5),
          label = node_labels,
          color = node_colors,
          x = node_x,
          y = node_y
        ),
        link = dict(
          source = source_indices,
          target = target_indices,
          value = link_values,
          color = link_colors
        )
    )])

    fig.update_layout(
        title_text="<b>Segmented HLA-DR Sequential Haplotype Flow</b><br><sup>Total → DRB1 Lineage → Ordered Genomic Architecture (Strict Horizontal Flow)</sup>",
        font_size=13,
        height=850,
        width=1300 
    )

    html_out = os.path.join(args.outdir, "Segmented_DR_Sankey.html")
    fig.write_html(html_out)
    print(f"Interactive HTML saved to {html_out}")
    
    try:
        png_out = os.path.join(args.outdir, "Segmented_DR_Sankey.png")
        fig.write_image(png_out, scale=3) 
        print(f"High-Res PNG saved to {png_out}")
        
        pdf_out = os.path.join(args.outdir, "Segmented_DR_Sankey.pdf")
        fig.write_image(pdf_out)
        print(f"Vector PDF saved to {pdf_out}")
    except Exception:
        print("\n Static image export skipped. Open the HTML file to view or save!")

if __name__ == "__main__":
    main()
