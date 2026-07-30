import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches

# ==========================================
# 0. SETUP COMMAND LINE ARGUMENTS
# ==========================================
parser = argparse.ArgumentParser(description="Generate LOH plots for a specific tolerance level.")
parser.add_argument("--tol", type=int, required=True, help="The tolerance level number (e.g., enter 2 for 2per500kb)")
args = parser.parse_args()

# Store the provided tolerance level
tol_level = args.tol

# ==========================================
# 1. DEFINE YOUR MARKERS & BOUNDARIES HERE
# ==========================================
gene_markers = {
    'HLA-G': 29794320,
    'HLA-A': 29910247,
    'MUC22': 31003992,
    'HLA-C': 31236526,
    'HLA-B': 31321649,
    'MICA': 31402581,
    'C4': 31982057,
    'HLA-DRA':32437894, 
    'HLA-DRB1': 32578775,
    'HLA-DQ': 32627241,
    'HLA-DP': 33043444
}

# The absolute boundaries of targeted region
first_snp_pos = 28524330
last_snp_pos = 33413682
target_length = last_snp_pos - first_snp_pos  # 4,889,352

# 2. Load the data dynamically based on the --tol flag
csv_file = f"LOH_Tol_{tol_level}per500kb_Summary.csv"
try:
    df = pd.read_csv(csv_file)
except FileNotFoundError:
    print(f"❌ Error: Could not find {csv_file}. Make sure it is in the same directory.")
    exit()

# 3. Clean the Data
df = df[df['LOH_Start_Position'] != 'ERROR'].copy()

# Exclude BSH sample from analysis
df = df[df['Sample'] != 'BSH']

df['LOH_Start_Position'] = pd.to_numeric(df['LOH_Start_Position'])
df['LOH_End_Position'] = pd.to_numeric(df['LOH_End_Position'])
df['LOH_Length'] = pd.to_numeric(df['LOH_Length'])

# Sort by Start Position (earliest first), then by Length (shortest first)
df = df.sort_values(by=['LOH_Start_Position', 'LOH_Length'], ascending=[True, True])
df = df.reset_index(drop=True)

# 4. Initialize the Plot
fig, ax = plt.subplots(figsize=(12, 8))

# 5. Draw the Range Bars
y_positions = range(len(df), 0, -1)

# Dynamic Color Assignment
bar_colors = ['darkblue' if length >= target_length else 'royalblue' for length in df['LOH_Length']]

ax.barh(y=y_positions, 
        width=df['LOH_Length'], 
        left=df['LOH_Start_Position'], 
        color=bar_colors,  
        edgecolor=bar_colors, 
        linewidth=0.5,
        height=0.6,
        zorder=2)

# Lock the bottom of the y-axis strictly to 0
ax.set_ylim(bottom=0)

# ==========================================
# 6. Draw the Vertical Gene & Boundary Markers
# ==========================================
y_min, y_max = ax.get_ylim()

# First and Last SNP Boundaries (Gray Dashed)
boundary_color = 'dimgray'

ax.axvline(x=first_snp_pos, color=boundary_color, linestyle='--', linewidth=1.5, alpha=0.8, zorder=3)
ax.text(x=first_snp_pos, y=y_max, s='First SNP', color=boundary_color, fontsize=10, fontweight='bold', 
        rotation=90, verticalalignment='bottom', horizontalalignment='center')

ax.axvline(x=last_snp_pos, color=boundary_color, linestyle='--', linewidth=1.5, alpha=0.8, zorder=3)
ax.text(x=last_snp_pos, y=y_max, s='Last SNP', color=boundary_color, fontsize=10, fontweight='bold', 
        rotation=90, verticalalignment='bottom', horizontalalignment='center')

# Existing Gene Markers (Red Dotted)
for gene, pos in gene_markers.items():
    ax.axvline(x=pos, color='red', linestyle=':', alpha=0.6, zorder=3)
    ax.text(x=pos, 
            y=y_max, 
            s=gene, 
            color='black', 
            fontsize=10, 
            rotation=90, 
            verticalalignment='bottom', 
            horizontalalignment='center')

# 7. Formatting & Aesthetics 
ax.set_title(f'Homozygous tracts of samples ({tol_level} SNPs per 500kb Tol)', fontsize=14, pad=60)
ax.set_xlabel('chr6 (Mb)', fontsize=12)
ax.set_ylabel('Sample Number', fontsize=12)

# --- Custom Legend ---
partial_patch = mpatches.Patch(color='royalblue', label='Partial Homozygous')
full_patch = mpatches.Patch(color='darkblue', label='Full Homozygous')
ax.legend(handles=[partial_patch, full_patch], 
          loc='lower left', 
          framealpha=0.9, 
          edgecolor='black')

# Force ticks every 500,000 base pairs (0.5 Mb)
ax.xaxis.set_major_locator(ticker.MultipleLocator(500000))

# Converts base pairs to Megabases
formatter = ticker.FuncFormatter(lambda x, pos: f'{x / 1_000_000:.2f}')
ax.xaxis.set_major_formatter(formatter)

ax.grid(axis='y', linestyle='-', alpha=0.2, zorder=1)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.set_xlim(28400000, 33600000)

# 8. Save and Display
plt.tight_layout()
output_filename = f"LOH_Tol_{tol_level}per500kb_Graph.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"Graph successfully saved as {output_filename}")

# plt.show() # Commented this out so it doesn't pause a bash loop
