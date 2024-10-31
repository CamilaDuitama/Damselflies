import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import numpy as np
import os
import argparse

# Fixed file paths
RESEQ_FILE = "Zenodo/reseq_coverage_norepeat_500_window.bed"
NANO_FILE = "Zenodo/Ifem_nano_coverage_norepeat_500_window.bed"
POPMAP_FILE = "Zenodo/SwD_popmap"

# Fig. 3a: Unitigs mapped to reference
def plot_unitigs(blast_file, assembly_type, comparison):
    print(f"Processing unitigs mapped to {assembly_type} reference for {comparison}...")
    blast_out = pd.read_csv(blast_file, sep="\t", header=None, usecols=[1,8])
    blast_out.columns = ["contig", "start"]
    
    unloc_2 = blast_out[blast_out["contig"] == "SUPER_13_unloc_2_RagTag"]
    unloc_2["comp"] = comparison
    
    percentage = len(unloc_2) / len(blast_out) * 100
    print(f"Percentage of significant unitigs mapping to SUPER_13_unloc_2_RagTag: {percentage:.2f}%")
    
    fig = px.histogram(unloc_2, x=unloc_2["start"]/1000000, nbins=36, 
                       labels={"x": "Position (Mb)", "y": "Unitig count"},
                       title=f"Unitigs mapped to {assembly_type} reference for {comparison}")
    fig.update_layout(showlegend=False, xaxis_range=[3.496, 3.78], yaxis_range=[0, 4500])
    
    return fig

# Fig. 3c: Read-depth coverage
def plot_read_depth():
    print("Processing read-depth coverage...")
    cov = pd.read_csv(RESEQ_FILE, sep="\t", header=None)
    cov.columns = ["contig", "start", "end", "coverage"]
    
    popmap = pd.read_csv(POPMAP_FILE, sep="\t")
    
    cov = cov[cov["contig"].isin(["1776_1", "28_1"])]
    
    cov["sample"] = np.tile(popmap["ind"].values, len(cov) // len(popmap))
    cov["morph"] = np.tile(popmap["pop"].values, len(cov) // len(popmap))
    cov["midpoint"] = (cov["end"] + cov["start"]) / 2
    
    S1 = pd.DataFrame({"sample": popmap["ind"]})
    S1["avg"] = S1["sample"].map(cov[cov["contig"] == "28_1"].groupby("sample")["coverage"].mean())
    
    cov["relcov"] = cov.apply(lambda row: row["coverage"] / S1.loc[S1["sample"] == row["sample"], "avg"].values[0], axis=1)
    
    ncov = pd.read_csv(NANO_FILE, sep="\t", header=None)
    ncov.columns = ["contig", "start", "end", "coverage"]
    
    ncov = ncov[ncov["end"] < 15000000]
    ncov["morph"] = np.repeat(["A", "I", "O"], len(ncov) // 3)
    ncov["midpoint"] = (ncov["end"] + ncov["start"]) / 2
    
    S2 = pd.DataFrame({"sample": ["A", "I", "O"]})
    S2["avg"] = S2["sample"].map(ncov[ncov["contig"] == "28_1"].groupby("morph")["coverage"].mean())
    
    ncov["relcov"] = ncov.apply(lambda row: row["coverage"] / S2.loc[S2["sample"] == row["morph"], "avg"].values[0], axis=1)
    
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.02)
    
    colors = px.colors.qualitative.Set1[:3]
    
    for i, morph in enumerate(["A", "I", "O"]):
        cov_morph = cov[(cov["contig"] == "1776_1") & (cov["morph"] == morph)]
        ncov_morph = ncov[(ncov["contig"] == "1776_1") & (ncov["morph"] == morph)]
        
        fig.add_trace(
            go.Scatter(x=(cov_morph["midpoint"].max() - cov_morph["midpoint"])/1000000, y=cov_morph["relcov"],
                       mode='lines', name=f'{morph} (short read)', line=dict(color=colors[i])),
            row=i+1, col=1
        )
        fig.add_trace(
            go.Scatter(x=(ncov_morph["midpoint"].max() - ncov_morph["midpoint"])/1000000, y=ncov_morph["relcov"],
                       mode='lines', name=f'{morph} (long read)', line=dict(color='black', dash='dash')),
            row=i+1, col=1
        )
    
    fig.update_layout(height=600, width=800, title_text="Read-depth coverage")
    fig.update_xaxes(title_text="Window midpoint on unlocalized scaffold 2 (Chr 13) (Mb)")
    fig.update_yaxes(title_text="Read depth relative to background", range=[0, 5])
    
    return fig

# Main execution
def main(blast_file):
    # Create docs directory if it doesn't exist
    os.makedirs("docs", exist_ok=True)
    
    # Extract information from the blast file path
    path_parts = blast_file.split('/')
    filter_params = path_parts[1]
    assembly_type = path_parts[2].split('_')[1]
    comparison = path_parts[3].split('_')[2].split('.')[0]
    
    # Fig. 3a
    fig_3a = plot_unitigs(blast_file, assembly_type, comparison)
    output_file_3a = f"docs/fig_3a_{filter_params}_{assembly_type}_{comparison}.html"
    fig_3a.write_html(output_file_3a)
    print(f"Figure 3a has been generated and saved as {output_file_3a}")
    
    # Fig. 3c
    fig_3c = plot_read_depth()
    output_file_3c = f"docs/fig_3c_{filter_params}_{assembly_type}_{comparison}.html"
    fig_3c.write_html(output_file_3c)
    print(f"Figure 3c has been generated and saved as {output_file_3c}")
    
    print("All plots have been generated and saved as HTML files in the docs folder.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate plots for unitig mapping and read depth coverage.")
    parser.add_argument("blast_file", help="Path to the BLAST-like output file")
    args = parser.parse_args()

    main(args.blast_file)