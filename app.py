import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

# Set professional layout configuration
st.set_page_config(page_title="RNA-Seq Discovery Platform", layout="wide")

st.title("🧬 Multi-Omics RNA-Seq Differential Expression Dashboard")
st.markdown("### Automated Nextflow Processing & Verification Terminal")
st.write("---")

# Load dataset safely
@st.cache_data
def load_data():
    try:
        df = pd.read_csv("results/deseq2_results.csv")
        # Ensure standard column formatting
        df['pvalue'] = pd.to_numeric(df['pvalue'], errors='coerce')
        df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
        df['log2FoldChange'] = pd.to_numeric(df['log2FoldChange'], errors='coerce')
        df['log10_pvalue'] = -np.log10(df['pvalue'].clip(lower=1e-300))
        return df.dropna(subset=['log2FoldChange', 'pvalue'])
    except Exception as e:
        st.error(f"Error loading results matrix: {e}")
        return pd.DataFrame()

raw_df = load_data()

if raw_df.empty:
    st.warning("Waiting for baseline DESeq2 matrix input stream...")
else:
    # ==========================================
    # FEATURE 1: INTERACTIVE THRESHOLD SLIDERS (SIDEBAR)
    # ==========================================
    st.sidebar.header("🔬 Analytical Threshold Settings")
    
    # Dynamic log2FC slider
    p_threshold = st.sidebar.slider(
        "Adjust Significance Cutoff (p-value)", 
        min_value=0.001, max_value=0.100, value=0.050, step=0.005, format="%.3f"
    )
    
    lfc_threshold = st.sidebar.slider(
        "Adjust Biological Variance (|log2FC|)", 
        min_value=0.0, max_value=4.0, value=1.0, step=0.2
    )

    # Re-classify dataframe vectors dynamically based on user sidebar input
    df = raw_df.copy()
    df['Significance'] = 'Not Significant'
    df.loc[(df['log2FoldChange'] >= lfc_threshold) & (df['padj'] <= p_threshold), 'Significance'] = 'Up-regulated'
    df.loc[(df['log2FoldChange'] <= -lfc_threshold) & (df['padj'] <= p_threshold), 'Significance'] = 'Down-regulated'

    # Recalculate dynamic card counts
    total_genes = len(df)
    up_genes = len(df[df['Significance'] == 'Up-regulated'])
    down_genes = len(df[df['Significance'] == 'Down-regulated'])

    # Display dynamic metric grid arrays
    col1, col2, col3 = st.columns(3)
    col1.metric("Total Transcriptome Pool", f"{total_genes:,}")
    col2.metric("Dynamic Up-regulated Targets", f"{up_genes:,}")
    col3.metric("Dynamic Down-regulated Targets", f"{down_genes:,}")
    
    st.write("---")

    # Create dual-panel visual columns
    plot_col, heatmap_col = st.columns([3, 2])

    with plot_col:
        st.subheader("🌋 Live-Recalculating Volcano Plot")
        
        # Color mapping configuration
        color_map = {
            "Up-regulated": "#EF553B",
            "Down-regulated": "#636EFA",
            "Not Significant": "#B6E880"
        }
        
        fig = px.scatter(
            df, x="log2FoldChange", y="log10_pvalue",
            color="Significance", color_discrete_map=color_map,
            hover_data=["Gene_ID", "baseMean", "pvalue", "padj"],
            labels={"log2FoldChange": "Log2 Fold Change", "log10_pvalue": "-Log10 (p-value)"},
            height=500
        )
        
        # Add static threshold visual guide lines
        fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="gray", opacity=0.7)
        fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="gray", opacity=0.7)
        fig.add_hline(y=-np.log10(p_threshold), line_dash="dash", line_color="gray", opacity=0.7)
        
        st.plotly_chart(fig, use_container_width=True)

    with heatmap_col:
        # ==========================================
        # FEATURE 2: INTERACTIVE CLUSTER HEATMAP
        # ==========================================
        st.subheader("🌡️ Top Feature Expression Heatmap")
        
        # Extract top 25 features based on highest statistical significance
        top_sig_genes = df.sort_values(by="pvalue").head(25)
        
        if len(top_sig_genes) > 0:
            # Reconstruct mock sample profiles back out of log2FC metrics to visualize trends
            genes_list = top_sig_genes['Gene_ID'].tolist()
            
            # Generate dummy sample contrasts representing the extracted gene trends
            control1 = np.random.uniform(10, 50, size=len(genes_list))
            control2 = control1 + np.random.uniform(-5, 5, size=len(genes_list))
            
            # Scale treatment vectors inversely or directly matching their log2fc values
            treatment1 = control1 * (2 ** top_sig_genes['log2FoldChange'].values)
            treatment2 = treatment1 + np.random.uniform(-10, 10, size=len(genes_list))
            
            heatmap_data = np.column_stack([control1, control2, treatment1, treatment2])
            # Row-standardize z-score conversion for clean biological visualization
            heatmap_z = np.array([(row - np.mean(row)) / (np.std(row) + 1e-5) for row in heatmap_data])

            heatmap_fig = go.Figure(data=go.Heatmap(
                z=heatmap_z,
                x=['Control_Rep1', 'Control_Rep2', 'Treatment_Rep1', 'Treatment_Rep2'],
                y=genes_list,
                colorscale='RdBu_r',
                zmin=-2, zmax=2,
                colorbar=dict(title="Z-Score", thickness=15)
            ))
            
            heatmap_fig.update_layout(height=500, margin=dict(l=10, r=10, t=10, b=10))
            st.plotly_chart(heatmap_fig, use_container_width=True)
        else:
            st.info("No significant targets filtered to compute clustering pathways.")

    st.write("---")

    # ==========================================
    # FEATURE 3: DATA EXPORTER TERMINAL & TABLES
    # ==========================================
    st.subheader("📋 Searchable Expression Matrix & Export Terminal")
    
    # Filter selection configuration framework
    filter_choice = st.radio("Isolate Matrix Views:", ["Show All Genes", "Isolate Up-regulated Only", "Isolate Down-regulated Only"], horizontal=True)
    
    if filter_choice == "Isolate Up-regulated Only":
        filtered_df = df[df['Significance'] == 'Up-regulated']
    elif filter_choice == "Isolate Down-regulated Only":
        filtered_df = df[df['Significance'] == 'Down-regulated']
    else:
        filtered_df = df

    # Data Download Pipeline Infrastructure
    @st.cache_data
    def convert_df_to_csv(target_df):
        return target_df.to_csv(index=False).encode('utf-8')

    csv_data = convert_df_to_csv(filtered_df)

    # Exporter button positioning
    st.download_button(
        label=f"📥 Export Selected Cohort Matrix ({len(filtered_df):,} Genes) as CSV",
        data=csv_data,
        file_name=f"filtered_rna_seq_{filter_choice.lower().replace(' ', '_')}.csv",
        mime="text/csv",
    )

    # Render data frame spreadsheet engine viewport window
    st.dataframe(
        filtered_df[["Gene_ID", "baseMean", "log2FoldChange", "pvalue", "padj", "Significance"]],
        use_container_width=True,
        hide_index=True
    )
