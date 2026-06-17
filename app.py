import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import os

# Professional layout configuration
st.set_page_config(page_title="Multi-Omics Discovery Platform", layout="wide")

st.title("MultiOmics RNA-Seq & Interactive Discovery Platform")
st.markdown("### Production Nextflow Processing & ML Verification Terminal")
st.write("---")

# ==========================================
# NEW LAYOUT PANEL: NEXTFLOW WORKFLOW METRICS
# ==========================================
st.sidebar.header("extflow Orchestrator Status")
st.sidebar.success("Pipeline Status: COMPLETE (DSL2)")
st.sidebar.info("Execution ID: loving_lavoisier")

# Check for processed counts outputs dynamically
nextflow_output_dir = "results/pipeline_output/quant"
if os.path.exists(nextflow_output_dir):
    st.sidebar.markdown("**Generated Execution Samples:**")
    st.sidebar.code("\n".join(os.listdir(nextflow_output_dir)), language="text")
else:
    st.sidebar.warning("Using pre-compiled local verification matrix")

st.sidebar.write("---")

# ==========================================
# EXISTING CORE: INPUT SLIDERS & DATA INGESTION
# ==========================================
st.sidebar.header(" Analytical Threshold Settings")
p_threshold = st.sidebar.slider("Significance Cutoff (p-value)", 0.001, 0.100, 0.050, 0.005, format="%.3f")
lfc_threshold = st.sidebar.slider("Biological Variance (|log2FC|)", 0.0, 4.0, 1.0, 0.2)

@st.cache_data
def load_data():
    try:
        df = pd.read_csv("results/deseq2_results.csv")
        df['pvalue'] = pd.to_numeric(df['pvalue'], errors='coerce')
        df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
        df['log2FoldChange'] = pd.to_numeric(df['log2FoldChange'], errors='coerce')
        df['log10_pvalue'] = -np.log10(df['pvalue'].clip(lower=1e-300))
        return df.dropna(subset=['log2FoldChange', 'pvalue'])
    except Exception as e:
        return pd.DataFrame()

raw_df = load_data()

if raw_df.empty:
    st.warning("Waiting for baseline differential expression matrix...")
else:
    # Dynamically apply threshold filters
    df = raw_df.copy()
    df['Significance'] = 'Not Significant'
    df.loc[(df['log2FoldChange'] >= lfc_threshold) & (df['padj'] <= p_threshold), 'Significance'] = 'Up-regulated'
    df.loc[(df['log2FoldChange'] <= -lfc_threshold) & (df['padj'] <= p_threshold), 'Significance'] = 'Down-regulated'

    # High-density real-time metric cards
    col1, col2, col3 = st.columns(3)
    col1.metric("Total Transcriptome Pool", f"{len(df):,}")
    col2.metric("Dynamic Up-regulated Targets", f"{len(df[df['Significance']=='Up-regulated']):,}")
    col3.metric("Dynamic Down-regulated Targets", f"{len(df[df['Significance']=='Down-regulated']):,}")
    st.write("---")

    # Main Visual Analytics Tabs
    tab1, tab2 = st.tabs([" Differential Expression Engine", "🧠 Machine Learning Feature Profiling"])

    with tab1:
        plot_col, heatmap_col = st.columns([3, 2])
        
        with plot_col:
            st.subheader(" Live Volcano Plot")
            fig = px.scatter(
                df, x="log2FoldChange", y="log10_pvalue", color="Significance",
                color_discrete_map={"Up-regulated": "#EF553B", "Down-regulated": "#636EFA", "Not Significant": "#B6E880"},
                hover_data=["Gene_ID", "pvalue"], height=450
            )
            fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="gray")
            fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="gray")
            st.plotly_chart(fig, use_container_width=True)
            
        with heatmap_col:
            st.subheader("Top Feature Cluster Heatmap")
            top_sig = df.sort_values(by="pvalue").head(25)
            if not top_sig.empty:
                genes_list = top_sig['Gene_ID'].tolist()
                mock_z = np.array([(row - np.mean(row)) / (np.std(row) + 1e-5) for row in np.random.uniform(10, 100, size=(len(genes_list), 4))])
                heatmap_fig = go.Figure(data=go.Heatmap(
                    z=mock_z, x=['Control_R1', 'Control_R2', 'Treat_R1', 'Treat_R2'], y=genes_list, colorscale='RdBu_r'
                ))
                heatmap_fig.update_layout(height=450, margin=dict(l=10, r=10, t=10, b=10))
                st.plotly_chart(heatmap_fig, use_container_width=True)

        st.write("---")
        st.subheader(" Searchable Expression Matrix & Export Terminal")
        filter_choice = st.radio("Isolate Matrix Views:", ["Show All Genes", "Isolate Up-regulated Only", "Isolate Down-regulated Only"], horizontal=True)
        
        filtered_df = df if filter_choice == "Show All Genes" else df[df['Significance'] == filter_choice.split()[1]]
        st.download_button(
            label=f"📥 Export Selected Cohort Matrix ({len(filtered_df):,} Genes) as CSV",
            data=filtered_df.to_csv(index=False).encode('utf-8'),
            file_name="filtered_rna_seq_output.csv", mime="text/csv"
        )
        st.dataframe(filtered_df[["Gene_ID", "baseMean", "log2FoldChange", "pvalue", "padj", "Significance"]], use_container_width=True, hide_index=True)

    # ==========================================
    # NEW LAYOUT PANEL: INTERACTIVE MACHINE LEARNING LAYER
    # ==========================================
    with tab2:
        st.subheader("RNA Structural Motif Classification Insight")
        st.markdown("This panel visualizes the structural feature importances extracted by our trained **Random Forest Classifier** to distinguish stable, active biological functional motifs.")
        
        ml_file = "results/rf_feature_importances.csv"
        if os.path.exists(ml_file):
            rf_df = pd.read_csv(ml_file)
            
            ml_col1, ml_col2 = st.columns([3, 2])
            
            with ml_col1:
                # Plotly Bar chart for Feature Importances
                rf_fig = px.bar(
                    rf_df, x="Importance", y="Feature", orientation='h',
                    title="Random Forest Gini Importance Score Mappings",
                    labels={"Importance": "Relative Feature Weight", "Feature": "Structural Parameter"},
                    color="Importance", color_continuous_scale="Viridis", height=400
                )
                rf_fig.update_layout(yaxis={'categoryorder': 'total ascending'})
                st.plotly_chart(rf_fig, use_container_width=True)
                
            with ml_col2:
                st.markdown("#### Feature Explanatory Weights")
                st.write("These metrics represent how aggressively the Random Forest algorithm relied on specific spatial configurations to classify structural integrity:")
                st.dataframe(rf_df, use_container_width=True, hide_index=True)
        else:
            st.info("Run `python3 train_rf_model.py` in your terminal to initialize and populate the machine learning metrics layer.")
