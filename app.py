import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import os

st.set_page_config(page_title="Multi Omics NGS Dashboard", layout="wide")

st.title("Multi Omics RNA Seq Differential Expression Dashboard")
st.markdown("Automated Nextflow Processing & Verification Terminal")

results_path = "results/deseq2_results.csv"

if os.path.exists(results_path):
    # Read matrix and fix unnamed column to reflect Gene ID
    df = pd.read_csv(results_path)
    df.rename(columns={df.columns[0]: "Gene_ID"}, inplace=True)
    
    # Calculate -log10(p-value) for standard volcano plotting
    df['-log10_pvalue'] = -np.log10(df['pvalue'].replace(0, np.nan))
    df['-log10_pvalue'] = df['-log10_pvalue'].fillna(df['-log10_pvalue'].max())

    # Metric summary rows
    sig_threshold = 0.05
    lfc_threshold = 1.0
    
    up_regulated = df[(df['pvalue'] < sig_threshold) & (df['log2FoldChange'] > lfc_threshold)].shape[0]
    down_regulated = df[(df['pvalue'] < sig_threshold) & (df['log2FoldChange'] < -lfc_threshold)].shape[0]

    col1, col2, col3 = st.columns(3)
    col1.metric("Total Genes Analyzed", f"{df.shape[0]}")
    col2.metric("Up regulated Genes", f"{up_regulated}")
    col3.metric("Down regulated Genes", f"{down_regulated}")

    st.markdown("---")
    
    # Create the Interactive Volcano Plot
    st.subheader(" Volcano Plot Visualization")
    
    df['Significance'] = 'Not Significant'
    df.loc[(df['pvalue'] < sig_threshold) & (df['log2FoldChange'] > lfc_threshold), 'Significance'] = 'Up regulated'
    df.loc[(df['pvalue'] < sig_threshold) & (df['log2FoldChange'] < -lfc_threshold), 'Significance'] = 'Down regulated'
    
    fig = px.scatter(
        df, x="log2FoldChange", y="-log10_pvalue",
        color="Significance",
        hover_data=["Gene_ID", "pvalue", "baseMean"],
        color_discrete_map={'Not Significant': 'grey', 'Up regulated': 'red', 'Down regulated': 'blue'},
        labels={"log2FoldChange": "Log2 Fold Change", "-log10_pvalue": "-log10(p-value)"},
        template="plotly_dark"
    )
    st.plotly_chart(fig, use_container_width=True)

    # Data Table Viewer
    st.subheader(" Interactive Expression Matrix Data")
    st.dataframe(df, use_container_width=True)
else:
    st.warning(f"Waiting for pipeline results... Could not find a matrix file at `{results_path}`.")
