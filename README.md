#  MultiOmics Automated RNASeq & Interactive Discovery Platform

An end to end parallelized transcriptomics workflow framework that integrates high throughput sequencing secondary analysis, adaptive statistical differential expression testing, and downstream machine learning classification.

---

##  System Architecture Overview

The platform is designed in two primary analytical tracks to handle data execution and predictive insights independently:

1. Track 1 (Backend Pipeline Orchestration): A scalable Nextflow (DSL2) pipeline that parallelizes high throughput quality control filtering (`FastQC`) and transcript expression map estimation.
2. Track 2 (Downstream Statistical & ML Engine): An adaptive DESeq2 core in R that leverages gene wise dispersion estimations to robustly handle high variance inputs. Downstream data is fed into a **Random Forest Classifier** (`scikit-learn`) to predict and map RNA structural motif stability based on spatial configurations.
3. Interactive Viewport: A full stack Streamlit dashboard delivering interactive visualization engines powered by Plotly.

---

##  Features Built & Deployed

###  Core Analytics & Statistical Resilience
* Automated Asynchronous Parallelization: Leverages Nextflow channels to stream fastq inputs through parallel processing nodes.
* Resilient Dispersion Modeling:The R backend dynamically overrides standard parametric curve fitting limits on high variance inputs using gene wise shrinkages to effortlessly process over 6,000 yeast genes.
* Dynamic Volcano Visualizer: High density Plotly visualization coupled with real time sidebar sliders allowing users to dynamically manipulate $\log_2\text{Fold Change}$ and p value cutoffs.
* Row Standardized Z-Score Heatmaps: Extracts the top 25 most statistically significant genetic features to render dynamic expression contrasts between control and treatment replicates.

###  Predictive Machine Learning Integration
* Structural Motif Classification: Employs a tuned Random Forest model to calculate Gini feature importance across key spatial variables (e.g., Free Energy $\Delta G$, GC Content, Stem/Loop length configurations).
* On The Fly Cohort Exporter:A custom data streaming terminal allows researchers to filter expression matrices dynamically and export localized target gene cohorts directly as CSV files.

---


