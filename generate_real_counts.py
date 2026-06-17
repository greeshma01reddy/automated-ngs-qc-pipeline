import os
import pandas as pd
import numpy as np

print("Generating high-density production count matrix (6,000+ Yeast Genes)...")

# 1. Create a realistic list of ~6,000 Systematic Yeast Gene IDs
yeast_chromosomes = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
gene_ids = []
for chrom in yeast_chromosomes:
    for g in range(1, 380): 
        strand = 'C' if g % 2 == 0 else 'W'
        gene_ids.append(f"Y{chrom}{g:03d}{strand}")

# 2. Simulate raw sequencing featureCounts metadata structure
df = pd.DataFrame({
    'Geneid': gene_ids,
    'Chr': [f"chr{np.random.choice(yeast_chromosomes)}" for _ in gene_ids],
    'Start': np.random.randint(100, 1000000, size=len(gene_ids)),
    'End': np.random.randint(1001, 1002000, size=len(gene_ids)),
    'Strand': [np.random.choice(['+', '-']) for _ in gene_ids],
    'Length': np.random.randint(500, 5000, size=len(gene_ids))
})

# 3. Simulate deep RNA-seq count columns with biological variance
df['Control_rep1'] = np.random.negative_binomial(n=10, p=0.01, size=len(gene_ids))
df['Control_rep2'] = df['Control_rep1'] + np.random.randint(-15, 15, size=len(gene_ids))

df['Treatment_rep1'] = df['Control_rep1'].copy()
df['Treatment_rep2'] = df['Control_rep2'].copy()

# Inject significant biological differential expression spikes
up_regulated_indices = np.random.choice(len(gene_ids), size=500, replace=False)
down_regulated_indices = np.random.choice(list(set(range(len(gene_ids))) - set(up_regulated_indices)), size=500, replace=False)

df.loc[up_regulated_indices, ['Treatment_rep1', 'Treatment_rep2']] *= np.random.randint(4, 12, size=500)[:, None]
df.loc[down_regulated_indices, ['Treatment_rep1', 'Treatment_rep2']] = (df.loc[down_regulated_indices, ['Treatment_rep1', 'Treatment_rep2']] / np.random.randint(4, 8, size=500)[:, None]).astype(int)

# Ensure no negative count values slip in
for col in ['Control_rep1', 'Control_rep2', 'Treatment_rep1', 'Treatment_rep2']:
    df[col] = df[col].clip(lower=0)

# 4. Save to match exactly what featureCounts outputs (Safe isolated file blocks)
output_path = "sample_counts.txt"

# First write the header line and close the file
with open(output_path, 'w') as f:
    f.write("# featureCounts v2.0.3; Command line invoked simulated transcriptomics matrix\n")

# Append the dataframe cleanly via its path string
df.to_csv(output_path, sep='\t', index=False, mode='a')

print(f"Success! Master counts matrix exported containing {len(gene_ids)} unique genes.")

