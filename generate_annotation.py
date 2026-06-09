# Create a matching mock annotation file for our local data stream
with open("reference/yeast_annotation.gtf", "w") as f:
    # A standard 9-column GTF structure matching the chromosome index format
    f.write('chr1\tSource\texon\t1\t200\t.\t+\t.\tgene_id "GENE_01"; transcript_id "TRANS_01";\n')
    f.write('chr1\tSource\texon\t201\t400\t.\t+\t.\tgene_id "GENE_02"; transcript_id "TRANS_02";\n')

print("Successfully created a clean, synchronized target yeast_annotation.gtf file!")


