import gzip
import random

def generate_fastq(filename, pair_suffix):
    bases = ['A', 'C', 'G', 'T']
    # Simulated gene sequences to give the aligner something real to map
    target_seq = "ATGGCGTACGCGTAATCAAAGCAGCGATTGCAAAAGCTCAGTTTTG"
    
    with gzip.open(filename, 'wt') as f:
        for i in range(1000):
            f.write(f"@ERR123456.{i} SEQ_{i}/{pair_suffix}\n")
            # Mix some matching target sequences with random sequencing noise
            if random.random() > 0.3:
                seq = target_seq + "".join(random.choices(bases, k=54))
            else:
                seq = "".join(random.choices(bases, k=100))
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write("I" * 100 + "\n") # High quality Phred scores

generate_fastq("data/sample_1.fastq.gz", "1")
generate_fastq("data/sample_2.fastq.gz", "2")
print("Successfully generated valid paired-end FASTQ datasets locally!")


