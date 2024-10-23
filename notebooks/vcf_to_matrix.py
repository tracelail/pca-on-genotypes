from pysam import VariantFile
import numpy as np
from sklearn import decomposition
import pandas as pd

vcf_filename = "data/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf"

genotypes = []
samples = []
variant_ids = []

with VariantFile(vcf_filename) as vcf_reader:
    counter = 0
    for record in vcf_reader:
        counter += 1
        if counter % 100 == 0:
            alleles = [record.samples[x].allele_indices for x in record.samples]
            samples = [sample for sample in record.samples]
            genotypes.append(alleles)
            variant_ids.append(record.id)
        if counter % 4943 == 0:
            print(f"{round(100*counter/494328)}%")

genotypes = np.array(genotypes)
matrix = np.count_nonzero(genotypes, axis=2)
matrix = matrix.T
print(matrix.shape)

df = pd.DataFrame(matrix, columns=variant_ids, index=samples)
df.to_csv("matrix.csv")
