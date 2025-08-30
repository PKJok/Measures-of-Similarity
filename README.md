Measures of Similarity
similarity_index Function
The similarity_index function calculates a similarity score between two genotypes based on their genetic data, implementing a scoring rule where:

Score 2: Both alleles at a locus are identical.
Score 1: Only one allele is shared.
Score 0: No alleles are shared.

The score is averaged across all loci to produce a normalized similarity value between 0 and 1.
Function Definition
similarity_index <- function(x, y) {
  sum(x * y) / (2 * length(x))
}

This function takes two vectors, x and y, representing the genotypic scores (0, 1, or 2) for two genotypes across multiple loci, and computes their similarity using the formula sum(x * y) / (2 * length(x)).
Pairwise Similarity Matrix
The function is typically used within a nested loop to compute a similarity matrix for all pairs of genotypes in a dataset:
similarity_matrix <- matrix(0, nrow = n_genotypes, ncol = n_genotypes)
for (i in 1:n_genotypes) {
  for (j in 1:n_genotypes) {
    similarity_matrix[i, j] <- similarity_index(genotype_data[i, ], genotype_data[j, ])
  }
}

Here, genotype_data is a matrix where each row represents a genotype and each column represents a genetic marker (e.g., SNP).
Processing of x and y
When similarity_matrix[i, j] <- similarity_index(genotype_data[i, ], genotype_data[j, ]) is executed, the following steps occur:

Data Extraction:

genotype_data[i, ]: Extracts all columns for row i, returning a vector of genotypic scores (e.g., (0, 1, 2, ..., 1)) for genotype i. This vector is assigned to x.
genotype_data[j, ]: Similarly extracts the row for genotype j, assigned to y.
Both x and y are vectors of equal length, corresponding to the number of markers (e.g., 100 SNPs).


Element-wise Multiplication (x * y):

R performs element-wise multiplication between x and y. For example, if x = (0, 1, 2) and y = (0, 2, 2), the result is:
0 * 0 = 0
1 * 2 = 2
2 * 2 = 4
Resulting vector: (0, 2, 4).




Scoring Rule Encoding:

The input encoding (0, 1, 2) represents allelic states:
0: Homozygous reference (e.g., AA).
1: Heterozygous (e.g., Aa).
2: Homozygous alternate (e.g., aa).


The multiplication x * y produces values (0, 1, 2, 4) that act as a proxy for similarity:
4: High similarity (e.g., both homozygous alternate).
1: High similarity (e.g., both heterozygous).
2: Intermediate similarity (e.g., heterozygous vs. homozygous alternate).
0: No similarity (e.g., no shared alleles).


This does not directly match the textbook scores (2, 1, 0) but reflects the same relative similarity due to the encoding.


Summation and Normalization:

sum(x * y) sums the products across all loci (e.g., 0 + 2 + 4 = 6).
2 * length(x) computes the maximum possible score for identical genotypes (e.g., 2 * 3 = 6 for 3 markers).
The normalized similarity is sum(x * y) / (2 * length(x)) (e.g., 6 / 6 = 1), yielding a value between 0 and 1.



Example
Consider a dataset with 3 genotypes and 3 markers:



Genotype
SNP1
SNP2
SNP3



1
0
1
2


2
0
0
1


3
2
1
0


For genotypes 1 and 2:

x = (0, 1, 2), y = (0, 0, 1).
x * y = (0 * 0, 1 * 0, 2 * 1) = (0, 0, 2).
sum(x * y) = 0 + 0 + 2 = 2.
2 * length(x) = 2 * 3 = 6.
Similarity = 2 / 6 = 0.333.

The resulting similarity_matrix is a symmetric 3x3 matrix containing similarity scores for all genotype pairs.
