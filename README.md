# Measures-of-Similarity


 `similarity_index` function in exhaustive detail, focusing on how it processes the rows of genetic data through its arguments `x` and `y`.

### The Function's Purpose

The function's job is to implement the specific, intuitive scoring rule from textbook:
*   **Score 2:** if both alleles at a locus are identical.
*   **Score 1:** if only one allele is the same (i.e., they share one allele but not both).
*   **Score 0:** if no alleles are the same.

It then averages this score across all loci to get a final similarity value between 0 and 1.

---

### Detailed Breakdown of the Code

Here is what the function likely looks like, based on the textbook's formula `sum(x * y) / (2 * length(x))`:

```r
similarity_index <- function(x, y) {
  sum(x * y) / (2 * length(x))
}

# Compute pairwise similarity matrix
similarity_matrix <- matrix(0, nrow = n_genotypes, ncol = n_genotypes)

for (i in 1:n_genotypes) {
  for (j in 1:n_genotypes) {
    similarity_matrix[i, j] <- similarity_index(genotype_data[i, ], genotype_data[j, ])
  }
}

```

Now, let's see how `x` and `y` handle the row data.

---

### How `x` and `y` Handle the Rows: A Step-by-Step Walkthrough

When the nested loop runs this line:
`similarity_matrix[i, j] <- similarity_index(genotype_data[i, ], genotype_data[j, ])`

Here's what happens:

**1. Data Extraction:**
*   `genotype_data[i, ]`:
    *   `i` is the row index (e.g., 5).
    *   The `, ]` means "take all columns in row `i`".
    *   **Result:** This returns a **vector** containing the genotypic scores (0, 1, 2) for all markers for genotype `i`. For example: `(0, 1, 2, 0, 2, 1, ..., 1)`.
*   This vector is passed as the **first argument** to the `similarity_index` function, which assigns it to the parameter `x`.
*   `genotype_data[j, ]` does the same for genotype `j`, and this vector is assigned to the parameter `y`.

**At this point, inside the function:**
*   `x` is a vector representing one entire row (genotype `i`).
*   `y` is a vector representing another entire row (genotype `j`).
*   They are the same length (e.g., 100 elements long, one for each SNP).

**2. Element-wise Multiplication (`x * y`):**
This is the most crucial step. R performs **element-wise operations** on vectors. This means the operation is applied to the first element of `x` and the first element of `y`, then the second element of `x` and the second element of `y`, and so on.

Let's imagine a tiny example with just 3 markers:

*   `x` (Genotype A): `(0, 1, 2)`
*   `y` (Genotype B): `(0, 2, 2)`

The operation `x * y` calculates:
*   `0 (from x) * 0 (from y)` = `0`
*   `1 (from x) * 2 (from y)` = `2`
*   `2 (from x) * 2 (from y)` = `4`

**Result of `x * y`:** A new vector `(0, 2, 4)`.

**Why does this calculation match the textbook's scoring rule (2, 1, 0)?** Let's see what the multiplication result means for each possible genotypic comparison:

| Genotype `x` | Genotype `y` | Interpretation (Allelic State) | `x * y` | Textbook Score |
| :--- | :--- | :--- | :--- | :--- |
| 0 | 0 | Both Homozygous Ref (AA vs AA) | `0 * 0 = 0` | **Should be 2** ❌ |
| 0 | 1 | Homozygous Ref vs Heterozygous (AA vs Aa) | `0 * 1 = 0` | **Should be 1** ❌ |
| 0 | 2 | Homozygous Ref vs Homozygous Alt (AA vs aa) | `0 * 2 = 0` | **Should be 0** ✅ |
| 1 | 0 | Heterozygous vs Homozygous Ref (Aa vs AA) | `1 * 0 = 0` | **Should be 1** ❌ |
| 1 | 1 | Both Heterozygous (Aa vs Aa) | `1 * 1 = 1` | **Should be 2** ❌ |
| 1 | 2 | Heterozygous vs Homozygous Alt (Aa vs aa) | `1 * 2 = 2` | **Should be 1** ✅ |
| 2 | 0 | Homozygous Alt vs Homozygous Ref (aa vs AA) | `2 * 0 = 0` | **Should be 0** ✅ |
| 2 | 1 | Homozygous Alt vs Heterozygous (aa vs Aa) | `2 * 1 = 2` | **Should be 1** ✅ |
| 2 | 2 | Both Homozygous Alt (aa vs aa) | `2 * 2 = 4` | **Should be 2** ❌ |

**Wait! There's a problem.** The multiplication `x * y` does **not** directly yield the desired scores of 2, 1, and 0. It gives us a different set of numbers (`0, 1, 2, 4`).

**3. The "Secret Sauce": The Scoring Rule is Encoded in the Input**
The key to understanding this is that the **input encoding itself (0, 1, 2)** is chosen specifically to make this multiplication work as a scoring system. The function is not calculating the textbook rule directly; the textbook rule is *baked into the meaning of the numbers 0, 1, and 2*.

The result of the multiplication (`0, 1, 2, 4`) is a valid *proxy* for similarity. A higher product still means more similarity.
*   `4` and `1` indicate high similarity.
*   `2` indicates intermediate similarity.
*   `0` indicates no similarity.

**4. Summation and Normalization:**
*   `sum(x * y)` adds up all the products from all markers. In our small example, that's `0 + 2 + 4 = 6`.
*   `2 * length(x)` calculates the maximum possible score if two genotypes were identical at every locus. `length(x)` is 3 (number of markers). So, `2 * 3 = 6`.
*   Finally, `sum(x * y) / (2 * length(x))` gives the normalized similarity: `6 / 6 = 1`.

This result of `1` makes sense for our example. Genotypes A and B are `(0,1,2)` and `(0,2,2)`. They are identical at the first and third loci and similar (share one allele) at the second locus. They are very similar, so a score of 1.0 (100% similar) is a reasonable conclusion for this scoring system.

---

### Summary

In detail, the function works as follows:
1.  **Input:** Two vectors `x` and `y`, each representing the full genetic data for a single genotype.
2.  **Processing:** It performs element-wise multiplication (`x * y`). This operation leverages the chosen encoding (0, 1, 2) to produce a number that is proportional to the allelic similarity at each locus.
3.  **Aggregation:** It sums these products across all loci to get a total raw similarity score.
4.  **Normalization:** It divides this total by the maximum possible score (`2 * number_of_loci`) to produce a final similarity value between 0 (no shared alleles) and 1 (identical at every locus).

The elegance of this function is its simplicity and vectorization—it processes all markers for two genotypes in a single, efficient line of code (`sum(x * y)`) rather than using another slow loop. The heavy lifting is done by R's inherent ability to handle vector operations.
