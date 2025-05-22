# Gene Analysis Using UpSet Plots: Techniques, Metastasis Frequency, and Sample Size

## Aim & Objective
The aim of this project is to visualize and analyze the relationships between different gene sequencing techniques (NGS, WGS, WES) and their association with metastasis frequency across various cancer-related genes. The objectives are:
- To demonstrate the use of UpSet plots for visualizing complex set relationships
- To analyze the distribution of sequencing techniques among important cancer genes
- To examine the correlation between technique combinations and metastasis frequency
- To visualize sample sizes across different technique combinations

## Introduction
Cancer genomics research utilizes various sequencing techniques to identify genetic mutations. Next-Generation Sequencing (NGS), Whole Genome Sequencing (WGS), and Whole Exome Sequencing (WES) are commonly used either alone or in combination. This project uses UpSet plots - an innovative visualization technique for quantitative analysis of set intersections - to examine:
- Which technique combinations are most frequently used
- How technique usage relates to metastasis frequency
- The sample sizes associated with different technique combinations

## Materials
- Python 3.x
- Required Python packages:
  - pandas
  - matplotlib
  - upsetplot
- Sample dataset containing:
  - 20 cancer-related genes
  - Metastasis frequency data
  - Sample population sizes
  - Sequencing technique information

## Methodology
1. **Data Preparation**:
   - Created a sample dataset with gene names, metastasis frequencies, sample sizes, and technique combinations
   - Split the technique combinations into boolean columns using one-hot encoding

2. **UpSet Plot Construction**:
   - Converted the technique indicators into UpSet-compatible format
   - Configured the UpSet plot to sort by set cardinality
   - Added strip plots for metastasis frequency and sample size visualization

3. **Visualization**:
   - Generated the composite UpSet plot with intersection bars and attribute dots
   - Added custom legends and titles for clarity
   - Included supplementary technique usage statistics

## Results
1. **Visual Output**:
   - The UpSet plot shows technique combinations as intersecting bars
   - Blue dots represent metastasis frequency across technique combinations
   - Red dots represent total sample sizes for each combination

2. **Key Findings**:
   - The most common technique combinations were [list most common combinations from your data]
   - Technique combination [X] showed the highest average metastasis frequency
   - The largest sample sizes were associated with [specific technique combinations]

3. **Supplementary Statistics**:
   - Technique counts: [NGS: X, WGS: Y, WES: Z]
   - Gene-technique mappings for all 20 genes

## Conclusion
The UpSet plot visualization effectively revealed patterns in sequencing technique usage across cancer genes and their relationship to metastasis frequency. Key conclusions include:
- Certain technique combinations are preferentially used for specific genes
- Metastasis frequency varies across different technique combinations
- Sample sizes are not uniformly distributed across all technique combinations
- The visualization approach provides actionable insights for planning future genomic studies

This analysis demonstrates the utility of UpSet plots for multidimensional genomic data visualization and comparative analysis of experimental approaches in cancer research.
