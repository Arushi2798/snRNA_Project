# snRNA_Project
codes used for my thesis project

sample database: GEO WEBSITE<br>
programming language used is: R<br>
sample datset: GSE174367 <br>
objective: Single-Nucleus Exploration of Cellular Mechanisms and Gene Expression Patterns in Alzheimer's Disease<br>
Methodology: <ol>
<li> Data retrival <a href="/snRNA_Project/code_to_load.R"></li>
<li>quality control <a href="/snRNA_Project/QC.R" ></li>
<li>normalization to cell annotation<a href="/snRNA_Project/Norm_anno.R">
    <ul>
        <li>normalization</li>
        <li>finding variable gene</li>
        <li>scaling</li>
        <li>Dimension Reduction by: PCA, UMAP, t-SNE</li>
    </ul>
</li>
<li>feature extraction and dimensionality reduction</li>
<li>clustering</li>
<li>cell type inference and annotation</li>
<li>differential expression for disease gene identification</li>
</ol>
