# snRNA_Project
codes used for my thesis project

sample database: GEO WEBSITE<br>
programming language used is: R<br>
sample datset: GSE174367 <br>
objective: Single-Nucleus Exploration of Cellular Mechanisms and Gene Expression Patterns in Alzheimer's Disease<br>
Methodology: <ol>
<li><a href="code_to_load.R">Data retrieval</li>
<li><a href="QC.R" >Quality Control </li>
<li><a href="batch effect.R" >Batch correction and Integration </li>
<li><a href="Norm_anno.R">normalization to cell annotation</a>
    <ul>
        <li>normalization</li>
        <li>finding variable gene</li>
        <li>scaling</li>
        <li>Dimension Reduction by: PCA, UMAP, t-SNE</li>
        <li>clustering</li>
    </ul>
</li>
<li><a href="DEG.R" >differential expression for disease gene identification</li>
</ol>       
