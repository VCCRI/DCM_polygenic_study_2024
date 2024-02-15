# DCM_polygenic_study_2024
A 2024 study on polygenic risk scores and familial inheritance of dilated cardiomyopathy

This repository contains the code for a 2024 study on the role of polygenic risk scores in the inheritance of dilated cardiomyopathy. Briefly, the methodology for this study identified odds ratios for disease outcomes using logistic mixed effects models, accounting for sample relatedness.

This study is currently under review as: Thompson JM, Johnson R, Troup M, Rath EM, Young PE, Soka MJ, Ohanian M, Tarr IS, Giannoulatou E, Fatkin D. Polygenic risk in families with dilated cardiomyopathy. (2024). Circulation Genomics and Precision Medicine.

Genomic data for this study were accessed through Australian Genomics, funded by the NHMRC (1113531, 2000001) and the MRFF, administered by the Murdoch Children's Research Institute. The Australian Genomics Cardiovascular Disorders Flagship was funded by the MRFF Genomic Health Futures Mission (EPCD000028). 

The DCM_families_prs_analysis.R script takes scoring input produced by plink with the following command:

plink --bcf DCM_patients.bcf --geno 0 --double-id --allow-no-sex --vcf-filter --score score_file.txt header sum --out output_dir --allow-extra-chr --write-snplist 


The output has the following columns from plink in the header:
"NCI_ID"	"FID"	"IID"	"PHENO"	"CNT"	"CNT2"	"SCORESUM"

As well as added columns with family/phenotypic information:
"affected_relative"  "unaffected_relative"  "control"  "proband"  "monogenic"

This file is the basis of the "dcm_families_prs.csv" imported by the DCM_families_prs_analysis.R script.


The required "rel_mat.rel.id_edited" and "rel_mat.rel" files are also produced by plink with the following command:

plink --double-id --vcf DCM_patients.vcf --vcf-half-call m --maf 0.05 --geno 0.03 --make-rel square --out output_dir/rel_mat
