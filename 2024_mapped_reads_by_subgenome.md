# Quantifying mapped reads in each subgenome

An interesting question raised by a reviewer asks whether there is a siginficiantly different proportion of mapped reads in each subgenome. Because the subgenomes are different sizes I can create an expectation based on the size differences and use a chi square test to evaluate whether the actual number in each subgenome is different from the expectation.

First subset the bam file by subgenome:
```
samtools view -bo muel_Z23765_M_CAGA_sorted_Lsubgenome.bam muel_Z23765_M_CAGA_sorted.bam --region-file XL_L_subgenome.bed
samtools view -bo muel_Z23765_M_CAGA_sorted_sorted_Ssubgenome.bam muel_Z23765_M_CAGA_sorted.bam --region-file XL_S_subgenome.bed
```
now count the mapped reads in each subgenome:
```
samtools flagstat muel_Z23765_M_CAGA_sorted_Lsubgenome.bam
samtools flagstat muel_Z23765_M_CAGA_sorted_sorted_Ssubgenome.bam
```
