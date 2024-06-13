
cut -f 1-11 ../../2.impute/Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.LIMS1 > LIMS1.SNPs

ml tabix/0.2.6

ml samtools/1.9

bgzip ForLocusZoom.tab

tabix -s 1 -b 2 -e 2 -S 1 -f ForLocusZoom.tab.gz

zcat ForLocusZoom.tab.gz | awk 'BEGIN{printf("MarkerName\tP-value\n")} (NR>1){printf("chr%s:%s\t%s\n", $1, $2, $5)}' > ForLocusZoom.app.tab

ml locuszoom/1.4

#locuszoom --metal ForLocusZoom.app.tab --source 1000G_Nov2014 --build hg19 --pop AMR+ASN+AFR+EUR --refgene LIMS1 

locuszoom --metal ForLocusZoom.app.tab --source 1000G_Nov2014 --build hg19 --pop EUR --refgene LIMS1 

locuszoom --metal ForLocusZoom.app.tab --source 1000G_Nov2014 --build hg19 --pop EUR --refsnp rs200106875 --chr 2 --start 109065017 --end 109306702 width=15 showRefsnpAnnot=FALSE ymax=5 


## SNP30

zcat /sc/arion/projects/GOCAR/zhangz05/raw_data/imputation_2nd_round/MSSM/MSSM_GoNL_1KG_chr2:105-110Mb_info.gz | head -1 > SNP30.info

awk '(NR==FNR){tt[$3] = $0; next} {if(tt[$3] != ""){printf("%s\t%s\t%s\t%s\t%s\n", tt[$3], $2, $7, $5, $6)}}' 30SNPs <(zcat /sc/arion/projects/GOCAR/zhangz05/raw_data/imputation_2nd_round/MSSM/MSSM_GoNL_1KG_chr2:105-110Mb_info.gz) >> SNP30.info

#awk '(NR==FNR){tt[$3] = 1; next} {if(tt[$3] == ""){print $0}}' SNP30.info 30SNPs

#awk '(NR==FNR){tt[$3] = $2"\t"$7"\t"$5"\t"$6; next} {if(tt[$3]!=""){printf("%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $5, $6, tt[$3])}}' SNP30.info SNP30.vcf.bk > tt && mv tt SNP30.info

#awk '(NR==FNR){tt[$3] = $2"\t"$7"\t"$5"\t"$6; next} {if(tt[$3]==""){printf("%s\t%s\t%s\t%s\t%s\t\n", $1, $2, $3, $5, $6)}}' SNP30.info SNP30.vcf > tt && cat tt >> SNP30.info

awk '(NR==FNR){tt[$3] = $0; next} {split($1, tt2, ":"); if(tt[tt2[2]]!=""){printf("%s\t%s\t%s\t%s\n", tt[tt2[2]], $1, $8, $7)}}' SNP30.info <(zcat /sc/arion/projects/GOCAR/zhangz05/CTOT/data/genotype_imputation/imputation_KGP_phase3/Imputation_results/chr2.info.gz) > tt && mv tt SNP30.info


awk '(NR==FNR){tt[$1] = 1; next} {if(tt[$1]== ""){print $1}}' SNP30.info /sc/arion/projects/GOCAR/Sun/5.DRWork/2.impute/Diff_homo/SNP30/SNP30.info

#awk '(NR==FNR){tt[$1]=1; next} {if(tt[$2]!=""){print $0}}' 30SNPs <(zcat ../../2.impute/all.tped.gz) > SNP30.tped

awk '(NR==FNR){tt[$3]=1; next} {if(tt[$3]!=""){print $0}}' 30SNPs.bk <(zcat ../../2.impute/Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.gz) | awk '$2 == 2'> SNP30.vcf

cut -f 1,59- SNP30.vcf > SNP30.tab


zcat /sc/arion/projects/GOCAR/zhangz05/CTOT/data/genotype_imputation/imputation_KGP_phase3/Imputation_results/chr2.dose.vcf.gz | head -19 > SNP30.CTOT.vcf

awk '(NR==FNR){tt[$3] = 1; next} {if(tt[$2]!=""){print $0}}' 30SNPs <(zcat /sc/arion/projects/GOCAR/zhangz05/CTOT/data/genotype_imputation/imputation_KGP_phase3/Imputation_results/chr2.dose.vcf.gz) >> SNP30.CTOT.vcf

gzip SNP30.CTOT.vcf 

module load vcftools/0.1.15 && vcftools --gzvcf SNP30.CTOT.vcf.gz --extract-FORMAT-info GT --out SNP30.CTOT.vcf

        head -1 SNP30.CTOT.vcf.GT.FORMAT | awk '{printf("ID\t%s\n", $0)}' > SNP30.CTOT.vcf.GT.FORMAT.withID 

        awk 'BEGIN{FS="\t"; OFS = "\t"}(NR>1){printf("%s:%s\t%s\n", $1, $2, $0)}' SNP30.CTOT.vcf.GT.FORMAT >> SNP30.CTOT.vcf.GT.FORMAT.withID 

awk '(NR==FNR){tt[$3]=1;next}{if(tt[$3]==""){print $3}}' <(bcftools view -H SNP30.CTOT.vcf.gz) <(bcftools view -H /sc/arion/projects/GOCAR/Sun/5.DRWork/4.CTOT2/5.SNP30.2/chr2.dose.vcf.gz)


awk '{printf("chr%s\t%s\t%s\n", $2, $3, $3)}' 30SNPs > 30SNPs.bed


bcftools view -h /sc/arion/projects/GOCAR/zhangz05/1000G_30X/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz > SNP30.1000G.vcf

awk '(NR==FNR){tt[$2] = 1; next} {if(tt[$2]!=""){print $0}}' 30SNPs.hg38.bed <(zcat /sc/arion/projects/GOCAR/zhangz05/1000G_30X/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz) >> SNP30.1000G.vcf

gzip SNP30.1000G.vcf

module load vcftools/0.1.15 && vcftools --gzvcf SNP30.1000G.vcf.gz --extract-FORMAT-info GT --out SNP30.1000G.vcf

        head -1 SNP30.1000G.vcf.GT.FORMAT | awk '{printf("ID\t%s\n", $0)}' > SNP30.1000G.vcf.GT.FORMAT.withID

        awk 'BEGIN{FS="\t"; OFS = "\t"}(NR>1){printf("%s:%s\t%s\n", $1, $2, $0)}' SNP30.1000G.vcf.GT.FORMAT >> SNP30.1000G.vcf.GT.FORMAT.withID


awk '{printf("chr%s\t%s\t%s\n", $2, $3, $3)}' 30SNPs > 30SNPs.bed



#bcftools view -R 30SNPs.hg38.bed.bk /sc/arion/projects/GOCAR/Sun/5.DRWork/PaperLIMS1/FigNEJ/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -O z -o SNP30.1000G.vcf.gz

#bcftools view -R 30SNPs.bed /sc/arion/projects/GOCAR/Sun/5.DRWork/PaperLIMS1/FigNEJ/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -O z -o SNP30.1000G.vcf.gz

awk '($0!~/^#/){print $23}' 1000G_2504_high_coverage.sequence.index > Sam2504


## UCSC customized track

awk '{printf("chr%s\t%s\t%s\t%s\n", $2, $3, $4+1, $1)}' SNP30.vcf > CusTracks

