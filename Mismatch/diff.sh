
set -e

head -1 ./Split/xaaVariant_Matrix.xls > Variant_Matrix.xls

for xx in ./Split/x*.xls
do

	sed '1d' $xx >> Variant_Matrix.xls

done

cat Split2/*.naRate.xls > NARate.xls

awk '$2<=0.05' NARate.xls > NARate.xls.filtered

perl monomorphic.pl ../all.tped.gz > Monomorphic 

awk '(NR==FNR){tt[$1] = 1; next} {if(tt[$1]==""){print $0}}' Monomorphic NARate.xls.filtered > NARate.xls.filtered.all

zcat Variant_Matrix.xls.gz | head -1 | cut -f2- > All.freq.xls

for ff in Split2/*freq.xls
do

        cat $ff >> All.freq.xls

done


zcat Variant_Matrix.xls.gz | head -1 > Variant_Matrix.xls.All

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.gz) <(zcat Variant_Matrix.xls.gz) >> Variant_Matrix.xls.All

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]==""){print $0}}' <(zcat Variant_Matrix.xls.gz) <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.gz) > MissingDiff


zcat Variant_Matrix.xls.All.gz | head -1 > Variant_Matrix.xls.All.filt

awk '(NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' NARate.xls.filtered.all <(zcat Variant_Matrix.xls.All.gz) >> Variant_Matrix.xls.All.filt

gzip Variant_Matrix.xls.All.filt


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.All.filt.comm

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR == FNR){if($2 >= 0.05){tt[$1] = 1}; next} {if(tt[$1] != ""){print $0;}}' All.mutfreq.xls <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.All.filt.comm

gzip Variant_Matrix.xls.All.filt.comm


#zcat Variant_Matrix.xls.All.filt.gz | sort | uniq > Variant_Matrix.xls.All.filt.uniq
#
#gzip Variant_Matrix.xls.All.filt.uniq

## GenomeWide Mismatch

zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.All.filt.Exon

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.Exonic.gz) <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.All.filt.Exon

gzip Variant_Matrix.xls.All.filt.Exon

#Rscript diff_sum_homo.r Variant_Matrix.xls.All.filt.Exon.gz 

Rscript diff_sum_heta.r Variant_Matrix.xls.All.filt.Exon.gz


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.All.filt.syno

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.syno.gz) <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.All.filt.syno

gzip Variant_Matrix.xls.All.filt.syno

#Rscript diff_sum_homo.r Variant_Matrix.xls.All.filt.syno.gz 

Rscript diff_sum_heta.r Variant_Matrix.xls.All.filt.syno.gz


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.All.filt.syno.NoSyAll.trmem

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.syno.NoSyAll.trmem.gz) <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.All.filt.syno.NoSyAll.trmem

gzip Variant_Matrix.xls.All.filt.syno.NoSyAll.trmem

#Rscript diff_sum_homo.r Variant_Matrix.xls.All.filt.syno.gz 

Rscript diff_sum_heta.r Variant_Matrix.xls.All.filt.syno.NoSyAll.trmem.gz


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.All.filt.syno.NoSyAll.trmemRev

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.syno.NoSyAll.trmemRev.gz) <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.All.filt.syno.NoSyAll.trmemRev

gzip Variant_Matrix.xls.All.filt.syno.NoSyAll.trmemRev

#Rscript diff_sum_homo.r Variant_Matrix.xls.All.filt.syno.gz 

Rscript diff_sum_heta.r Variant_Matrix.xls.All.filt.syno.NoSyAll.trmemRev.gz


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.All.filt.frShift

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.syno.frShift.gz) <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.All.filt.frShift

gzip Variant_Matrix.xls.All.filt.frShift

awk '(NR==FNR){tt[$1]=1;next}{if(tt[$1]==""){print $1}}' <(zcat Variant_Matrix.xls.All.filt.frShift.gz) <(zcat /sc/arion/projects/GOCAR/Sun/5.DRWork/2.impute/Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.syno.frShift.gz) > tt

#Rscript diff_sum_homo.r Variant_Matrix.xls.All.filt.frShift.gz

Rscript diff_sum_heta.r Variant_Matrix.xls.All.filt.frShift.gz


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.All.filt.NoSy

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.syno.NoSy.gz) <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.All.filt.NoSy

gzip Variant_Matrix.xls.All.filt.NoSy

#Rscript diff_sum_homo.r Variant_Matrix.xls.All.filt.NoSy.gz

Rscript diff_sum_heta.r Variant_Matrix.xls.All.filt.NoSy.gz


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.All.filt.SyOther

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.syno.SyOther.gz) <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.All.filt.SyOther

gzip Variant_Matrix.xls.All.filt.SyOther

#Rscript diff_sum_homo.r Variant_Matrix.xls.All.filt.SyOther.gz

Rscript diff_sum_heta.r Variant_Matrix.xls.All.filt.SyOther.gz


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.All.filt.Stop

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.syno.Stop.gz) <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.All.filt.Stop

gzip Variant_Matrix.xls.All.filt.Stop

#Rscript diff_sum_homo.r Variant_Matrix.xls.All.filt.Stop.gz

Rscript diff_sum_heta.r Variant_Matrix.xls.All.filt.Stop.gz


## Gene Level Mismatch

zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.LIMS1

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.LIMS1 <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.LIMS1

Rscript diff_sum_homo.r Variant_Matrix.xls.LIMS1

Rscript diff_sum_heta.r Variant_Matrix.xls.LIMS1

awk '(NR==FNR){tt[$1]=1;next}{if(tt[$1]!=""){print $0}}' Variant_Matrix.xls.LIMS1.comm.heta ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.LIMS1 > Variant_Matrix.xls.LIMS1.comm.heta.ann


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.GCC2

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.GCC2 <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.GCC2

Rscript diff_sum_homo.r Variant_Matrix.xls.GCC2

Rscript diff_sum_heta.r Variant_Matrix.xls.GCC2

awk '(NR==FNR){tt[$1]=1;next}{if(tt[$1]!=""){print $0}}' Variant_Matrix.xls.GCC2.comm.heta ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.GCC2 > Variant_Matrix.xls.GCC2.comm.heta.ann


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.APOL1

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.APOL1 <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.APOL1

Rscript diff_sum_homo.r Variant_Matrix.xls.APOL1

Rscript diff_sum_heta.r Variant_Matrix.xls.APOL1


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.APOL1.exon

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.APOL1.exon <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.APOL1.exon

Rscript diff_sum_homo.r Variant_Matrix.xls.APOL1.exon

Rscript diff_sum_heta.r Variant_Matrix.xls.APOL1.exon


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.NATD1

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.NATD1 <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.NATD1

Rscript diff_sum_homo.r Variant_Matrix.xls.NATD1

Rscript diff_sum_heta.r Variant_Matrix.xls.NATD1


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.rs893403

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.rs893403 <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.rs893403

Rscript diff_sum_homo.r Variant_Matrix.xls.rs893403

Rscript diff_sum_heta.r Variant_Matrix.xls.rs893403


zcat Variant_Matrix.xls.All.filt.gz | head -1 > Variant_Matrix.xls.rs10205370

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.rs10205370 <(zcat Variant_Matrix.xls.All.filt.gz) >> Variant_Matrix.xls.rs10205370

Rscript diff_sum_homo.r Variant_Matrix.xls.rs10205370

Rscript diff_sum_heta.r Variant_Matrix.xls.rs10205370


head -1 Variant_Matrix.xls > Variant_Matrix.xls.Genes33

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.Genes33 Variant_Matrix.xls >> Variant_Matrix.xls.Genes33

Rscript diff_sum_homo.r Variant_Matrix.xls.Genes33


awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR == FNR){tt[$1] = $7"\t"$8; next} {if(tt[$1] != ""){printf("%s\t%s\n", tt[$1], $0)}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.Genes33 Variant_Matrix.xls.Genes33.mutfreq.xls > Variant_Matrix.xls.Genes33.mutfreq.xls.ann.xls

Rscript plot.r


#Rscript diff_sum_homo.r Variant_Matrix.xls.LIMS.ud1M


head -1 Variant_Matrix.xls > Variant_Matrix.xls.TF

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.TF Variant_Matrix.xls >> Variant_Matrix.xls.TF

Rscript diff_sum_homo.r Variant_Matrix.xls.TF


head -1 Variant_Matrix.xls > Variant_Matrix.xls.Kidney

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.Kidney Variant_Matrix.xls >> Variant_Matrix.xls.Kidney

Rscript diff_sum_homo.r Variant_Matrix.xls.Kidney


head -1 Variant_Matrix.xls > Variant_Matrix.xls.Whole_Blood

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.Whole_Blood Variant_Matrix.xls >> Variant_Matrix.xls.Whole_Blood

Rscript diff_sum_homo.r Variant_Matrix.xls.Whole_Blood



head -1 Variant_Matrix.xls > Variant_Matrix.xls.missense.transmem

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.missense.transmem Variant_Matrix.xls >> Variant_Matrix.xls.missense.transmem

Rscript diff_sum_homo.r Variant_Matrix.xls.missense.transmem


head -1 Variant_Matrix.xls > Variant_Matrix.xls.missense

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.missense Variant_Matrix.xls >> Variant_Matrix.xls.missense

Rscript diff_sum_homo.r Variant_Matrix.xls.missense

head -1 Variant_Matrix.xls > Variant_Matrix.xls.all

awk 'BEGIN{FS = "\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$1]!=""){print $0}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID Variant_Matrix.xls >> Variant_Matrix.xls.all

head -1 Variant_Matrix.xls.all > Variant_Matrix.xls.all.freq

cat Split3/*.xls >> Variant_Matrix.xls.all.freq

Rscript diff_sum_homo.r Variant_Matrix.xls.all.freq


#sed -i 's/^X//g' Variant_Matrix.xls

#sed -i 's/./:/g' Variant_Matrix.xls

#awk '{gsub("^X", "", $1); print $0}' Variant_Matrix.xls | le

#awk 'BEGIN{FS="\t"; OFS="\t"} {gsub("\\.", ":", $1); print $0}' Variant_Matrix.xls > tt && mv tt Variant_Matrix.xls

awk 'BEGIN{FS="\t"; OFS = "\t";} (NR == FNR){tt[$1] = $8; next} {if(tt[$1] != ""){print($1, tt[$1], $0)}}' <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.gz) <(zcat Variant_Matrix.xls.All.filt.comm.gz) | cut -f 2,4-388 > Variant_Matrix.xls.Gene.Ori

awk 'BEGIN{FS="\t"; OFS = "\t";} (NR == FNR){tt[$1] = $8; next} {if(tt[$1] != ""){print($1, tt[$1], $0)}}' <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.gz) <(zcat Variant_Matrix.xls.All.filt.comm.gz) | cut -f 1,2 > Variant_Matrix.xls.Gene.Ori.map

awk 'BEGIN{FS="\t"; OFS = "\t";} (NR == FNR){tt[$1] = 1; next} {if(tt[$1] != ""){print($1, $2, $3, $4, $5, $6, $7, $8 )}}' <(zcat Variant_Matrix.xls.All.filt.comm.gz) <(zcat ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.gz) > Variant_Matrix.xls.Gene.Ori.map

zcat Variant_Matrix.xls.Gene.Ori.gz | cut -f 1 | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $1,$2}'> SNP.count

#cat Split2/*.Gene > Variant_Matrix.xls.Gene.1
#
#module load R/3.4.3 && Rscript ./sum.r Variant_Matrix.xls.Gene.1 
#
#head -1 Variant_Matrix.xls > tt && cat tt Variant_Matrix.xls.Gene.1.Gene > Variant_Matrix.xls.Gene


awk 'BEGIN{FS="\t"; OFS = "\t";} (NR == FNR){tt[$1] = $8; next} {if(tt[$1] != ""){print($1, tt[$1], $0)}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.exonic.noSym Variant_Matrix.xls | cut -f 2,4-388 > Variant_Matrix.xls.Gene.exonic.noSym.Ori

module load R/3.4.3 && Rscript ./sum_homo.r Variant_Matrix.xls.Gene.exonic.noSym.Ori

head -1 Variant_Matrix.xls > tt && cat tt Variant_Matrix.xls.Gene.exonic.noSym.Ori.Gene > ttt && mv ttt Variant_Matrix.xls.Gene.exonic.noSym.Ori.Gene


awk 'BEGIN{FS="\t"; OFS = "\t";} (NR == FNR){tt[$1] = $8; next} {if(tt[$1] != ""){print($1, tt[$1], $0)}}' ../Ann/All.tped.out.2.step2.vcf.ann.hg19_multianno.txt.withID.noExonic Variant_Matrix.xls | cut -f 2,4-388 > Variant_Matrix.xls.Gene.noExonic.Ori

#module load R/3.4.3 && Rscript ./sum_homo.r Variant_Matrix.xls.Gene.noExonic.Ori

cat Split4/*.Gene > Variant_Matrix.xls.Gene.noExonic.1

module load R/3.4.3 && Rscript ./sum.r Variant_Matrix.xls.Gene.noExonic.1

head -1 Variant_Matrix.xls > tt && cat tt Variant_Matrix.xls.Gene.noExonic.1.Gene > Variant_Matrix.xls.Gene.noExonic.Ori.Gene






cut -f 1 Variant_Matrix.xls.Gene.Ori | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $1,$2}'> SNP.count

cut -f 3 ../../1.SNP/ann/humandb/hg19_refGene.txt | sort | uniq > Chrs

#awk 'BEGIN{FS="\t"; OFS = "\t"} (NR==FNR){tt[$1] = 1; next} {if(tt[$3]!=""){print $0}}' Chrs ../../1.SNP/ann/humandb/hg19_refGene.txt > hg19_refGene.txt.clean

#cut -f 3,13 hg19_refGene.txt.clean | sort | uniq | cut -f 2 | sort | uniq -d

perl cal_len.pl ../../1.SNP/ann/humandb/hg19_refGene.txt > Gene_length

