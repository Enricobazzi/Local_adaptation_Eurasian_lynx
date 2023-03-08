---
title: "Preparing_Genetic_Data"
author: "Enrico"
date: "2022-10-17"
output: html_document
---

### 1. Raw Reads Preprocessing 

Raw reads were trimmed using the software [Trimmomatic](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu170) in different scripts using the following structure:

```{bash}
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE \
  -threads $THREADS \
  ${INfastq_PATH}/${INfastqR1} ${INfastq_PATH}/${INfastqR2} \
  ${OUT}/${OUTfastqPREFIX}_s1_pe ${OUT}/${OUTfastqPREFIX}_s1_se \
  ${OUT}/${OUTfastqPREFIX}_s2_pe ${OUT}/${OUTfastqPREFIX}_s2_se \
  ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
```

Quality control of trimmed reads was conducted using the software [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)

### 2. Alignment of Reads to Reference Genome

Alignment of trimmed reads to a Canada lynx [reference genome](http://www.nature.com/articles/s41586-021-03451-0) was conducted with multiple script using the softwares [BWA](http://arxiv.org/abs/1303.3997), [samtools](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btp352), [Picard tools](http://broadinstitute.github.io/picard), and [GATK](https://genome.cshlp.org/content/20/9/1297) whose general structure was the following:

```{bash}
# Mapping
echo " - Mapping ${fastqid} of ${SAMPLE} -"
bwa mem ${REF} \
  ${INfastq_PATH}/${fastqid}_1.fastq.gz \
  ${INfastq_PATH}/${fastqid}_2.fastq.gz \
  -t ${THREADS} | 
  samtools view -hbS -@ ${THREADS} - \
  -o ${OUT}/${SAMPLE}_${fastqid}.LyCa_ref.bam

# Sorting
echo " - Sorting ${SAMPLE} -"
samtools sort -@ $THREADS $OUT/${SAMPLE}_${fastqid}.LyCa_ref.bam \
  -o ${OUT}/${SAMPLE}_${fastqid}.LyCa_ref.sorted.bam \
  && rm ${OUT}/${SAMPLE}_${fastqid}.LyCa_ref.bam

# Adding READ Groups
echo " - Adding READ Groups of ${fastqid} of ${SAMPLE} -"
java -jar ${EBROOTPICARD}/picard.jar AddOrReplaceReadGroups \
  I=${OUT}/${SAMPLE}_${fastqid}.LyCa_ref.sorted.bam \
  O=${OUT}/${SAMPLE}_${fastqid}_LyCa_ref_sorted_rg.bam \
  RGID=${fastqid} RGLB=${SAMPLE}_lib \
  RGPL=Illumina RGPU=${fastqid} RGSM=${SAMPLE} \
  VALIDATION_STRINGENCY=SILENT \
  && rm $OUT/${SAMPLE}_${fastqid}.LyCa_ref.sorted.bam

echo " - Merging all ${SAMPLE} BAMs and Re-Sorting -"

ls ${OUT}/${SAMPLE}_*_sorted_rg.bam  > ${OUT}/${SAMPLE}.bam.list

if [[ $(wc -l < ${OUT}/${SAMPLE}.bam.list) -ge 2 ]]
 # if sample has more than one bam -> merge all sample's bams and sort again
 then

  samtools merge  -@ $THREADS \
  -r ${OUT}/${SAMPLE}_merged.bam \
  -b ${OUT}/${SAMPLE}.bam.list

  # BAMARRAY=($(cat ${OUT}/${SAMPLE}.bam.list))
  # for k in ${BAMARRAY[@]}
  #  do
  #   echo " - Removing ${k} -"
  #   rm ${k}
  # done

  echo " - Re-Sorting ${SAMPLE} -"
  samtools sort  -@ $THREADS $OUT/${SAMPLE}_merged.bam \
  -o $OUT/${SAMPLE}_LyCa_ref_sorted_rg.bam \
  && rm $OUT/${SAMPLE}_merged.bam;

 else
 
 # if sample has only one bam -> change name
  mv $(cat $OUT/${SAMPLE}.bam.list) $OUT/${SAMPLE}_LyCa_ref_sorted_rg.bam

fi

# Marking Duplicates
echo " - Marking Duplicates of $fastqid of ${SAMPLE} -"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  METRICS_FILE=$OUT/${SAMPLE}_rmdup.txt \
  I=$OUT/${SAMPLE}_LyCa_ref_sorted_rg.bam \
  O=$OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup.bam \
  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800

rm $OUT/${SAMPLE}_LyCa_ref_sorted_rg.bam

# Re-Sorting
echo " - Re-Sorting ${SAMPLE} -"
samtools sort $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup.bam \
  -@ $THREADS \
  -o $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted.bam

rm $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup.bam

# Indexing for GATK
echo " - Indexing ${SAMPLE} -"
samtools index $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted.bam

# Realigning:
# RealignerTargetCreator
echo " - Realigner Target Creator on ${SAMPLE} -"
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator \
  -nt $THREADS -R $REF \
  -I $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted.bam \
  -o $OUT/${SAMPLE}_realignertargetcreator.intervals

# IndelRealigner
echo " - Realigning INDELS of ${SAMPLE} -"
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner \
  -R $REF \
  -targetIntervals $OUT/${SAMPLE}_realignertargetcreator.intervals \
  -I $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted.bam \
  -o $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted_indelrealigner.bam

rm $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted.bam

# Get stats
echo " - Getting stats of ${SAMPLE} -"
samtools flagstat \
  $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
  > $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted_indelrealigner.stats

```

### 3. Variant Calling

Variant calling was carried out using the software [GATK v4.1.4.1](https://genome.cshlp.org/content/20/9/1297). To perform variant calling in a more efficient manner the reference genome was split into different regions, each of which was called using a script including the the following command:

```{bash}
# INbamARRAY = list of bams to be included in calling
# BED = bed of the region of which to do the calling
# REF = reference genome
# OUTvcf = output vcf file

echo "HaplotypeCaller of region ${bed}"
gatk HaplotypeCaller \
     -R $REF \
     $(for bam in ${INbamARRAY[@]}; do echo "-I ${bam}";done) \
     -L $BED \
     -O $OUTvcf
```

### 4. Variant Filtering

Filtering of variants based on their quality and characteristics was conducted using a custom [script](code/variant_filtering_1to5.sh) that used the softwares [GATK v4.1.1.4](https://genome.cshlp.org/content/20/9/1297) [BED tools v2.28.0](https://academic.oup.com/bioinformatics/article/26/6/841/244688) and
[BCFtools 1.9](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giab008/6137722). Additional filtering based on missingness was conducted with another custom [script](code/variant_filtering_missing.sh).

In order to filter variants based on read depth, removing SNPs with either too few or too many reads aligned, was conducted as a three step process.

(1) Mean read depth along the genome was calculated from a random sample of 200 windows stretching 100 kbp using bedtools random and samtools depth using a custom [script](code/depth_of_random_200_w_subsample.sh).

(2) Read depth limits are calculated using an [R script](code/depth_calculate_limits.R) and saved in a table.

(3) Filtering based on the limits set in the table is applied using another custom [script](code/variant_filtering_depth.sh).

Variants from desired samples (based on the population of origin - see manuscript) and with a minimum MAF of 0.05 are then extracted using the following commands:

```{bash}
cd /home/ebazzicalupo/BayPass/VCF

popARRAY=($(grep -m1 "#CHROM" ll_wholegenome_LyCa_ref.sorted.filter7.vcf \
  | tr '\t' '\n' | grep "c_ll" | cut -d"_" -f3 | sort -u \
  | grep -vE "ba|og|no|po|cr"))
samplesARRAY=($(for pop in ${popARRAY[@]}; do grep -m1 "#CHROM" \
  ll_wholegenome_LyCa_ref.sorted.filter7.vcf | tr '\t' '\n' \
  | grep "c_ll_${pop}" | grep -vE "c_ll_vl_0137|c_ll_tu_0154"; done))

/opt/gatk-4.1.0.0/gatk IndexFeatureFile \
  -F ll_wholegenome_LyCa_ref.sorted.filter7.vcf

/opt/gatk-4.1.0.0/gatk SelectVariants \
  -R /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.fa \
  -V ll_wholegenome_LyCa_ref.sorted.filter7.vcf \
  $(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
  -L /GRUPOS/grupolince/reference_genomes/lynx_canadensis/autosomic_scaffolds.bed \
  -O ll_wholegenome_LyCa_ref.sorted.filter7.finalset.vcf

bcftools view -i 'MAF>0.05' \
  ll_wholegenome_LyCa_ref.sorted.filter7.finalset.vcf \
  > ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.vcf
```

### 5. Variant Pruning for RDA

Variants were pruned and raw files of all SNPs and neutral SNPs were generated using the software [PLINK](https://linkinghub.elsevier.com/retrieve/pii/S0002929707613524) with a custom [script](code/prune_and_convert_to_raw_plink.sh)
