#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// IEDB installation
iedb_chck_file_name = ".iedb_install_ok.chck"
iedb_chck_file = file("${params.outdir}iedb"+ "/" + iedb_chck_file_name)

    process 'install_IEDB' {

        input:
        val(iedb_MHCI_url)
        val(iedb_MHCII_url)

        output:
        path("${iedb_chck_file_name}")

        if(!iedb_chck_file.exists() || iedb_chck_file.isEmpty()) {

          log.warn "WARNING: IEDB yet not installed, starting installation. This may take a while..."

          script:
          def mhci_file = iedb_MHCI_url.split("/")[-1]
          def mhcii_file = iedb_MHCII_url.split("/")[-1]
          """
          CWD=`pwd`
          cd /opt/iedb/
          rm -f $mhci_file
          wget $iedb_MHCI_url
          tar -xzvf $mhci_file
          cd mhc_i
          bash -c "./configure"
          cd /opt/iedb/
          rm -f $mhci_file

          rm -f $mhcii_file
          wget $iedb_MHCII_url
          tar -xzvf $mhcii_file
          cd mhc_ii
          bash -c "python ./configure.py"
          cd /opt/iedb/
          rm $mhcii_file

          cd \$CWD
          echo "OK" > ${iedb_chck_file_name}
          """

        } else {
            script:
            """
            echo "OK" > "${iedb_chck_file_name}"
            """
        }
}

process 'pVACbind_class_I' {

  errorStrategy 'ignore'

  input: 
  tuple val(meta), path(peptides_I), path(optiTypeOutput), path(reads), path(HLA_types_I), path(HLA_types_II)
  path(iedb)

  output:
  path("test.txt")
    
  script:
    """
    pvacbind run --iedb-install-directory /opt/iedb $peptides_fasta $meta.ID $hla_alleles "NetMHCpan" . 
    """
}

process 'pVACbind_class_II' {

  errorStrategy 'ignore'

  input: 
  tuple val(meta), path (peptides_II), path(HLAHD_Output), path(reads), path(HLA_types_I), path(HLA_types_II)
  path(iedb)

  output:
  path("test.txt")
    
  script:
    """
    pvacbind run --iedb-install-directory /opt/iedb $peptides_fasta $meta.ID $hla_alleles "NetMHCIIpan" . 
    """
}

process 'Indices' {

  errorStrategy 'ignore'

  input: 
  path(genome)
  val(hisat_index)
  val(star_index)

  output:
  path("hisat_genome_index/"), optional: true
  path("star_index/"), optional: true
  val("${params.outdir}Indices/hisat_genome_index/hisat_index")
  val("${params.outdir}Indices/star_genome_index/star_index")
    
  script:
    """
    if [ -z $hisat_index ]; then
      mkdir hisat_genome_index
      hisat2-build -p ${task.cpus} $genome hisat_index > hisat_index_log
      mv hisat_index* hisat_genome_index
    fi

    if [ -z $star_index ]; then
      /usr/local/bioinf/rna-star/STAR-2.7.9a/STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles $genome   
    fi
    """
}

process 'OptiType' {

  errorStrategy 'ignore'

  input: 
  tuple val(meta), path(fastq), path(HLA_types_I), path(HLA_types_II)

  output:
  tuple val(meta), path ('*.tsv')
    
  script:
  reads = (meta.libType == "PE") ? fastq[0] + " " + fastq[1] : fastq
    """
    if [ $HLA_types_I == "HLA_default.txt" ]; then
       yara_mapper --version-check 0 -e 3 -t 8 -f bam /home/ausserh/nc_na/HLA_II_typing/yara_index/hla_reference_rna AK11_CRC01_R1_001.fastq.gz AK11_CRC01_R2_001.fastq.gz | samtools view -@ 8 -h -F 4 -b1 -o "AK11_CRC01_mapped_1.bam"
       python /usr/local/bioinf/optiType/latest/OptiTypePipeline.py -i $reads --rna -o . -v -p $meta.ID
    else
       echo "Not used" > optiType_not_used.tsv
    fi

    """
}

/*
hisat is used to align the reads
--outSAMstrandField intronMotif, needed for Stringtie, requires that the strand of the read can be determined from splice junctions - either from annotations, or from the intron motif.
${zcat} zipped input is automatically detected
If --two-pass is specified on the command line, hisat is run in two pass mode
*/ 

process 'alignment' {
  errorStrategy 'ignore'
  input:
      path (genome)
      val (aligner)
      val (self_hisat_index)
      val (self_star_index)
      val (given_hisat_index)
      val (given_star_index)
      tuple val(meta), path(fastq), val(hla_types), val(hla_types_II)
      val two_pass
      val riboseq

  output:
      tuple val(meta), path ("${meta.ID}_Aligned.sortedByCoord.out.bam"), path ("${meta.ID}_Aligned.sortedByCoord.out.bam.bai")

  script:
  // Check if genomeDir is provided and not empty
    def index_star = (given_star_index != "") ? given_star_index : self_star_index
    def index_hisat = (given_hisat_index != "") ? given_hisat_index : self_hisat_index
    def zcat = (fastq[0].getExtension() == "gz") ? "--readFilesCommand zcat" : ""
    reads_star = (meta.libType == "PE") ? fastq[0] + " " + fastq[1] : fastq
    reads_hisat = (meta.libType == "PE") ? "-1 " + fastq[0] + " -2 " + fastq[1] : fastq
    def two_pass =  two_pass ? "--twopassMode Basic" : ""
    def riboseq =  riboseq ? "--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2" : ""
    """
    if [ $aligner == "hisat" ]; then
      hisat2 -p ${task.cpus} -x $index_hisat $reads_hisat | samtools sort -O BAM > "${meta.ID}_Aligned.sortedByCoord.out.bam"
      samtools index "${meta.ID}_Aligned.sortedByCoord.out.bam"
    else
      /usr/local/bioinf/rna-star/STAR-2.7.9a/STAR \
      --runMode alignReads --runThreadN ${task.cpus} --genomeDir $index_star \
      --readFilesIn $reads_star ${two_pass} --outReadsUnmapped Fastx ${zcat} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif ${riboseq}

      samtools sort -@ 4 -m 4G -o "${meta.ID}_Aligned.sortedByCoord.out.bam" Aligned.out.bam
      samtools index "${meta.ID}_Aligned.sortedByCoord.out.bam"
    fi
    """
}

/*
StringTie
-j minimum junction coverage (default: 1)
-s minimum reads per bp coverage to consider for single-exon transcript
-o output path/file name for the assembled transcripts GTF (default: stdout)
-p number of threads (CPUs) to use (default: 1)
-L long reads processing; also enforces -s 1.5 -g 0 (default:false)
Clean_overlapers.py --BED "${meta.ID}_stringtie_inverted.bed" --Out "${meta.ID}_stringtie_inverted_clean.bed" 
bedtools intersect -s -a "${meta.ID}_stringtie_renamed.gtf" -b $inverted_bed > "${meta.ID}_stringtie_inverted_new.bed"
-f 0.001 -m 500 -a 1 -j 0.5 -t -c 0.1 -g 10 -M 0.8 
rename_stringtie_normal.py Custom script to add sample name to transcript ID
*/

process 'StringTie' {
  errorStrategy 'ignore'

  input:
      tuple val (meta), path (bam), path(bai)
      path reference_gtf
      val(longreads)

  output:
      tuple val (meta), path("${meta.ID}_stringtie_renamed.gtf")
      path ("${meta.ID}_stringtie_renamed.gtf")
      tuple val (meta), path ("${meta.ID}_vaf.tsv")


  script:
  def longreads =  longreads ? "-L" : ""
  """
  stringtie -p ${task.cpus} -o "${meta.ID}_stringtie.gtf" $bam -G $reference_gtf
  python3 bin/rename_stringtie.py --stringtie_out "${meta.ID}_stringtie.gtf" --gtf_out "${meta.ID}_stringtie_renamed.gtf" --vaf_out "${meta.ID}_vaf.tsv"
  """
}

/*
SAMTOOLS:
-h       include header in SAM output
-b       output BAM
-f INT   only include reads with all  of the FLAGs in INT present [0]

  samtools view -bq 60 -@ ${task.cpus} "${meta.ID}_hla.bam" > "${meta.ID}_hla_hq.bam"
  samtools rmdup "${meta.ID}_hla_hq.bam" "${meta.ID}_hla_hq_clean.bam"


BEDTOOLS: samtools rmdup "${meta.ID}_hla.bam" "${meta.ID}_hla_clean.bam"
*/

process 'HLA_extraction' {

  input: 
  tuple val(meta), path (bam), path(bai), path(reads), path(HLA_types_I), path(HLA_types_II)
    
  output:
  tuple val(meta), path("${meta.ID}_hlatmp.1.fastq"), path("${meta.ID}_hlatmp.2.fastq")
    
  script:
    """
  if [ $HLA_types_II == "HLA_II_default.txt" ]; then
    #Extract MHC region 
    samtools view -h -@ ${task.cpus} -f 2 -b $bam chr6:28,510,120-33,480,577 > "${meta.ID}_hla.bam"
    samtools rmdup "${meta.ID}_hla.bam" "${meta.ID}_hla_clean.bam"
    samtools sort -n "${meta.ID}_hla_clean.bam" > "${meta.ID}_hla_sorted.bam"
    bedtools bamtofastq -i "${meta.ID}_hla_sorted.bam" -fq "${meta.ID}_hlatmp.1.fastq" -fq2 "${meta.ID}_hlatmp.2.fastq"
  else
    touch "${meta.ID}_hlatmp.1.fastq"
    touch "${meta.ID}_hlatmp.2.fastq"
  fi
    """
}

/*
Running HLA-HD to type class II genes
[-t thread]
[-m minmum_tag_size]
*/  

process 'HLA_HD' {

  input: 
  tuple val(meta), path(fastq_1), path(fastq_2), path(reads), path(HLA_types_I), path(HLA_types_II)
    
  output:
  tuple val(meta), path ('HLA_HD_out/*/result/*final.result.txt')
    
  script:
    """
  if [ $HLA_types_II == "HLA_II_default.txt" ]; then
    export PATH=$PATH:/usr/local/bioinf/HLA_HD/hlahd.1.4.0/bin/ 
    mkdir HLA_HD_out
    hlahd.sh -t ${task.cpus} -m 50 -f \
      /usr/local/bioinf/HLA_HD/hlahd.1.4.0/freq_data/ \
      $fastq_1 \
      $fastq_2 \
      /usr/local/bioinf/HLA_HD/hlahd.1.4.0/HLA_gene.split.txt \
      /usr/local/bioinf/HLA_HD/hlahd.1.4.0/dictionary $meta.ID HLA_HD_out
  else
    mkdir HLA_HD_out
    mkdir HLA_HD_out/test
    mkdir HLA_HD_out/test/result
    touch HLA_HD_out/test/result/test_final.result.txt
  fi      
    """
}

/*
Peptides are cleaved with a sliding window approach
 --window_shift 1
 --window_shift 1
*/ 

process 'Filtering' {
  errorStrategy 'ignore'

  input:
      tuple val (meta), path (peptides_class_I)
      path database_I
      val (split)

  output:
      tuple val(meta), path ("${meta.ID}_peptides_I_filtered.tsv")
      tuple val (meta), path ("${meta.ID}_peptides_MHCI_filtered.fasta.split/${meta.ID}_peptides_MHCI_filtered.part_0*")
      tuple val(meta), path ("${meta.ID}_peptides_I_filtered_paste.fasta")
      tuple val (meta), path ("${meta.ID}_peptides_MHCII_filtered.fasta.split/${meta.ID}_peptides_MHCII_filtered.part_0*"), optional: true
      

  script:
  """
  grep "Differential" $peptides_class_I || true > "${meta.ID}_Differential.fasta"
  if [ -s "${meta.ID}_Differential.fasta" ]; then
        # The file is not-empty.
        seqkit grep -j ${task.cpus} -n -v -f <(seqkit seq -n "${meta.ID}_Differential.fasta") $peptides_class_I > "${meta.ID}_peptides_I_reduced.fasta"
  else
        # The file is empty.
        cat $peptides_class_I > "${meta.ID}_peptides_I_reduced.fasta"
  fi
  
  seqkit common -j ${task.cpus} -s "${meta.ID}_peptides_I_reduced.fasta" $database_I > "${meta.ID}_common_I_peptides.fasta"
  seqkit grep -j ${task.cpus} -n -v -f <(seqkit seq -n "${meta.ID}_common_I_peptides.fasta") $peptides_class_I | seqkit rmdup -j ${task.cpus} > "${meta.ID}_peptides_I_filtered.fasta"
  grep -v ">" "${meta.ID}_peptides_I_filtered.fasta" > "${meta.ID}_peptides_I_filtered_paste.fasta"
  seqkit fx2tab "${meta.ID}_peptides_I_filtered.fasta" > "${meta.ID}_peptides_I_filtered.tsv"

  seqkit seq -m 1 -M 12 "${meta.ID}_peptides_I_filtered.fasta" > "${meta.ID}_peptides_MHCI_filtered.fasta"
  seqkit seq -m 13 -M 50 "${meta.ID}_peptides_I_filtered.fasta" > "${meta.ID}_peptides_MHCII_filtered.fasta"

  seqkit split2 --by-part $split "${meta.ID}_peptides_MHCI_filtered.fasta"
  seqkit split2 --by-part $split "${meta.ID}_peptides_MHCII_filtered.fasta"
  """
}

process 'Create_capture_bed' {

  input:
      path (GTF)
      val (tpm_max_diff)
      val (cov_max_diff)

  output:
      tuple path ("capture_bed.bed"), path ("All.gtf"), path ("capture_bed_differential.bed"), path ("capture_bed_novel.bed")

  script:
  """
  cat *.gtf > All.gtf

  # Novel capture bed

  # Exon
  cat All.gtf | awk 'BEGIN{OFS="\t";} (\$3=="exon" || \$3=="CDS") {print \$1,\$4-1,\$5,\$8,\$9,\$7}' | \
  bedtools sort | bedtools merge -s -c 6 -o distinct -i - | awk -v OFS='\t' '{print \$0, "EXON"}' | awk -v OFS='\t' '{print \$0, "1000"}' | awk -v OFS='\t' '{print \$1,\$2,\$3,\$5,\$6,\$4}' > All_exon.bed

  # Intron
  cat All.gtf | awk 'BEGIN{OFS="\t";} \$3=="transcript" {print \$1,\$4-1,\$5,\$8,\$9,\$7}' | bedtools sort | bedtools subtract -s -a stdin -b All_exon.bed | \
  bedtools sort | bedtools merge -s -c 6 -o distinct -i - | awk -v OFS='\t' '{print \$0, "INTRON"}' | awk -v OFS='\t' '{print \$0, "1000"}' | awk -v OFS='\t' '{print \$1,\$2,\$3,\$5,\$6,\$4}' > All_intron.bed

  # Intergenic
  grep "chr" All.gtf | grep -w "+" | grep -v "_alt" | grep -v "_random" | grep -v "chrUn" | bedtools sort -i - | \
  bedtools complement -i stdin -g /home/ausserh/myScratch/Annotation/chr_lengths_clean.genome | awk -v OFS='\t' '{print \$0, "INTERGENIC"}' | awk -v OFS='\t' '{print \$0, "1000"}' | awk -v OFS='\t' '{print \$0, "+"}' > All_intergenic_plus.bed
  
  grep "chr" All.gtf | grep -w "-" | grep -v "_alt" | grep -v "_random" | grep -v "chrUn" | bedtools sort -i - | \
  bedtools complement -i stdin -g /home/ausserh/myScratch/Annotation/chr_lengths_clean.genome | awk -v OFS='\t' '{print \$0, "INTERGENIC"}' | awk -v OFS='\t' '{print \$0, "1000"}' | awk -v OFS='\t' '{print \$0, "-"}' > All_intergenic_minus.bed

  cat All_intergenic_plus.bed All_intergenic_minus.bed > All_intergenic.bed

  bedtools subtract -s -a All_intergenic.bed -b All.gtf > All_intergenic_2.bed
  bedtools subtract -s -a All_intergenic_2.bed -b All_intron.bed > All_intergenic_3.bed
  
  cat All_intron.bed All_intergenic_3.bed > capture_bed_novel.bed

  # Differential capture bed

  capture_bed_differential.py --GTF All.gtf --Out All_capture_differential_1.bed

  cat All_capture_differential_1.bed | bedtools sort | bedtools merge -s -c 4,4,6,5,5 -o mean,max,distinct,mean,max -d 50 -i stdin > All_capture_differential_1_merged.bed

  capture_bed_differential_2.py --BED All_capture_differential_1_merged.bed --tpm_max_diff $tpm_max_diff --cov_max_diff $cov_max_diff --Out capture_bed_differential.bed

  cat capture_bed_novel.bed capture_bed_differential.bed > capture_bed.bed
  """
}


/*
grep -w -f IDs.txt $stringtie_gtf > "${meta.ID}_stringtie_short.gtf"
*/ 

process 'Translation' {
  errorStrategy 'ignore'

  input:
      tuple val (meta), path (stringtie_gtf), path (bed)
      path (genome)
      path (gtf)
      path (ref_prot)
      val peptide_length
    
  output:
      tuple val (meta), path ("${meta.ID}_peptide_regions.split.bed0*")
      tuple val (meta), path ("${meta.ID}_stringtie_short.gtf")
      tuple val (meta), path ("${meta.ID}_peptides_class_sorted.fasta")
      tuple val (meta), path ("${meta.ID}_peptide_regions_short_sorted.fasta")
      tuple val (meta), path ("${meta.ID}_regions.fasta")
      tuple val (meta), path ("${meta.ID}_stringtie.tsv")
      tuple val (meta), path ("${meta.ID}_stringtie.fasta")
      tuple val (meta), path ("${meta.ID}_peptide_regions.bed")
      tuple val (meta), path ("${meta.ID}_overlaps.tsv")
      tuple val (meta), path ("${meta.ID}_closest_reference.tsv")
      tuple val (meta), path ("${meta.ID}_translated.tsv")

  script:
  """
  cut -f 4 $bed | sed 's![^_]*\$!!' | sed 's/.\$//' | tail -n +2 > IDs.txt
  Cancer_transcripts.py --GTF $stringtie_gtf --ID IDs.txt --Out "${meta.ID}_stringtie_short_0.gtf"
  grep "chr" "${meta.ID}_stringtie_short_0.gtf" | grep -v "_alt" | grep -v "_random" | grep -v "chrUn" > "${meta.ID}_stringtie_short.gtf"
  grep "chr" $bed | grep -v "_alt" | grep -v "_random" | grep -v "chrUn" > "${meta.ID}_overlap_clean.bed"
  /home/ausserh/nc_na/gffread/gffread/./gffread -w "${meta.ID}_stringtie.fasta" -g $genome "${meta.ID}_stringtie_short.gtf"
  /home/ausserh/nc_na/gffread/gffread/./gffread -w "${meta.ID}_regions.fasta" -g $genome "${meta.ID}_overlap_clean.bed"

  seqkit fx2tab $ref_prot > gencode.v38.pc_translations.tsv
  seqkit fx2tab "${meta.ID}_stringtie.fasta" > "${meta.ID}_stringtie.tsv"

  translation_module_4.py --Short_GTF "${meta.ID}_stringtie_short.gtf" --Stringtie_tsv "${meta.ID}_stringtie.tsv" \
  --Gencode_GTF /data/genomes/hg38/annotation/gencode/gencode.v38.primary_assembly.annotation.gtf --Reference_tsv gencode.v38.pc_translations.tsv \
  --closest_out "${meta.ID}_closest_reference.tsv" --matching_out "${meta.ID}_translated.tsv" --Regions_bed "${meta.ID}_overlap_clean.bed" --peptide_length $peptide_length --Bed_out "${meta.ID}_peptide_regions.bed" --Overlaps_out "${meta.ID}_overlaps.tsv" --Regions_fasta "${meta.ID}_regions.fasta" --Peptides_out "${meta.ID}_peptides.fasta"
  
  seqkit sort -n -i "${meta.ID}_peptides.fasta" > "${meta.ID}_peptides_class_sorted.fasta"

  cut -f1,2,3,4,6,7 "${meta.ID}_peptide_regions.bed" > "${meta.ID}_peptide_regions_short.bed"
  /home/ausserh/nc_na/gffread/gffread/./gffread -w "${meta.ID}_peptide_regions_short.fasta" -g $genome "${meta.ID}_peptide_regions_short.bed"
  seqkit sort -n -i "${meta.ID}_peptide_regions_short.fasta" > "${meta.ID}_peptide_regions_short_sorted.fasta"

  split -l\$((`wc -l < "${meta.ID}_peptide_regions.bed"`/16)) "${meta.ID}_peptide_regions.bed" "${meta.ID}_peptide_regions.split.bed" -da 4
  """
}

process 'Annotation' {
  errorStrategy 'ignore'

  input:
      tuple val (meta), path (GTF)
      val (tpm_min_novel)
      val (cov_min_novel)
      val (tpm_min_diff)
      val (cov_min_diff)
      path(inverted_bed)

  output:
      tuple val (meta), path ("${meta.ID}_overlap.bed")
      tuple val (meta), path ("${meta.ID}_annotation.bed")
 

  script:
  """
  bedtools intersect -wo -s -a $GTF -b $inverted_bed > "${meta.ID}_stringtie_specific_new.bed"
  Clean_overlapers_2.5.py --BED "${meta.ID}_stringtie_specific_new.bed" --GTF $GTF --TPM_min_novel $tpm_min_novel --Cov_min_novel $cov_min_novel --Out_anno "${meta.ID}_annotation.bed" --Out_bed "${meta.ID}_overlap.bed" --TPM_min_diff $tpm_min_diff --Cov_min_diff $cov_min_diff
  """
}

process 'Annotation_2' {
  errorStrategy 'ignore'

  input:
      tuple val (meta), path (bed), path(bam), path(bai)
      val (aligner)
      val (BAM_cov)

  output:
      tuple val (meta), path ("${meta.ID}_${bed.extension}_selected_sequences.fasta")
      tuple val (meta), path ("${meta.ID}_${bed.extension}_coverage_reads.tsv")
      tuple val (meta), path ("${meta.ID}_${bed.extension}_selected_sequences.bed")
 

  script:
  def quality = (aligner == "hisat") ? 60 : 255
  """
  mv $bed "${meta.ID}_regions.bed"
  samtools view -bq $quality -@ ${task.cpus} -L "${meta.ID}_regions.bed" $bam > "${meta.ID}_survivor_coverage_1.bam"
  samtools rmdup "${meta.ID}_survivor_coverage_1.bam" "${meta.ID}_survivor_coverage.bam"
  java -jar /home/ausserh/myScratch/jvarkit/jvarkit/dist/sortsamrefname.jar --samoutputformat BAM "${meta.ID}_survivor_coverage.bam" | java -jar /home/ausserh/myScratch/jvarkit/jvarkit/dist/biostar154220.jar -n 100 --samoutputformat BAM > "${meta.ID}_survivor_coverage_limited.bam"
  samtools sort "${meta.ID}_survivor_coverage_limited.bam" > "${meta.ID}_survivor_coverage_limited_sorted.bam" 
  samtools index "${meta.ID}_survivor_coverage_limited_sorted.bam" 
  java -jar /home/ausserh/myScratch/jvarkit/dist/sam4weblogo.jar -F tabular -r "${meta.ID}_regions.bed" "${meta.ID}_survivor_coverage_limited_sorted.bam" | grep -v '>' | grep -v "\\--" > "${meta.ID}_${bed.extension}_coverage_reads.tsv"
  process_subreads.py --tsv_file "${meta.ID}_${bed.extension}_coverage_reads.tsv" --tsv_out "${meta.ID}_${bed.extension}_selected_sequences.bed" --BAM_cov $BAM_cov
  awk '{ printf ">%s\\n%s\\n",\$4"_"\$7,\$5 }' "${meta.ID}_${bed.extension}_selected_sequences.bed" | seqkit sort -n -i > "${meta.ID}_${bed.extension}_selected_sequences.fasta"
  """
}

process 'Combine' {

  input:
      tuple val (meta), path (selected_fasta), path (bed)

  output:
      tuple val (meta), path ("${meta.ID}_selected.fasta")
      tuple val (meta), path ("${meta.ID}_selected.bed") 
 

  script:
  """
  cat *.fasta > "${meta.ID}_selected.fasta"
  cat *.bed > "${meta.ID}_selected.bed"
  """
}

process 'Translation_2' {
  errorStrategy 'ignore' 

  input:
      tuple val (meta), path (selected), path (peptides), path (regions)

  output:
      tuple val (meta), path ("${meta.ID}_patient_specific_ncnas.tsv")
      tuple val (meta), path ("${meta.ID}_patient_specific_ncnas_meta.tsv")

  script:
  """
  translate_with_snps.py --fasta_in $selected --peptides $peptides --regions $regions --tsv_out "${meta.ID}_patient_specific_ncnas.tsv" --tsv_out_meta "${meta.ID}_patient_specific_ncnas_meta.tsv"
  """
}

process 'netMHCpan' {
  errorStrategy 'ignore' 

  input:
  tuple val(meta), path(peptides_I), path(optiTypeOutput), path(reads), path(HLA_types_I), path(HLA_types_II)

    
  output:
  tuple val(meta), path("${peptides_I.baseName}_peptides_I_binding.tsv")
    
  script:
    """
    mkdir netmhcpan_tmp
    export TMPDIR="netmhcpan_tmp"
    paste <(awk 'NR % 2 == 0' $peptides_I) <(awk 'NR % 2 == 1' $peptides_I) > "${meta.ID}_peptides_I_pasted.fasta"
    if [ $HLA_types_I == "HLA_default.txt" ]; then
    # Use OptiType output as input
      Class_I=\$(tr -s '\t'  '\n'<$optiTypeOutput | sed -n 11,16p | grep "*" | awk '{ gsub(/*/,"", \$1); print }' | awk '\$1="HLA-"\$1' | awk '{print \$1}' | paste -s -d, -)
    else
    # Use HLA types file as input
      Class_I=\$(< $HLA_types_I)
    fi
    /usr/local/bioinf/bin/netMHCpan -a \$Class_I -p -inptype 1 -BA -f "${meta.ID}_peptides_I_pasted.fasta" -xls -xlsfile "${peptides_I.baseName}_peptides_I_binding.tsv"
    """
}

/*
Using netMHCIIpan
*/ 

process 'netMHCIIpan' {
  errorStrategy 'ignore'

  input: 
  tuple val(meta), path (peptides_II), path(HLAHD_Output), path(reads), path(HLA_types_I), path(HLA_types_II)
    
  output:
  tuple val(meta), path ("${peptides_II.baseName}_peptides_II_binding.tsv")
    
  script:
    """
    mkdir netmhcIIpan_tmp
    export TMPDIR="netmhcIIpan_tmp"
    paste <(awk 'NR % 2 == 0' $peptides_II) <(awk 'NR % 2 == 1' $peptides_II) > "${meta.ID}_peptides_II_pasted.fasta"
    awk -F '\t' '{print \$1}' "${meta.ID}_peptides_II_pasted.fasta" > peptides_short.tsv
    if [ $HLA_types_II == "HLA_II_default.txt" ]; then
      HLA_HD_out_reshape.py --HLA_HD_out $HLAHD_Output --file_out HLA_HD_out_reshape.csv
      Class_II=\$(awk '{print \$1}' HLA_HD_out_reshape.csv |  paste -s -d, -)
      echo \$Class_II > Test.txt
    else
    # Use HLA types file as input
      Class_II=\$(< $HLA_types_II)
    fi
    /usr/local/bioinf/bin/netMHCIIpan -a \$Class_II -inptype 1 -tdir netmhcIIpan_tmp -BA -f peptides_short.tsv -xls -xlsfile "${peptides_II.baseName}_peptides_II_binding.tsv"
    """
}

/*
Creating final output table  \$Class_II
*/ 

process 'final_out_1' {
  errorStrategy 'ignore'

  input: 
  tuple val(meta), path(peptides_I_binding)
    
  output:
  tuple val(meta), path ("${meta.ID}_final_I_out.tsv")
    
  script:
    """
   head -2 *1_peptides_I_binding.tsv > "${meta.ID}_final_I_out.tsv"; tail -n +3 -q *_peptides_I_binding.tsv >> "${meta.ID}_final_I_out.tsv"
    """
}

process 'Metadata_MHCI' {

  input: 
  tuple val(meta), path(filtering), path(vaf), path(bed), path(translation), path(bed_2), path(specific), path(bind)
  path (anno_2)
    
  output:
  tuple val(meta), path ("${meta.ID}_final_out_combined.tsv")
  tuple val(meta), path ("${meta.ID}_score_novel_sb.tsv")
  tuple val(meta), path ("${meta.ID}_score_novel_wb.tsv")
  tuple val(meta), path ("${meta.ID}_score_diff_sb.tsv")
  tuple val(meta), path ("${meta.ID}_score_diff_wb.tsv")
    
  script:
    """
  final_out_6.py --Filtering $filtering --VAF $vaf --BED $bed --Translation $translation --BED_2 $bed_2 --Specific $specific --BIND $bind --Out "${meta.ID}_final_out_combined_1.tsv" --Out_header "${meta.ID}_final_out_combined_0.tsv"
  bedtools intersect -wao -s -a "${meta.ID}_final_out_combined_1.tsv" -b $anno_2 > "${meta.ID}_final_out_combined_2.tsv"
  final_out_5.py --BED "${meta.ID}_final_out_combined_2.tsv" --BED_old "${meta.ID}_final_out_combined_0.tsv" --Out "${meta.ID}_final_out_combined.tsv"
  Final_score.py --Final "${meta.ID}_final_out_combined.tsv" --Out_novel_sb "${meta.ID}_score_novel_sb.tsv" --Out_novel_wb "${meta.ID}_score_novel_wb.tsv" --Out_diff_sb "${meta.ID}_score_diff_sb.tsv" --Out_diff_wb "${meta.ID}_score_diff_wb.tsv"
    """
}

process 'Metadata_MHCII' {

  input: 
  tuple val(meta), path(filtering), path(vaf), path(bed), path(translation), path(bed_2), path(specific), path(bind)
  path (anno_2)
    
  output:
  tuple val(meta), path ("${meta.ID}_final_out_combined.tsv")
  tuple val(meta), path ("${meta.ID}_score_novel_sb.tsv")
  tuple val(meta), path ("${meta.ID}_score_novel_wb.tsv")
  tuple val(meta), path ("${meta.ID}_score_diff_sb.tsv")
  tuple val(meta), path ("${meta.ID}_score_diff_wb.tsv")
    
  script:
    """
  final_out_6.py --Filtering $filtering --VAF $vaf --BED $bed --Translation $translation --BED_2 $bed_2 --Specific $specific --BIND $bind --Out "${meta.ID}_final_out_combined_1.tsv" --Out_header "${meta.ID}_final_out_combined_0.tsv"
  bedtools intersect -wao -s -a "${meta.ID}_final_out_combined_1.tsv" -b $anno_2 > "${meta.ID}_final_out_combined_2.tsv"
  final_out_5.py --BED "${meta.ID}_final_out_combined_2.tsv" --BED_old "${meta.ID}_final_out_combined_0.tsv" --Out "${meta.ID}_final_out_combined.tsv"
  Final_score.py --Final "${meta.ID}_final_out_combined.tsv" --Out_novel_sb "${meta.ID}_score_novel_sb.tsv" --Out_novel_wb "${meta.ID}_score_novel_wb.tsv" --Out_diff_sb "${meta.ID}_score_diff_sb.tsv" --Out_diff_wb "${meta.ID}_score_diff_wb.tsv"
    """
}


process 'final_out_2' {
  errorStrategy 'ignore'

  input: 
  tuple val(meta), path(peptides_II_binding)
    
  output:
  tuple val(meta), path ("${meta.ID}_final_II_out.tsv")
    
  script:
    """
   head -2 *1_peptides_II_binding.tsv > "${meta.ID}_final_II_out.tsv"; tail -n +3 -q *_peptides_II_binding.tsv >> "${meta.ID}_final_II_out.tsv"
    """
}