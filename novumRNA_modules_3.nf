#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// IEDB installation
iedb_chck_file_name = "iedb_install_ok.chck"
iedb_chck_file = file("${params.outdir_ref}/iedb/" + iedb_chck_file_name)

    process 'install_IEDB' {

      input:
      val(iedb_MHCI_url)
      val(iedb_MHCII_url)
      val(iedb)

      output:
      path("${iedb_chck_file_name}")

      script:
      def mhci_file = iedb_MHCI_url.split("/")[-1]
      def mhcii_file = iedb_MHCII_url.split("/")[-1]
      """
      if [ -z $iedb ]; then
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
        echo "OK" > "${iedb_chck_file_name}"
      else
        echo "OK" > "${iedb_chck_file_name}"
      fi
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
      /hisat2-2.1.0/hisat2-build -p ${task.cpus} $genome hisat_index > hisat_index_log
      mv hisat_index* hisat_genome_index
    fi

    if [ -z $star_index ]; then
      /STAR-2.7.9a/source/STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles $genome   
    fi
    """
}

process 'OptiType' {

  errorStrategy 'ignore'

  input: 
  tuple val(meta), path(fastq), path(HLA_types_I), path(HLA_types_II)
  val(HLA_ref)

  output:
  tuple val(meta), path ('*.tsv')
    
  script:
  reads = (meta.libType == "PE") ? fastq[0] + " " + fastq[1] : fastq
    """
    if [ $HLA_types_I == "HLA_default.txt" ]; then
       /yara-build/bin/yara_mapper --version-check 0 -e 3 -t 8 -f bam $HLA_ref $fastq | samtools view -@ 8 -h -F 4 -b1 -o "${meta.ID}_mapped_1.bam"
       samtools sort "${meta.ID}_mapped_1.bam" > "${meta.ID}_mapped_1_sorted.bam"
       samtools index "${meta.ID}_mapped_1_sorted.bam"
       python3 /OptiType-1.3.3/OptiTypePipeline.py -i "${meta.ID}_mapped_1_sorted.bam" --rna -o . -v -p $meta.ID
    else
       echo "Not used" > optiType_not_used.tsv
    fi

    """
}

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
      /hisat2-2.1.0/hisat2 -p ${task.cpus} -x $index_hisat $reads_hisat | samtools sort -O BAM > "${meta.ID}_Aligned.sortedByCoord.out.bam"
      samtools index "${meta.ID}_Aligned.sortedByCoord.out.bam"
    else
      /STAR-2.7.9a/source/STAR \
      --runMode alignReads --runThreadN ${task.cpus} --genomeDir $index_star \
      --readFilesIn $reads_star ${two_pass} --outReadsUnmapped Fastx ${zcat} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif ${riboseq}

      samtools sort -@ ${task.cpus} -m 4G -o "${meta.ID}_Aligned.sortedByCoord.out.bam" Aligned.out.bam
      samtools index "${meta.ID}_Aligned.sortedByCoord.out.bam"
    fi
    """
}

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
  python3 /scripts/rename_stringtie.py --stringtie_out "${meta.ID}_stringtie.gtf" --gtf_out "${meta.ID}_stringtie_renamed.gtf" --vaf_out "${meta.ID}_vaf.tsv"
  """
}

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
    /bedtools2/bin/bedtools bamtofastq -i "${meta.ID}_hla_sorted.bam" -fq "${meta.ID}_hlatmp.1.fastq" -fq2 "${meta.ID}_hlatmp.2.fastq"
  else
    touch "${meta.ID}_hlatmp.1.fastq"
    touch "${meta.ID}_hlatmp.2.fastq"
  fi
    """
}

process 'HLA_HD' {

  input: 
  tuple val(meta), path(fastq_1), path(fastq_2), path(reads), path(HLA_types_I), path(HLA_types_II)

  output:
  tuple val(meta), path ('HLA_HD_out/*/result/*final.result.txt')
    
  script:
    """
  if [ $HLA_types_II == "HLA_II_default.txt" ]; then
    mkdir HLA_HD_out
    hlahd.sh -t ${task.cpus} -m 50 -f \
      /opt/hlahd/freq_data/ \
      $fastq_1 \
      $fastq_2 \
      /opt/hlahd/HLA_gene.split.txt \
      /opt/hlahd/dictionary $meta.ID HLA_HD_out
  else
    mkdir HLA_HD_out
    mkdir HLA_HD_out/test
    mkdir HLA_HD_out/test/result
    touch HLA_HD_out/test/result/test_final.result.txt
  fi      
    """
}

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
      path(genome_length)

  output:
      tuple path ("capture_bed.bed"), path ("All.gtf"), path ("capture_bed_differential.bed"), path ("capture_bed_novel.bed")

  script:
  """
  cat *.gtf > All.gtf

  # Novel capture bed

  # Exon
  cat All.gtf | awk 'BEGIN{OFS="\t";} (\$3=="exon" || \$3=="CDS") {print \$1,\$4-1,\$5,\$8,\$9,\$7}' | \
  /bedtools2/bin/bedtools sort | /bedtools2/bin/bedtools merge -s -c 6 -o distinct -i - | awk -v OFS='\t' '{print \$0, "EXON"}' | awk -v OFS='\t' '{print \$0, "1000"}' | awk -v OFS='\t' '{print \$1,\$2,\$3,\$5,\$6,\$4}' > All_exon.bed

  # Intron
  cat All.gtf | awk 'BEGIN{OFS="\t";} \$3=="transcript" {print \$1,\$4-1,\$5,\$8,\$9,\$7}' | /bedtools2/bin/bedtools sort | /bedtools2/bin/bedtools subtract -s -a stdin -b All_exon.bed | \
  /bedtools2/bin/bedtools sort | /bedtools2/bin/bedtools merge -s -c 6 -o distinct -i - | awk -v OFS='\t' '{print \$0, "INTRON"}' | awk -v OFS='\t' '{print \$0, "1000"}' | awk -v OFS='\t' '{print \$1,\$2,\$3,\$5,\$6,\$4}' > All_intron.bed

  # Intergenic
  grep "chr" All.gtf | grep -w "+" | grep -v "_alt" | grep -v "_random" | grep -v "chrUn" | /bedtools2/bin/bedtools sort -i - | \
  /bedtools2/bin/bedtools complement -i stdin -g $genome_length | awk -v OFS='\t' '{print \$0, "INTERGENIC"}' | awk -v OFS='\t' '{print \$0, "1000"}' | awk -v OFS='\t' '{print \$0, "+"}' > All_intergenic_plus.bed
  
  grep "chr" All.gtf | grep -w "-" | grep -v "_alt" | grep -v "_random" | grep -v "chrUn" | /bedtools2/bin/bedtools sort -i - | \
  /bedtools2/bin/bedtools complement -i stdin -g $genome_length | awk -v OFS='\t' '{print \$0, "INTERGENIC"}' | awk -v OFS='\t' '{print \$0, "1000"}' | awk -v OFS='\t' '{print \$0, "-"}' > All_intergenic_minus.bed

  cat All_intergenic_plus.bed All_intergenic_minus.bed > All_intergenic.bed

  /bedtools2/bin/bedtools subtract -s -a All_intergenic.bed -b All.gtf > All_intergenic_2.bed
  /bedtools2/bin/bedtools subtract -s -a All_intergenic_2.bed -b All_intron.bed > All_intergenic_3.bed
  
  cat All_intron.bed All_intergenic_3.bed > capture_bed_novel.bed

  # Differential capture bed

  /scripts/capture_bed_differential.py --GTF All.gtf --Out All_capture_differential_1.bed

  cat All_capture_differential_1.bed | /bedtools2/bin/bedtools sort | /bedtools2/bin/bedtools merge -s -c 4,4,6,5,5 -o mean,max,distinct,mean,max -d 50 -i stdin > All_capture_differential_1_merged.bed

  /scripts/capture_bed_differential_2.py --BED All_capture_differential_1_merged.bed --tpm_max_diff $tpm_max_diff --cov_max_diff $cov_max_diff --Out capture_bed_differential.bed

  cat capture_bed_novel.bed capture_bed_differential.bed > capture_bed.bed
  """
}

process 'Protein_to_peptides' {
  errorStrategy 'ignore'

  input:
      path (proteins)
      path (ref_pep)
      val length_1

  output:
      path ("Control_peptides_len_${length_1}_rmdup.fasta")

  script:
  """
  if [ $ref_pep == "Test_ref_pep.txt" ]; then
    /scripts/protein_2_peptide.py --fasta_in $proteins --fasta_out "Control_peptides_len_${length_1}.fasta" --pep_len $length_1 --window_shift 1
    seqkit rmdup -s "Control_peptides_len_${length_1}.fasta" > "Control_peptides_len_${length_1}_rmdup.fasta"
    rm "Control_peptides_len_${length_1}.fasta"
  elif [ $ref_pep == "Control_peptides_len_${length_1}_rmdup.fasta" ]; then
    echo "Fine" 
  else
    mv $ref_pep "Control_peptides_len_${length_1}_rmdup.fasta"
  fi  
  """
}


process 'Translation' {
  errorStrategy 'ignore'

  input:
      tuple val (meta), path (stringtie_gtf), path (bed)
      path (genome)
      path (gtf)
      path (ref_prot)
      val peptide_length
      val(split)
    
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
  /scripts/Cancer_transcripts.py --GTF $stringtie_gtf --ID IDs.txt --Out "${meta.ID}_stringtie_short_0.gtf"
  grep "chr" "${meta.ID}_stringtie_short_0.gtf" | grep -v "_alt" | grep -v "_random" | grep -v "chrUn" > "${meta.ID}_stringtie_short.gtf"
  grep "chr" $bed | grep -v "_alt" | grep -v "_random" | grep -v "chrUn" > "${meta.ID}_overlap_clean.bed"
  gffread -w "${meta.ID}_stringtie.fasta" -g $genome "${meta.ID}_stringtie_short.gtf"
  gffread -w "${meta.ID}_regions.fasta" -g $genome "${meta.ID}_overlap_clean.bed"

  seqkit fx2tab $ref_prot > gencode.v38.pc_translations.tsv
  seqkit fx2tab "${meta.ID}_stringtie.fasta" > "${meta.ID}_stringtie.tsv"

  /scripts/translate.py --Short_GTF "${meta.ID}_stringtie_short.gtf" --Stringtie_tsv "${meta.ID}_stringtie.tsv" \
  --Gencode_GTF /data/genomes/hg38/annotation/gencode/gencode.v38.primary_assembly.annotation.gtf --Reference_tsv gencode.v38.pc_translations.tsv \
  --closest_out "${meta.ID}_closest_reference.tsv" --matching_out "${meta.ID}_translated.tsv" --Regions_bed "${meta.ID}_overlap_clean.bed" --peptide_length $peptide_length --Bed_out "${meta.ID}_peptide_regions.bed" --Overlaps_out "${meta.ID}_overlaps.tsv" --Regions_fasta "${meta.ID}_regions.fasta" --Peptides_out "${meta.ID}_peptides.fasta"
  
  seqkit sort -n -i "${meta.ID}_peptides.fasta" > "${meta.ID}_peptides_class_sorted.fasta"

  cut -f1,2,3,4,6,7 "${meta.ID}_peptide_regions.bed" > "${meta.ID}_peptide_regions_short.bed"
  gffread -w "${meta.ID}_peptide_regions_short.fasta" -g $genome "${meta.ID}_peptide_regions_short.bed"
  seqkit sort -n -i "${meta.ID}_peptide_regions_short.fasta" > "${meta.ID}_peptide_regions_short_sorted.fasta"

  split -l\$((`wc -l < "${meta.ID}_peptide_regions.bed"`/$split)) "${meta.ID}_peptide_regions.bed" "${meta.ID}_peptide_regions.split.bed" -da 4
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
  /bedtools2/bin/bedtools intersect -wo -s -a $GTF -b $inverted_bed > "${meta.ID}_stringtie_specific_new.bed"
  /scripts/Clean_overlapers.py --BED "${meta.ID}_stringtie_specific_new.bed" --GTF $GTF --TPM_min_novel $tpm_min_novel --Cov_min_novel $cov_min_novel --Out_anno "${meta.ID}_annotation.bed" --Out_bed "${meta.ID}_overlap.bed" --TPM_min_diff $tpm_min_diff --Cov_min_diff $cov_min_diff
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
  java -jar /jvarkit/dist/sortsamrefname.jar --samoutputformat BAM "${meta.ID}_survivor_coverage.bam" > "${meta.ID}_survivor_coverage_refname.bam"
  java -jar /jvarkit/dist/biostar154220.jar -n 100 --samoutputformat BAM "${meta.ID}_survivor_coverage_refname.bam" > "${meta.ID}_survivor_coverage_limited.bam"
  samtools sort "${meta.ID}_survivor_coverage_limited.bam" > "${meta.ID}_survivor_coverage_limited_sorted.bam" 
  samtools index "${meta.ID}_survivor_coverage_limited_sorted.bam" 
  java -jar /jvarkit/dist/sam4weblogo.jar -F tabular -r "${meta.ID}_regions.bed" "${meta.ID}_survivor_coverage_limited_sorted.bam" | grep -v '>' | grep -v "\\--" > "${meta.ID}_${bed.extension}_coverage_reads.tsv"
  /scripts/process_subreads.py --tsv_file "${meta.ID}_${bed.extension}_coverage_reads.tsv" --tsv_out "${meta.ID}_${bed.extension}_selected_sequences.bed" --BAM_cov $BAM_cov
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
  /scripts/translate_with_snps.py --fasta_in $selected --peptides $peptides --regions $regions --tsv_out "${meta.ID}_patient_specific_ncnas.tsv" --tsv_out_meta "${meta.ID}_patient_specific_ncnas_meta.tsv"
  """
}

process 'pVACbind_class_I' {

  errorStrategy 'ignore'

  input: 
  tuple val(meta), path(peptides_I), path(optiTypeOutput), path(reads), path(HLA_types_I), path(HLA_types_II)
  path(iedb)

  output:
  tuple val(meta), path("${peptides_I.baseName}_peptides_I_binding.tsv")
    
  script:
    """
    if [ $HLA_types_I == "HLA_default.txt" ]; then
    # Use OptiType output as input
      Class_I=\$(tr -s '\t'  '\n'<$optiTypeOutput | sed -n 11,16p | awk '\$1="HLA-"\$1' | awk '{print \$1}' | paste -s -d, -)
    else
    # Use HLA types file as input
      Class_I=\$(< $HLA_types_I)
    fi

    mkdir split_fastas

    awk '/^>/{s=\$0; next} {print s ORS \$0 > ("split_fastas/" length(\$0) ".fasta")}' $peptides_I

    # Loop through each filtered fasta file and run pvacbind
    for filtered_file in split_fastas/*.fasta; do
      peptide_length=\$(basename "\$filtered_file" .fasta)
      mkdir "${meta.ID}_\$peptide_length"
      pvacbind run --iedb-install-directory /opt/iedb \$filtered_file "${meta.ID}_\$peptide_length" \$Class_I -e1 \$peptide_length "NetMHCpan" "${meta.ID}_\$peptide_length"
    done

    awk 'FNR==1 && NR!=1 {next} {print}' */MHC_Class_I/*all_epitopes.tsv > "${peptides_I.baseName}_peptides_I_binding.tsv"
    """
}

process 'pVACbind_class_II' {

  errorStrategy 'ignore'

  input: 
  tuple val(meta), path (peptides_II), path(HLAHD_Output), path(reads), path(HLA_types_I), path(HLA_types_II)
  path(iedb)

  output:
  tuple val(meta), path("${peptides_II.baseName}_peptides_II_binding.tsv")
    
  script:
    """
    if [ $HLA_types_II == "HLA_II_default.txt" ]; then
      /scripts/HLA_HD_out_reshape.py --HLA_HD_out $HLAHD_Output --file_out HLA_HD_out_reshape.csv
      Class_II=\$(awk '{print \$1}' HLA_HD_out_reshape.csv |  paste -s -d, -)
      echo \$Class_II > Test.txt
    else
    # Use HLA types file as input
      Class_II=\$(< $HLA_types_II)
    fi

    mkdir split_fastas

    awk '/^>/{s=\$0; next} {print s ORS \$0 > ("split_fastas/" length(\$0) ".fasta")}' $peptides_II

    # Loop through each filtered fasta file and run pvacbind
    for filtered_file in split_fastas/*.fasta; do
      peptide_length=\$(basename "\$filtered_file" .fasta)
      mkdir "${meta.ID}_\$peptide_length"
      pvacbind run --iedb-install-directory /opt/iedb \$filtered_file "${meta.ID}_\$peptide_length" \$Class_II -e1 \$peptide_length "NetMHCIIpan" "${meta.ID}_\$peptide_length"
    done

    awk 'FNR==1 && NR!=1 {next} {print}' */MHC_Class_II/*all_epitopes.tsv > "${peptides_II.baseName}_peptides_II_binding.tsv"
    """
}

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
  errorStrategy 'ignore'

  input: 
  tuple val(meta), path(filtering), path(vaf), path(bed), path(translation), path(bed_2), path(specific), path(bind)
  path (anno_2)
    
  output:
  tuple val(meta), path ("${meta.ID}_final_out_combined.tsv")
    
  script:
    """
  /scripts/final_out_2.py --Filtering $filtering --VAF $vaf --BED $bed --Translation $translation --BED_2 $bed_2 --Specific $specific --BIND $bind --Out "${meta.ID}_final_out_combined_1.tsv" --Out_header "${meta.ID}_final_out_combined_0.tsv"
  /bedtools2/bin/bedtools intersect -wao -s -a "${meta.ID}_final_out_combined_1.tsv" -b $anno_2 > "${meta.ID}_final_out_combined_2.tsv"
  /scripts/final_out_1.py --BED "${meta.ID}_final_out_combined_2.tsv" --BED_old "${meta.ID}_final_out_combined_0.tsv" --Out "${meta.ID}_final_out_combined.tsv"
    """
}

process 'Metadata_MHCII' {
  errorStrategy 'ignore'

  input: 
  tuple val(meta), path(filtering), path(vaf), path(bed), path(translation), path(bed_2), path(specific), path(bind)
  path (anno_2)
    
  output:
  tuple val(meta), path ("${meta.ID}_final_out_combined.tsv")
    
  script:
    """
  /scripts/final_out_2.py --Filtering $filtering --VAF $vaf --BED $bed --Translation $translation --BED_2 $bed_2 --Specific $specific --BIND $bind --Out "${meta.ID}_final_out_combined_1.tsv" --Out_header "${meta.ID}_final_out_combined_0.tsv"
  /bedtools2/bin/bedtools intersect -wao -s -a "${meta.ID}_final_out_combined_1.tsv" -b $anno_2 > "${meta.ID}_final_out_combined_2.tsv"
  /scripts/final_out_1.py --BED "${meta.ID}_final_out_combined_2.tsv" --BED_old "${meta.ID}_final_out_combined_0.tsv" --Out "${meta.ID}_final_out_combined.tsv"
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

process 'Rerun_samplesheet' {
  errorStrategy 'ignore'

  input: 
  tuple val(meta), path(gtf), path(vaf), path(bam), path(bai), path(opti), path(hlahd), path(hla_I), path(hla_II)
  val(outdir) 
    
  output:
  tuple val(meta), path ("${meta.ID}_rerun_samplesheet.csv")
    
  script:
  def batch = file(params.input_fastq).splitCsv(header:true)
  def Read1 = file(batch.Read1)
  def Read2 = file(batch.Read2)
  def HLA_I = file(batch.HLA_types)
  def HLA_II = file(batch.HLA_types_II)
    """
    echo "ID,Read1,Read2,GTF,VAF,BAM,BAI,OPTI,HLAHD,HLA_types,HLA_types_II" > "${meta.ID}_rerun_samplesheet.csv"
    echo "$Read1,$Read2,${outdir}StringTie/$gtf,${outdir}StringTie/$gtf,${outdir}alignment/$bam,${outdir}alignment/$bai,${outdir}OptiType/$opti,${outdir}HLA_HD/$hlahd,${outdir}/$hla_I,${outdir}/$hla_II"
    """
}