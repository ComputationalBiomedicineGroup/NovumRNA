#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {install_IEDB; Indices; OptiType; alignment; HLA_extraction; HLA_HD; StringTie; Create_capture_bed; Protein_to_peptides; Filtering; Translation; Translation_2; Annotation; Annotation_2; pVACbind_class_I; pVACbind_class_II; Metadata_MHCI; Metadata_MHCII; Combine; final_out_1; final_out_2} from "/home/ausserh/projects/2021/CRCnoncanonical/NovumRNA/novumRNA_modules_3.nf"


workflow analysis {
  def batchCSV = file(params.input_fastq).splitCsv(header:true)

  def raw_data = []

  for ( row in batchCSV ) {
      def meta  = [:]
      def reads = []
      
      meta.ID = row.ID
      meta.libType = "SE"
      if (row.Read1) { reads.add(file(row.Read1, checkIfExists: true)) }
          if (row.Read2) {
              reads.add(file(row.Read2, checkIfExists: true))
              meta.libType = "PE"
          }
      def hlaTypesFile = row.HLA_types
      if (hlaTypesFile) {
        def hlaTypesFilePath = file(hlaTypesFile, checkIfExists: true)
        raw_data.add([meta, reads, hlaTypesFilePath])
    } else {
        raw_data.add([meta, reads, "/home/ausserh/projects/2021/CRCnoncanonical/nonCanonicalNeoAG/HLA_default.txt"])
    }
      def hlaTypesFile_II = row.HLA_types_II
      if (hlaTypesFile_II) {
        def hlaTypesFilePath_II = file(hlaTypesFile_II, checkIfExists: true)
        raw_data.last() << hlaTypesFilePath_II
    } else {
        raw_data.last() << "/home/ausserh/projects/2021/CRCnoncanonical/nonCanonicalNeoAG/HLA_II_default.txt"
    }
  }
  def file = new File("/home/ausserh/raw_data.txt")
  batch_raw_data_ch = Channel.fromList(raw_data)
  install_IEDB(params.IEDB_MHCI_url, params.IEDB_MHCII_url, params.IEDB_check)
  Indices(params.genome, params.hisat_index, params.star_index)
  Protein_to_peptides(params.ref_proteome, params.Ref_pep, params.peptide_length)
  OptiType(batch_raw_data_ch, params.HLA_ref)
  alignment(params.genome, params.aligner, Indices.out[2], Indices.out[3], params.hisat_index, params.star_index, batch_raw_data_ch, params.two_pass, params.riboseq)  
  HLA_extraction(alignment.out[0].join(batch_raw_data_ch).transpose())
  HLA_HD(HLA_extraction.out.join(batch_raw_data_ch).transpose())
  StringTie(alignment.out[0], params.reference_GTF, params.longreads)
  Annotation(StringTie.out[0], params.tpm_min_novel, params.cov_min_novel, params.tpm_min_diff, params.cov_min_diff, params.inverted_bed)
  Translation(StringTie.out[0].join(Annotation.out[0]).transpose(), params.genome, params.reference_GTF, params.ref_proteome, params.peptide_length, params.split_anno_2)
  Annotation_2(Translation.out[0].join(alignment.out[0]).transpose(), params.aligner, params.BAM_cov)
  Combine(Annotation_2.out[0].join(Annotation_2.out[2]).groupTuple())
  Translation_2(Combine.out[0].join(Translation.out[2]).join(Translation.out[3]).transpose())
  Filtering(Translation_2.out[0], Protein_to_peptides.out, params.split_netMHCpan)
  pVACbind_class_I(Filtering.out[1].join(OptiType.out).join(batch_raw_data_ch).transpose(), install_IEDB.out)
  pVACbind_class_II(Filtering.out[3].join(HLA_HD.out).join(batch_raw_data_ch).transpose(), install_IEDB.out)
  final_out_1(pVACbind_class_I.out.groupTuple())
  final_out_2(pVACbind_class_II.out.groupTuple())
  Metadata_MHCI(Filtering.out[0].join(StringTie.out[2]).join(Annotation.out[1]).join(Translation.out[8]).join(Combine.out[1]).join(Translation_2.out[1]).join(final_out_1.out).groupTuple().transpose(), params.Annotation_2)
  Metadata_MHCII(Filtering.out[0].join(StringTie.out[2]).join(Annotation.out[1]).join(Translation.out[8]).join(Combine.out[1]).join(Translation_2.out[1]).join(final_out_2.out).groupTuple().transpose(), params.Annotation_2)
}

workflow analysis_short {
  def batchCSV_gtf = file(params.input_gtf).splitCsv(header:true)

  def raw_data_GTF = []
  def raw_data_BAM = []
  def raw_data_HLA = []

  for ( row in batchCSV_gtf ) {
      def meta  = [:]
      def GTF = []
      def BAM = []
      def BAI = []
      def HLA = []

      meta.ID = row.ID
      if (row.GTF) { GTF.add(file(row.GTF, checkIfExists: true)) }
      if (row.BAM) { BAM.add(file(row.BAM, checkIfExists: true)) }
      if (row.BAI) { BAI.add(file(row.BAI, checkIfExists: true)) }
      if (row.HLA) { HLA.add(file(row.HLA, checkIfExists: true)) }
      raw_data_GTF.add([meta, GTF])
      raw_data_BAM.add([meta, BAM, BAI])
      raw_data_HLA.add([meta, HLA])
  }

  batch_raw_data_ch_GTF = Channel.fromList(raw_data_GTF)
  batch_raw_data_ch_BAM = Channel.fromList(raw_data_BAM)
  batch_raw_data_ch_HLA = Channel.fromList(raw_data_HLA)

  HLA_extraction(batch_raw_data_ch_BAM.transpose())
  HLA_HD(HLA_extraction.out)
  Annotation(batch_raw_data_ch_GTF, params.tpm_min_novel, params.cov_min_novel, params.tpm_min_diff, params.cov_min_diff, params.inverted_bed)
  Translation(batch_raw_data_ch_GTF.join(Annotation.out[0]).transpose(), params.genome, params.reference_GTF, params.ref_proteome, params.peptide_length)
  Annotation_2(Translation.out[0].join(batch_raw_data_ch_BAM.transpose()).transpose(), params.quality, params.BAM_cov)
  Combine(Annotation_2.out[0].join(Annotation_2.out[2]).groupTuple())
  Translation_2(Combine.out[0].join(Translation.out[2]).join(Translation.out[3]).transpose())
  Filtering(Translation_2.out[0], params.database_I, params.split_netMHCpan)
  netMHCpan(Filtering.out[1].join(batch_raw_data_ch_HLA.transpose()).transpose())
  netMHCIIpan(Filtering.out[3].join(HLA_HD.out).transpose())
  final_out_1(netMHCpan.out.groupTuple())
  final_out_2(netMHCIIpan.out.groupTuple())
  Metadata_MHCI(Filtering.out[0].join(StringTie.out[2]).join(Annotation.out[1]).join(Translation.out[8]).join(Combine.out[1]).join(Translation_2.out[1]).join(final_out_1.out).groupTuple().transpose(), params.Annotation_2)
  Metadata_MHCII(Filtering.out[0].join(StringTie.out[2]).join(Annotation.out[1]).join(Translation.out[8]).join(Combine.out[1]).join(Translation_2.out[1]).join(final_out_2.out).groupTuple().transpose(), params.Annotation_2)
}

workflow capture_bed {

  def batchCSV = file(params.input_fastq).splitCsv(header:true)

  def raw_data = []

  for ( row in batchCSV ) {
      def meta  = [:]
      def reads = []
      
      meta.ID = row.ID
      meta.libType = "SE"
      if (row.Read1) { reads.add(file(row.Read1, checkIfExists: true)) }
          if (row.Read2) {
              reads.add(file(row.Read2, checkIfExists: true))
              meta.libType = "PE"
          }
      raw_data.add([meta, reads])
  }

  batch_raw_data_ch = Channel.fromList(raw_data)
  hisat(params.genome, params.hisat_index, batch_raw_data_ch, params.two_pass)
  StringTie(hisat.out[0], params.reference_GTF, params.longreads)
  Create_capture_bed(StringTie.out[1].collect(), params.tpm_max_diff, params.cov_max_diff)
}

workflow capture_bed_align {

  def batchCSV = file(params.input_fastq).splitCsv(header:true)

  def raw_data = []

  for ( row in batchCSV ) {
      def meta  = [:]
      def reads = []
      
      meta.ID = row.ID
      meta.libType = "SE"
      if (row.Read1) { reads.add(file(row.Read1, checkIfExists: true)) }
          if (row.Read2) {
              reads.add(file(row.Read2, checkIfExists: true))
              meta.libType = "PE"
          }
      raw_data.add([meta, reads])
  }

  batch_raw_data_ch = Channel.fromList(raw_data)
  OptiType(batch_raw_data_ch)
  hisat(params.genome, params.hisat_index, batch_raw_data_ch, params.two_pass)
  StringTie(hisat.out[0], params.reference_GTF, params.longreads)
}

workflow capture_bed_gtf {

  input_ch_capture_gtf = Channel.fromPath(params.input_capture_gtf)
              .splitCsv(header: false, sep:',')
              
  Create_capture_bed(input_ch_capture_gtf.collect(), params.tpm_max_diff, params.cov_max_diff, params.genome_length)
}