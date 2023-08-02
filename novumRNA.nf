#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {install_IEDB; Indices; OptiType; alignment; HLA_extraction; HLA_HD; StringTie; Create_capture_bed; Protein_to_peptides; Filtering; Translation; Translation_2; Annotation; Annotation_2; pVACbind_class_I; pVACbind_class_II; Metadata_MHCI; Metadata_MHCII; Combine; final_out_1; final_out_2; Rerun_samplesheet} from "./novumRNA_modules.nf"


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
        raw_data.add([meta, reads, file(params.default_hla)])
    }
      def hlaTypesFile_II = row.HLA_types_II
      if (hlaTypesFile_II) {
        def hlaTypesFilePath_II = file(hlaTypesFile_II, checkIfExists: true)
        raw_data.last() << hlaTypesFilePath_II
    } else {
        raw_data.last() << file(params.default_hla_II)
    }
  }
  batch_raw_data_ch = Channel.fromList(raw_data)
  install_IEDB(params.IEDB_MHCI_url, params.IEDB_MHCII_url, params.IEDB_check)
  Indices(params.genome, params.hisat_index, params.star_index)
  Protein_to_peptides(params.ref_proteome, params.Ref_pep, params.peptide_length)
  OptiType(batch_raw_data_ch, params.HLA_ref)
  alignment(params.genome, params.aligner, Indices.out[2], Indices.out[3], params.hisat_index, params.star_index, batch_raw_data_ch, params.two_pass, params.riboseq)  
  HLA_extraction(alignment.out[0].join(batch_raw_data_ch))
  HLA_HD(HLA_extraction.out.join(batch_raw_data_ch))
  StringTie(alignment.out[0], params.reference_GTF)
  Annotation(StringTie.out[0], params.tpm_min_novel, params.cov_min_novel, params.tpm_min_diff, params.cov_min_diff, params.capture_bed)
  Translation(StringTie.out[0].join(Annotation.out[0]).transpose(), params.genome, params.reference_GTF, params.ref_proteome, params.peptide_length, params.split_anno_2)
  Annotation_2(Translation.out[0].join(alignment.out[0]).transpose(), params.aligner, params.BAM_cov)
  Combine(Annotation_2.out[0].join(Annotation_2.out[2]).groupTuple())
  Translation_2(Combine.out[0].join(Translation.out[2]).join(Translation.out[3]).transpose())
  Filtering(Translation_2.out[0], Protein_to_peptides.out, params.split_netMHCpan)
  pVACbind_class_I(Filtering.out[1].join(OptiType.out).join(batch_raw_data_ch), install_IEDB.out)
  pVACbind_class_II(Filtering.out[3].join(HLA_HD.out).join(batch_raw_data_ch), install_IEDB.out)
  final_out_1(pVACbind_class_I.out.groupTuple())
  final_out_2(pVACbind_class_II.out.groupTuple())
  Metadata_MHCI(Filtering.out[0].join(StringTie.out[2]).join(Annotation.out[1]).join(Translation.out[8]).join(Combine.out[1]).join(Translation_2.out[1]).join(final_out_1.out).groupTuple().transpose(), params.Annotation_2)
  Metadata_MHCII(Filtering.out[0].join(StringTie.out[2]).join(Annotation.out[1]).join(Translation.out[8]).join(Combine.out[1]).join(Translation_2.out[1]).join(final_out_2.out).groupTuple().transpose(), params.Annotation_2)
  Rerun_samplesheet(batch_raw_data_ch.join(alignment.out[0]).join(StringTie.out[0]).join(StringTie.out[2]).join(OptiType.out).join(HLA_HD.out), params.outdir)
}

workflow analysis_short {
  def batchCSV_gtf = file(params.input_gtf).splitCsv(header:true)

  def raw_data = []
  def raw_data_GTF = []
  def raw_data_VAF = []
  def raw_data_BAM = []
  def raw_data_OPTI = []
  def raw_data_HLAHD = []
  def raw_data_HLA_types = []
  def raw_data_HLA_types_II = []

  for ( row in batchCSV_gtf ) {
      def meta  = [:]
      def reads = []
      def GTF = []
      def VAF = []
      def BAM = []
      def BAI = []
      def OPTI = []
      def HLAHD = []
      def HLA_types = []
      def HLA_types_II = []

      meta.ID = row.ID
      meta.libType = "SE"

      if (row.Read1) { reads.add(file(row.Read1, checkIfExists: true)) }
          if (row.Read2) {
              reads.add(file(row.Read2, checkIfExists: true))
              meta.libType = "PE"
          }

      if (row.GTF) { GTF.add(file(row.GTF, checkIfExists: true)) }
      if (row.VAF) { VAF.add(file(row.VAF, checkIfExists: true)) }
      if (row.BAM) { BAM.add(file(row.BAM, checkIfExists: true)) }
      if (row.BAI) { BAI.add(file(row.BAI, checkIfExists: true)) }
      if (row.OPTI) { OPTI.add(file(row.OPTI, checkIfExists: true)) }
      if (row.HLAHD) { HLAHD.add(file(row.HLAHD, checkIfExists: true)) }
      if (row.HLA_types) { HLA_types.add(file(row.HLA_types, checkIfExists: true)) }
      if (row.HLA_types_II) { HLA_types_II.add(file(row.HLA_types_II, checkIfExists: true)) }
      raw_data.add([meta, reads])
      raw_data_GTF.add([meta, GTF])
      raw_data_VAF.add([meta, VAF])
      raw_data_BAM.add([meta, BAM, BAI])
      raw_data_OPTI.add([meta, OPTI])
      raw_data_HLAHD.add([meta, HLAHD])
      raw_data_HLA_types.add([meta, HLA_types])
      raw_data_HLA_types_II.add([meta, HLA_types_II])
  }

  batch_raw_data_ch = Channel.fromList(raw_data)
  batch_raw_data_ch_GTF = Channel.fromList(raw_data_GTF)
  batch_raw_data_ch_VAF = Channel.fromList(raw_data_VAF)
  batch_raw_data_ch_BAM = Channel.fromList(raw_data_BAM)
  batch_raw_data_ch_OPTI = Channel.fromList(raw_data_OPTI)
  batch_raw_data_ch_HLAHD = Channel.fromList(raw_data_HLAHD)
  batch_raw_data_ch_HLA_types = Channel.fromList(raw_data_HLA_types)
  batch_raw_data_ch_HLA_types_II = Channel.fromList(raw_data_HLA_types_II)

  install_IEDB(params.IEDB_MHCI_url, params.IEDB_MHCII_url, params.IEDB_check)
  Protein_to_peptides(params.ref_proteome, params.Ref_pep, params.peptide_length)
  Annotation(batch_raw_data_ch_GTF, params.tpm_min_novel, params.cov_min_novel, params.tpm_min_diff, params.cov_min_diff, params.capture_bed)
  Translation(batch_raw_data_ch_GTF.join(Annotation.out[0]).transpose(), params.genome, params.reference_GTF, params.ref_proteome, params.peptide_length, params.split_anno_2)
  Annotation_2(Translation.out[0].join(batch_raw_data_ch_BAM.transpose()).transpose(), params.aligner, params.BAM_cov)
  Combine(Annotation_2.out[0].join(Annotation_2.out[2]).groupTuple())
  Translation_2(Combine.out[0].join(Translation.out[2]).join(Translation.out[3]).transpose())
  Filtering(Translation_2.out[0], Protein_to_peptides.out, params.split_netMHCpan)
  pVACbind_class_I(Filtering.out[1].join(batch_raw_data_ch_OPTI).join(batch_raw_data_ch).join(batch_raw_data_ch_HLA_types).join(batch_raw_data_ch_HLA_types_II).transpose(), install_IEDB.out)
  pVACbind_class_II(Filtering.out[3].join(batch_raw_data_ch_HLAHD).join(batch_raw_data_ch).join(batch_raw_data_ch_HLA_types).join(batch_raw_data_ch_HLA_types_II).transpose(), install_IEDB.out)
  final_out_1(pVACbind_class_I.out.groupTuple())
  final_out_2(pVACbind_class_II.out.groupTuple())
  Metadata_MHCI(Filtering.out[0].join(batch_raw_data_ch_VAF).join(Annotation.out[1]).join(Translation.out[8]).join(Combine.out[1]).join(Translation_2.out[1]).join(final_out_1.out).groupTuple().transpose(), params.Annotation_2)
  Metadata_MHCII(Filtering.out[0].join(batch_raw_data_ch_VAF).join(Annotation.out[1]).join(Translation.out[8]).join(Combine.out[1]).join(Translation_2.out[1]).join(final_out_2.out).groupTuple().transpose(), params.Annotation_2)
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
  Indices(params.genome, params.hisat_index, params.star_index)
  alignment(params.genome, params.aligner, Indices.out[2], Indices.out[3], params.hisat_index, params.star_index, batch_raw_data_ch, params.two_pass, params.riboseq)  
  StringTie(alignment.out[0], params.reference_GTF)
  Create_capture_bed(StringTie.out[1].collect(), params.tpm_max_diff, params.cov_max_diff, params.genome_length)
}

workflow capture_bed_short {

  input_ch_capture_gtf = Channel.fromPath(params.input_capture_gtf)
              .splitCsv(header: false, sep:',')
              
  Create_capture_bed(input_ch_capture_gtf.collect(), params.tpm_max_diff, params.cov_max_diff, params.genome_length)
}