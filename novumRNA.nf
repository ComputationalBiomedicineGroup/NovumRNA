#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {install_IEDB; Indices; OptiType; alignment; HLA_extraction; HLA_HD; StringTie; Create_capture_bed; Protein_to_peptides; Filtering; Translation; Translation_2; Capture_regions; BAM_coverage; pVACbind_class_I; pVACbind_class_II; Final_MHCI; Final_MHCII; Combine; Collect_binding; Collect_binding_II; Rerun_samplesheet} from "./novumRNA_modules.nf"


workflow analysis {

// Check if License(s) were accepted
params.accept_license = false
params.license        = "${params.input_ref}/LICENSE"

if (params.accept_license) {
    acceptLicense()
} else {
    checkLicense()
}

def raw_data = []
if (params.input_fastq == "${params.input_ref}/samplesheet_CRC_fastq_sub.csv") {
  raw_data.add([[ID:"Test_CRC01",libType:"PE"],["${params.input_ref}/AK11_CRC01_R1_combined_clean_rmdup.fastq", "${params.input_ref}/AK11_CRC01_R2_combined_clean_rmdup.fastq"],params.default_hla,"${params.input_ref}/HLA_class_II_default_alleles.txt"])
}
else {
  def batchCSV = file(params.input_fastq).splitCsv(header:true)

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
}

batch_raw_data_ch = Channel.fromList(raw_data)
install_IEDB(params.IEDB_MHCI_url, params.IEDB_MHCII_url, params.IEDB_check)
Indices(params.genome, params.hisat_index, params.star_index)
Protein_to_peptides(params.ref_proteome, params.Ref_pep, params.peptide_length)
OptiType(batch_raw_data_ch, params.HLA_ref)
alignment(params.genome, params.aligner, Indices.out[0], Indices.out[1], params.hisat_index, params.star_index, batch_raw_data_ch, params.two_pass, params.riboseq)  
HLA_extraction(alignment.out[0].join(batch_raw_data_ch))
HLA_HD(HLA_extraction.out.join(batch_raw_data_ch))
StringTie(alignment.out[0], params.reference_GTF)
Capture_regions(StringTie.out[0], params.tpm_min_novel, params.cov_min_novel, params.tpm_min_diff, params.cov_min_diff, params.capture_bed)
Translation(StringTie.out[0].join(Capture_regions.out[0]).transpose(), params.genome, params.reference_GTF, params.ref_proteome, params.peptide_length, params.split_BAM_coverage, params.ref_range)
BAM_coverage(Translation.out[0].join(alignment.out[0]).transpose(), params.aligner, params.BAM_cov)
Combine(BAM_coverage.out[0].join(BAM_coverage.out[2]).groupTuple())
Translation_2(Combine.out[0].join(Translation.out[2]).join(Translation.out[3]).transpose())
Filtering(Translation_2.out[0], Protein_to_peptides.out, params.split_pVACbind)
pVACbind_class_I(Filtering.out[1].join(OptiType.out).join(batch_raw_data_ch).transpose(by: 1), install_IEDB.out)
pVACbind_class_II(Filtering.out[3].join(HLA_HD.out).join(batch_raw_data_ch).transpose(by: 1), install_IEDB.out)
Collect_binding(pVACbind_class_I.out.groupTuple())
Collect_binding_II(pVACbind_class_II.out.groupTuple())
Final_MHCI(Filtering.out[0].join(StringTie.out[2]).join(Capture_regions.out[1]).join(Translation.out[8]).join(Combine.out[1]).join(Translation_2.out[1]).join(Collect_binding.out).groupTuple().transpose(), params.Annotation_2)
Final_MHCII(Filtering.out[0].join(StringTie.out[2]).join(Capture_regions.out[1]).join(Translation.out[8]).join(Combine.out[1]).join(Translation_2.out[1]).join(Collect_binding_II.out).groupTuple().transpose(), params.Annotation_2)
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
  Capture_regions(batch_raw_data_ch_GTF, params.tpm_min_novel, params.cov_min_novel, params.tpm_min_diff, params.cov_min_diff, params.capture_bed)
  Translation(batch_raw_data_ch_GTF.join(Capture_regions.out[0]).transpose(), params.genome, params.reference_GTF, params.ref_proteome, params.peptide_length, params.split_BAM_coverage, params.ref_range)
  BAM_coverage(Translation.out[0].join(batch_raw_data_ch_BAM).transpose(by: 1), params.aligner, params.BAM_cov)
  Combine(BAM_coverage.out[0].join(BAM_coverage.out[2]).groupTuple())
  Translation_2(Combine.out[0].join(Translation.out[2]).join(Translation.out[3]).transpose())
  Filtering(Translation_2.out[0], Protein_to_peptides.out, params.split_pVACbind)
  pVACbind_class_I(Filtering.out[1].join(batch_raw_data_ch_OPTI).join(batch_raw_data_ch).join(batch_raw_data_ch_HLA_types).join(batch_raw_data_ch_HLA_types_II).transpose(by: 1), install_IEDB.out)
  pVACbind_class_II(Filtering.out[3].join(batch_raw_data_ch_HLAHD).join(batch_raw_data_ch).join(batch_raw_data_ch_HLA_types).join(batch_raw_data_ch_HLA_types_II).transpose(by: 1), install_IEDB.out)
  Collect_binding(pVACbind_class_I.out.groupTuple())
  Collect_binding_II(pVACbind_class_II.out.groupTuple())
  Final_MHCI(Filtering.out[0].join(batch_raw_data_ch_VAF).join(Capture_regions.out[1]).join(Translation.out[8]).join(Combine.out[1]).join(Translation_2.out[1]).join(Collect_binding.out).groupTuple().transpose(), params.Annotation_2)
  Final_MHCII(Filtering.out[0].join(batch_raw_data_ch_VAF).join(Capture_regions.out[1]).join(Translation.out[8]).join(Combine.out[1]).join(Translation_2.out[1]).join(Collect_binding_II.out).groupTuple().transpose(), params.Annotation_2)
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
              meta.libType = "PE" }
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
  Indices(params.genome, params.hisat_index, params.star_index)
  alignment(params.genome, params.aligner, Indices.out[0], Indices.out[1], params.hisat_index, params.star_index, batch_raw_data_ch, params.two_pass, params.riboseq)  
  StringTie(alignment.out[0], params.reference_GTF)
  Create_capture_bed(StringTie.out[1].collect(), params.tpm_max_diff, params.cov_max_diff, params.genome_length)
}

workflow capture_bed_short {

  input_ch_capture_gtf = Channel.fromPath(params.input_capture_gtf)
              .splitCsv(header: false, sep:',')
              
  Create_capture_bed(input_ch_capture_gtf.collect(), params.tpm_max_diff, params.cov_max_diff, params.genome_length)
}

def showLicense() {

    licenseFile = file(params.license)
    log.info licenseFile.text

    log.info ""
    log.warn "To accept the license terms, please rerun with '--accept_license'"
    log.info ""

    exit 1
}

def acceptLicense() {
    log.info ""
    log.warn "I have read and accept the license terms"
    log.info ""

    licenseChckFile = file("${params.input_ref}/.license_accepted.chck")
    licenseChckFile.text = "License accepted"

    return true
}

def checkLicense() {
    licenseChckFile = file("${params.input_ref}/.license_accepted.chck")

    if(!licenseChckFile.exists()) {
        showLicense()
    } else {
        return true
    }
}
