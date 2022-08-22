cwlVersion: v1.0
class: CommandLineTool
id: kfdrc-vep105
label: VEP v105
doc: |
  Simplified description of what this tool does:
    1. Install needed plugins
    2. Untar cache if it is provided
    3. Run VEP on input VCF
    4. BGZIP output VCF
    5. TABIX output VCF

  VEP Parameters:
    1. input_file: Path to input file
    2. output_file: Path for output VCF or STDOUT
    3. stats_file: Path for output stats file
    4. warning_file: Path for output warnings file
    5. vcf: Writes output in VCF format
    6. offline: No database connections will be made, and a cache file or GFF/GTF file is required for annotation
    7. fork: Number of threads to run on
    8. ccds: Adds the CCDS transcript identifier
    9. uniprot: Adds best match accessions for translated protein products from three UniProt-related databases (SWISSPROT, TREMBL and UniParc)
    10. symbol: Adds the gene symbol (e.g. HGNC)
    11. numbers: Adds affected exon and intron numbering to to output. Format is Number/Total
    12. canonical: Adds a flag indicating if the transcript is the canonical transcript for the gene
    13. protein: Add the Ensembl protein identifier to the output where appropriate
    14. assembly: Select the assembly version to use if more than one available. If using the cache, you must have the appropriate assembly cache file installed
    15. dir_cache: Cache directory to use
    16. cache: Enables use of the cache
    17. merged: Use the merged Ensembl and RefSeq cache
    18. hgvs: Add HGVS nomenclature based on Ensembl stable identifiers to the output
    19. fasta: Specify a FASTA file or a directory containing FASTA files to use to look up reference sequence
    20. check_existing: Checks for the existence of known variants that are co-located with your input
    21. af_1kg: Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output
    22. af_esp: Include allele frequency from NHLBI-ESP populations
    23. af_gnomad: Include allele frequency from Genome Aggregation Database (gnomAD) exome populations
    24. plugin: Use named plugin
    25. custom: Add custom annotation to the output

  An example run of this tool will use a command like this:
    /bin/bash -c
    set -eo pipefail
    perl /opt/vep/src/ensembl-vep/INSTALL.pl
      --NO_TEST
      --NO_UPDATE
      --AUTO p
      --PLUGINS LoF,ExAC,gnomADc,CADD,dbNSFP,dbscSNV &&
    tar -xzf /path/to/cache.ext &&
    /opt/vep/src/ensembl-vep/vep
      --input_file /path/to/input_vcf.ext
      --output_file STDOUT
      --stats_file output_basename-string-value_stats.tool_name-string-value.html
      --warning_file output_basename-string-value_warnings.tool_name-string-value.txt
      --vcf
      --offline
      --fork $(inputs.cores)
      --ccds
      --uniprot
      --symbol
      --numbers
      --canonical
      --protein
      --assembly GRCh38
      --dir_cache $PWD
      --cache
      --merged
      --check_existing
      --af_1kg
      --af_esp
      --af_gnomad
      --hgvs
      --fasta /path/to/reference.ext
      --plugin CADD,/path/to/cadd_snvs.ext,/path/to/cadd_indels.ext
      --plugin dbNSFP,/path/to/dbnsfp.ext,ALL
      --plugin dbscSNV,/path/to/dbscsnv.ext
      --custom /path/to/phylop.ext,PhyloP,bigwig |
    bgzip -c > output_basename-string-value.tool_name-string-value.vep.vcf.gz &&
    tabix output_basename-string-value.tool_name-string-value.vep.vcf.gz

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: ${return inputs.ram * 1000}
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'ensemblorg/ensembl-vep:release_105.0'
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 1
    shellQuote: true
    valueFrom: >-
      set -eo pipefail

      ${
        var plugins = ["LoF","gnomADc"];
        if (inputs.cadd_indels) {
          plugins.push("CADD")
        }
        if (inputs.dbnsfp) {
          plugins.push("dbNSFP")
        }
        if (inputs.dbscsnv) {
          plugins.push("dbscSNV")
        }
        return "perl /opt/vep/src/ensembl-vep/INSTALL.pl --NO_TEST --NO_UPDATE --AUTO p --PLUGINS "+plugins.join(',')+" &&"
      }
      ${if(inputs.cache) {return "tar -xzf "+inputs.cache.path} else {return "echo 'No cache'"}} &&
      /opt/vep/src/ensembl-vep/vep
      --input_file $(inputs.input_vcf.path)
      --output_file STDOUT
      ${
        if (inputs.run_stats){
          var arg = " --stats_file " + inputs.output_basename + "_stats." + inputs.tool_name + ".html ";
          return arg;
        }
        else{
          return " --no_stats ";
        }
      }
      --warning_file $(inputs.output_basename)_warnings.$(inputs.tool_name).txt
      --buffer_size $(inputs.buffer_size)
      --vcf
      --domains
      --failed 1
      --gene_phenotype
      --no_escape
      --polyphen b
      --pubmed
      --regulatory
      --shift_hgvs 1
      --sift b
      --total_length
      --tsl
      --variant_class
      --xref_refseq
      --offline
      --fork $(inputs.cores)
      --ccds
      --uniprot
      --symbol
      --numbers
      --canonical
      --pick_order canonical,tsl,biotype,rank,ccds,length
      --flag_pick_allele
      --protein
      --assembly GRCh38
      --allele_number
      --dont_skip
      --allow_non_variant
      ${if(inputs.reference) {return "--hgvs --hgvsg --fasta " + inputs.reference.path} else {return ""}}
      ${if(inputs.cache) {return "--cache --dir_cache ."} else {return ""}}
      ${if(inputs.species) {return "--species " + inputs.species} else {return ""} }
      ${if(inputs.merged) {return "--merged"} else {return ""} }
      ${if(inputs.run_cache_existing) {return "--check_existing"} else {return ""} }
      ${if(inputs.cadd_indels && inputs.cadd_snvs) {return "--plugin CADD,"+inputs.cadd_snvs.path+","+inputs.cadd_indels.path} else {return ""}}
      ${if(inputs.run_cache_af) {return "--af_1kg --af_esp --af_gnomad"} else {return ""}}
      ${if(inputs.dbnsfp) {return "--plugin dbNSFP," + inputs.dbnsfp.path + "," + inputs.dbnsfp_fields} else {return ""}}
      ${if(inputs.dbscsnv) {return "--plugin dbscSNV,"+inputs.dbscsnv.path} else {return ""}}
      ${if(inputs.intervar) {return "--custom "+inputs.intervar.path+",Intervar,vcf,exact,0,STATUS"} else {return ""}}
      | bgzip -c -@ inputs.cores > $(inputs.output_basename).$(inputs.tool_name).vep.vcf.gz &&
      tabix $(inputs.output_basename).$(inputs.tool_name).vep.vcf.gz

inputs:
  input_vcf: { type: File, secondaryFiles: [.tbi], doc: "VCF file (with associated index) to be annotated" }
  ram: {type: 'int?', default: 32, doc: "In GB, may need to increase this value depending on the size/complexity of input"}
  cores: {type: 'int?', default: 16, doc: "Number of cores to use. May need to increase for really large inputs"}
  species: {type: 'string?', default: "homo_sapiens", doc: "Refer to the cache dir structure to set this"}
  buffer_size: {type: 'int?', default: 5000, doc: "Increase or decrease to balance speed and memory usage"}
  reference: { type: 'File?',  secondaryFiles: [.fai], doc: "Fasta genome assembly with indexes" }
  cache: { type: 'File?', doc: "tar gzipped cache from ensembl/local converted cache" }
  merged: { type: 'boolean?', doc: "Set to true if merged cache used", default: true }
  run_cache_existing: { type: boolean, doc: "Run the check_existing flag for cache" }
  run_cache_af: { type: boolean, doc: "Run the allele frequency flags for cache" }
  run_stats: { type: boolean, doc: "Create stats file? Disable for speed", default: true }
  cadd_indels: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD indel annotations" }
  cadd_snvs: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD SNV annotations" }
  dbnsfp: { type: 'File?', secondaryFiles: [.tbi,^.readme.txt], doc: "VEP-formatted plugin file, index, and readme file containing dbNSFP annotations" }
  dbnsfp_fields: { type: 'string?', doc: "csv string with desired fields to annotate. Use ALL to grab all",
    default: 'hg19_chr,hg19_pos\(1-based\),SIFT4G_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,VEST4_score,VEST4_rankscore,MetaSVM_pred,MetaLR_pred,MetaRNN_pred,M-CAP_pred,REVEL_score,REVEL_rankscore,PrimateAI_pred,DEOGEN2_pred,BayesDel_noAF_pred,ClinPred_pred,LIST-S2_pred,Aloft_pred,fathmm-MKL_coding_pred,fathmm-XF_coding_pred,Eigen-phred_coding,Eigen-PC-phred_coding,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,TWINSUK_AC,TWINSUK_AF,ALSPAC_AC,ALSPAC_AF,UK10K_AC,UK10K_AF,gnomAD_exomes_controls_AC,gnomAD_exomes_controls_AN,gnomAD_exomes_controls_AF,gnomAD_exomes_controls_nhomalt,gnomAD_exomes_controls_POPMAX_AC,gnomAD_exomes_controls_POPMAX_AN,gnomAD_exomes_controls_POPMAX_AF,gnomAD_exomes_controls_POPMAX_nhomalt,gnomAD_genomes_flag,gnomAD_genomes_AC,gnomAD_genomes_AN,gnomAD_genomes_AF,gnomAD_genomes_nhomalt,gnomAD_genomes_POPMAX_AC,gnomAD_genomes_POPMAX_AN,gnomAD_genomes_POPMAX_AF,gnomAD_genomes_POPMAX_nhomalt,gnomAD_genomes_controls_and_biobanks_AC,gnomAD_genomes_controls_and_biobanks_AN,gnomAD_genomes_controls_and_biobanks_AF,gnomAD_genomes_controls_and_biobanks_nhomalt,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,Interpro_domain,GTEx_V8_gene,GTEx_V8_tissue'
    }
  dbscsnv: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing dbscSNV annotations" }
  intervar: { type: 'File?', doc: "Intervar vcf-formatted file. See docs for custom build instructions", secondaryFiles: [.tbi] }
  output_basename: { type: string, doc: "String that will be used in the output filenames" }
  tool_name: { type: string, doc: "Tool name to be used in output filenames" }

outputs:
  output_vcf: { type: File, outputBinding: { glob: '*.vcf.gz' }, secondaryFiles: ['.tbi'] }
  output_html: { type: 'File?', outputBinding: { glob: '*.html' }}
  warn_txt: { type: 'File?', outputBinding: { glob: '*.txt' }}
