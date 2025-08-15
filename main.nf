// DEFAULT PARAMS

// mandatory
params.input = false // --input "data/*_R{1,2}.fastq.gz"
params.tool = ""

// optional
params.chunksize = false
params.nchunks = false

// PROCESSES, partially based on nf-core modules (https://nf-co.re/modules/)

process SEQKIT_SPLIT2 {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0' :
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("**/*.gz"), emit: reads
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = "${meta.args}"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqkit \\
        split2 \\
        $args \\
        --threads $task.cpus \\
        $reads \\
        --out-dir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

process BAKTA_BAKTA {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.10.4--pyhdfd78af_0' :
        'biocontainers/bakta:1.10.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("${prefix}.embl")             , emit: embl
    tuple val(meta), path("${prefix}.faa")              , emit: faa
    tuple val(meta), path("${prefix}.ffn")              , emit: ffn
    tuple val(meta), path("${prefix}.fna")              , emit: fna
    tuple val(meta), path("${prefix}.gbff")             , emit: gbff
    tuple val(meta), path("${prefix}.gff3")             , emit: gff
    tuple val(meta), path("${prefix}.hypotheticals.tsv"), emit: hypotheticals_tsv
    tuple val(meta), path("${prefix}.hypotheticals.faa"), emit: hypotheticals_faa
    tuple val(meta), path("${prefix}.tsv")              , emit: tsv
    tuple val(meta), path("${prefix}.txt")              , emit: txt
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ""
    prefix   = task.ext.prefix ?: "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_tf_cmd = prodigal_tf ? "--prodigal-tf ${prodigal_tf[0]}" : ""
    """
    bakta \\
        $fasta \\
        $args \\
        --threads $task.cpus \\
        --prefix $prefix \\
        $proteins_opt \\
        $prodigal_tf_cmd \\
        --db $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """
}

process PROKKA {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3af46b047c8fe84112adeaecf300878217c629b97f111f923ecf327656ddd141/data' :
        'community.wave.seqera.io/library/prokka_openjdk:10546cadeef11472' }"

    input:
    tuple val(meta), path(fasta)
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("${prefix}/*.gff"), emit: gff
    tuple val(meta), path("${prefix}/*.gbk"), emit: gbk
    tuple val(meta), path("${prefix}/*.fna"), emit: fna
    tuple val(meta), path("${prefix}/*.faa"), emit: faa
    tuple val(meta), path("${prefix}/*.ffn"), emit: ffn
    tuple val(meta), path("${prefix}/*.sqn"), emit: sqn
    tuple val(meta), path("${prefix}/*.fsa"), emit: fsa
    tuple val(meta), path("${prefix}/*.tbl"), emit: tbl
    tuple val(meta), path("${prefix}/*.err"), emit: err
    tuple val(meta), path("${prefix}/*.log"), emit: log
    tuple val(meta), path("${prefix}/*.txt"), emit: txt
    tuple val(meta), path("${prefix}/*.tsv"), emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args   ?: ''
    prefix               = task.ext.prefix ?: "${meta.id}"
    def input            = fasta.toString() - ~/\.gz$/
    def decompress       = fasta.getExtension() == "gz" ? "gunzip -c ${fasta} > ${input}" : ""
    def cleanup          = fasta.getExtension() == "gz" ? "rm ${input}" : ""
    def proteins_opt     = proteins ? "--proteins ${proteins}" : ""
    def prodigal_tf_in   = prodigal_tf ? "--prodigaltf ${prodigal_tf}" : ""
    """
    ${decompress}

    prokka \\
        ${args} \\
        --cpus ${task.cpus} \\
        --prefix ${prefix} \\
        ${proteins_opt} \\
        ${prodigal_tf_in} \\
        ${input}

    ${cleanup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}

process GFFREAD {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), path(gff)
    val args

    output:
    tuple val(meta), path("*_cleaned.gtf")  , emit: gtf             , optional: true
    tuple val(meta), path("*.gff3") , emit: gffread_gff     , optional: true
    tuple val(meta), path("*.fasta"), emit: gffread_fasta   , optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix           ?: "${meta.id}"+".${params.tool}"
    """
    # first remove all ? because gffread complains
    TAB=\$'\\t'
    sed "s/\$TAB?\$TAB/\$TAB.\$TAB/g" $gff > ${prefix}.noquestionmark.gtf
    
    gffread \\
        ${prefix}.noquestionmark.gtf \\
        $args \\
        -o ${prefix}.gtf
    
    # finally make sure all transcripts have gene_id
    cat ${prefix}.gtf | sed '/gene_id/!s/transcript_id\\([^;]*\\)/&; \\0/' | sed '/transcript_id.*transcript_id/s/transcript_id/gene_id/' > ${prefix}_cleaned.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}

process CAT_GTF {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    path(gff)
    path(input)

    output:
    path("*.gtf")       , emit: gtf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = "${input.baseName}."+params.tool
    """
	cat $gff > ${prefix}.gtf
    """
}

process CAT_FASTA {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    path(fasta)
    path(input)
    val(file_ending)

    output:
    path("*.${file_ending}")       , emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = "${input.baseName}."+params.tool
    """
	cat $fasta > ${prefix}.${file_ending}
    """
}

process CAT_TSV {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    path(tsv)
    path(input)
    val(file_ending)

    output:
    path("*.${file_ending}")       , emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = "${input.baseName}."+params.tool
    // make sure to print headeronly once (BAKTA: starts with hash, PROKKA: starts with "locus_tag" )
    def pattern     = params.tool == "bakta" ? "^#" : "^locus_tag"
    """
    grep \'$pattern\' ${tsv[0]} > ${prefix}.${file_ending}
	grep -v \'$pattern\' $tsv >> ${prefix}.${file_ending}
    """
}

// PIPELINE

workflow {
    // Check params and set channels from params
    if( !params.input ){ error("Cannot find any input files. Example: --input \"data/*.fasta\"") }
    ch_input = Channel.fromPath("${params.input}", checkIfExists: true)
    val_chunksize = params.chunksize ?: 1000
    val_nchunks = params.nchunks ?: 2
    if(params.chunksize){log.warn "Splitting into chunks of $val_chunksize records!"}
    if(!params.chunksize){log.warn "Splitting into $val_nchunks chunks !"}
    if( params.tool !in ["bakta","prokka"] ) { exit 1, '\nERROR: --tool must be in either "bakta" or "prokka"\n' }
    if( params.tool == "bakta" ) {
        if(!params.baktadb){exit 1, "\nERROR: --baktadb is required when using bakta\n"}
        ch_baktadb = Channel.fromPath("${params.baktadb}", checkIfExists: true)
    } else { ch_baktadb = Channel.empty() }

    // Set up input channel with metadata
    ch_input
        .map { fasta ->
                def meta = [:]
                meta.id           = "input"
                meta.single_end   = true
                meta.args         = params.chunksize ? "--by-size $val_chunksize -f -e .gz" : "--by-part $val_nchunks -f -e .gz"
                [ meta, fasta ] }
        .set { ch_fasta }

    // Splitting
    SEQKIT_SPLIT2 ( ch_fasta )

    // add chunk id
    SEQKIT_SPLIT2.out.reads
        .transpose()
        .map{ meta_old, fasta ->
            def meta = [:]
            meta.id           = "${fasta.baseName}"
            meta.part         = "${fasta.baseName}".split("\\.").find { it =~ /^part_/ } //the can be changed in SEQKIT_SPLIT2 with --by-part-prefix & --by-size-prefix, see https://bioinf.shenwei.me/seqkit/usage/#split2
            [ meta, fasta ] }
        .set{ ch_chunks_for_annotation }

    // Annotation
	if ( params.tool == "bakta" ) {
		BAKTA_BAKTA ( ch_chunks_for_annotation, ch_baktadb.collect(), [], [] )
		ch_annot = BAKTA_BAKTA.out
	}
	if ( params.tool == "prokka" ) {
		PROKKA ( ch_chunks_for_annotation, [], [] )
		ch_annot = PROKKA.out
	}
	
	// convert large gff3 file to tiny gtf file
	// see https://nfcore.slack.com/archives/CE8SSJV3N/p1737572156295019?thread_ts=1737463873.672569&cid=CE8SSJV3N
	// gffread results_test/bakta/MEGAHIT-group-all.contigs.head100.part_002.fa.gff3 --keep-exon-attrs -F -T --force-exons | sed '/gene_id/!s/transcript_id\([^;]*\)/&; \0/;/transcript_id.*transcript_id/s/transcript_id/gene_id/' > genomic_gffread_force-exons_geneid.gtf
	GFFREAD ( ch_annot.gff, '--keep-exon-attrs -F -T --force-exons' )

	/* THIS IS DONE ON THE HEADNODE; BAD IDEA! DONT DO THAT!
	GFFREAD.out.gtf
        .map{ meta, gff -> gff }
	    .collectFile(name: "${params.input}.gtf", newLine: false, cache: false, keepHeader: false, skip:0, sort: 'deep') //cache: false when this is last step only!
        .set { ch_gtf_total }
	ch_gtf_total.subscribe{ file(it).copyTo("${params.outdir}/") }
	*/

    // Summarize annotation in gtf
	GFFREAD.out.gtf
        .map{ meta, gff -> gff }
        .collect( sort: { it.baseName } )
        .set { ch_gtfs }
	CAT_GTF ( ch_gtfs, ch_input.collect() )
	
    // Summarize protein sequences in faa
	ch_annot.faa
        .map{ meta, fasta -> fasta }
        .collect( sort: { it.baseName } )
        .set { ch_fastas }
	CAT_FASTA ( ch_fastas, ch_input.collect(), "faa" )
	
    // Summarize annotation table as tsv
	ch_annot.tsv
        .map{ meta, fasta -> fasta }
        .collect( sort: { it.baseName } )
        .set { ch_fastas }
	CAT_TSV ( ch_fastas, ch_input.collect(), "tsv" )

}
