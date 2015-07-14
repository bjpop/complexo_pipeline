'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from stages import Stages


def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='complexo')
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')
    # Find the path to the reference genome
    reference_file = state.config.get_option('ref_grch37')
    # Stages are dependent on the state
    stages = Stages(state)

    # The original reference file
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_reference,
        name='original_reference',
        output=reference_file)

    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_fastqs,
        name='original_fastqs',
        output=fastq_files)

    # Index the reference using BWA 
    pipeline.transform(
        task_func=stages.index_reference_bwa,
        name='index_reference_bwa',
        input=output_from('original_reference'),
        filter=suffix('.fasta'),
        output=['.fasta.amb', '.fasta.ann', '.fasta.pac', '.fasta.sa', '.fasta.bwt'])
    
    # Index the reference using samtools 
    pipeline.transform(
        task_func=stages.index_reference_samtools,
        name='index_reference_samtools',
        input=output_from('original_reference'),
        filter=suffix('.fasta'),
        output='.fasta.fai')

    # Align paired end reads in FASTQ to the reference producing a BAM file
    (pipeline.transform(
        task_func=stages.align_bwa,
        name='align_bwa',
        input=output_from('original_fastqs'),
        # Match the R1 (read 1) FASTQ file and grab the path and sample name. 
        # This will be the first input to the stage.
        # We assume the sample name may consist of only alphanumeric
        # characters.
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9]+)_R1.fastq.gz'),
        # Add two more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        #    2. The reference genome file
        add_inputs=add_inputs(['{path[0]}/{sample[0]}_R2.fastq.gz', reference_file]),
        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the sample name. This is needed within the stage for finding out
        # sample specific configuration options
        extras=['{sample[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='{path[0]}/{sample[0]}.bam')
        # Ensure the reference is indexed before we run this stage
        .follows('index_reference_bwa')
        .follows('index_reference_samtools'))

    # Sort the BAM file using Picard 
    pipeline.transform(
        task_func=stages.sort_bam_picard,
        name='sort_bam_picard',
        input=output_from('align_bwa'),
        filter=suffix('.bam'),
        output='.sort.bam')

    # Mark duplicates in the BAM file using Picard 
    pipeline.transform(
        task_func=stages.mark_duplicates_picard,
        name='mark_duplicates_picard',
        input=output_from('sort_bam_picard'),
        filter=suffix('.sort.bam'),
        output=['.sort.dedup.bam', '.metricsdup'])

    # Generate chromosome intervals using GATK 
    pipeline.transform(
        task_func=stages.chrom_intervals_gatk,
        name='chrom_intervals_gatk',
        input=output_from('mark_duplicates_picard'),
        filter=suffix('.sort.dedup.bam'),
        add_inputs=add_inputs([reference_file]),
        output='.chr.intervals')

    return pipeline
