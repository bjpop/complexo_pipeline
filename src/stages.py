'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such 
as config, options, DRMAA and the logger.
'''

from utils import safe_make_dir
from runner import run_stage
import os

def java_command(jar_path, heap_mem_in_gb, command_args):
    '''Build a string for running a java command'''
    return 'java -Xmx{mem}g -jar {jar_path} {command_args}'.format(
        jar_path=jar_path, mem=heap_mem_in_gb, command_args=command_args)

class Stages(object):
    def __init__(self, state):
        self.state = state


    def original_reference(self, output):
        '''Original reference file'''
        pass


    def original_fastqs(self, output):
        '''Original fastq files'''
        pass


    def index_reference_bwa(self, reference_in, index_file_out):
        '''Index the reference genome using BWA'''
        command = "bwa index -a bwtsw {ref}".format(ref=reference_in)
        run_stage(self.state, 'index_reference_bwa', command)
    

    def index_reference_samtools(self, reference_in, index_file_out):
        '''Index the reference genome using samtools'''
        command = 'samtools faidx {ref}'.format(ref=reference_in)
        run_stage(self.state, 'index_reference_samtools', command)


    def align_bwa(self, inputs, bam_out, sample):
        '''Align the paired end fastq files to the reference genome using bwa'''
        fastq_read1_in, [fastq_read2_in, reference_in] = inputs
        # Get the number of cores to request for the job, this translates into the
        # number of threads to give to bwa's -t option
        cores = self.state.config.get_stage_option('align_bwa', 'cores')
        read_group = '"@RG\tID:{sample}\tSM:{sample}\tPL:Illumina"'.format(sample=sample)
        # Run bwa and pipe the output through samtools view to generate a BAM file
        command = 'bwa mem -t {cores} -R {read_group} {reference} {fastq_read1} {fastq_read2} ' \
                  '| samtools view -b -h -o {bam} -' \
                  .format(cores=cores,
                      read_group=read_group,
                      fastq_read1=fastq_read1_in,
                      fastq_read2=fastq_read2_in,
                      reference=reference_in,
                      bam=bam_out)
        run_stage(self.state, 'align_bwa', command)


    def sort_bam_picard(self, bam_in, sorted_bam_out):
        '''Sort the BAM file using Picard'''
        # Maximum memory requested in GB for the task
        mem = int(self.state.config.get_stage_option('sort_bam_picard', 'mem'))
        picard_args = 'SortSam INPUT={bam_in} OUTPUT={sorted_bam_out} ' \
                      'VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate ' \
                      'MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=True'.format(
                          bam_in=bam_in, sorted_bam_out=sorted_bam_out)
        # Bit of room between Java's max heap memory and what was requested.
        # Allows for other Java memory usage, such as stack.
        java_heap_mem = mem - 2
        command = java_command('$PICARD_HOME/lib/picard.jar', java_heap_mem, picard_args)
        run_stage(self.state, 'sort_bam_picard', command)


    def mark_duplicates_picard(self, bam_in, outputs):
        '''Mark duplicate reads using Picard'''
        dedup_bam_out, metrics_out = outputs
        # Maximum memory requested in GB for the task
        mem = int(self.state.config.get_stage_option('mark_duplicates_picard', 'mem'))
        picard_args = 'MarkDuplicates INPUT={bam_in} OUTPUT={dedup_bam_out} ' \
                      'METRICS_FILE={metrics_out} VALIDATION_STRINGENCY=LENIENT ' \
                      'MAX_RECORDS_IN_RAM=5000000 ASSUME_SORTED=True ' \
                      'CREATE_INDEX=True'.format(bam_in=bam_in, dedup_bam_out=dedup_bam_out,
                          metrics_out=metrics_out)
        # Bit of room between Java's max heap memory and what was requested.
        # Allows for other Java memory usage, such as stack.
        java_heap_mem = mem - 2
        command = java_command('$PICARD_HOME/lib/picard.jar', java_heap_mem, picard_args)
        run_stage(self.state, 'mark_duplicates_picard', command)


    def chrom_intervals_gatk(self, inputs, intervals_out):
        '''Generate chromosome intervals using GATK'''
        [bam_in, _metrics_dup], [reference_in] = inputs
        mem = int(self.state.config.get_stage_option('chrom_intervals_gatk', 'mem'))
        cores = self.state.config.get_stage_option('chrom_intervals_gatk', 'cores')
        mills_grch37 = self.state.config.get_option('mills_grch37')
        one_k_g_grch37_indels = self.state.config.get_option('one_k_g_grch37_indels')
        interval_grch37 = self.state.config.get_option('interval_grch37')
        gatk_args = '-T RealignerTargetCreator -R {reference} -I {bam} ' \
                    '--num_threads {threads} --known {mills_grch37} ' \
                    '--known {one_k_g_grch37_indels} -L {interval_grch37} ' \
                    '-o {out}'.format(reference=reference_in, bam=bam_in,
                            threads=cores, mills_grch37=mills_grch37,
                            one_k_g_grch37_indels=one_k_g_grch37_indels,
                            interval_grch37=interval_grch37,
                            out=intervals_out)
        java_heap_mem = mem - 2
        command = java_command('$GATK_HOME/GenomeAnalysisTK.jar', java_heap_mem, gatk_args)
        run_stage(self.state, 'chrom_intervals_gatk', command)


    def local_realignment_gatk(self, inputs, bam_out):
        '''Local realign reads using GATK'''
        target_intervals_in, [bam_in, reference_in] = inputs
        mem = int(self.state.config.get_stage_option('local_realignment_gatk', 'mem'))
        mills_grch37 = self.state.config.get_option('mills_grch37')
        one_k_g_grch37_indels = self.state.config.get_option('one_k_g_grch37_indels')
        interval_grch37 = self.state.config.get_option('interval_grch37')
        gatk_args = "-T IndelRealigner -R {reference} -I {bam} -L {interval_grch37} " \
                    "-targetIntervals {target_intervals} -known {mills_grch37} " \
                    "-known {one_k_g_grch37_indels} " \
                    "-o {out}".format(reference=reference_in, bam=bam_in,
                            mills_grch37=mills_grch37,
                            one_k_g_grch37_indels=one_k_g_grch37_indels,
                            interval_grch37=interval_grch37,
                            target_intervals=target_intervals_in,
                            out=bam_out)
        java_heap_mem = mem - 2
        command = java_command('$GATK_HOME/GenomeAnalysisTK.jar', java_heap_mem, gatk_args)
        run_stage(self.state, 'local_realignment_gatk', command)


    # XXX I'm not sure that --num_cpu_threads_per_data_thread has any benefit here
    def base_recalibration_gatk(self, inputs, outputs):
        '''Base recalibration using GATK'''
        bam_in, [reference_in] = inputs
        csv_out, log_out = outputs
        mem = int(self.state.config.get_stage_option('base_recalibration_gatk', 'mem'))
        mills_grch37 = self.state.config.get_option('mills_grch37')
        one_k_g_grch37_indels = self.state.config.get_option('one_k_g_grch37_indels')
        dbsnp_grch37 = self.state.config.get_option('dbsnp_grch37')
        gatk_args = "-T BaseRecalibrator -R {reference} -I {bam} " \
                    "--num_cpu_threads_per_data_thread 4 --knownSites {dbsnp_grch37} " \
                    "--knownSites {mills_grch37} --knownSites {one_k_g_grch37_indels} " \
                    "-log {log} -o {out}".format(reference=reference_in, bam=bam_in,
                            mills_grch37=mills_grch37, dbsnp_grch37=dbsnp_grch37,
                            one_k_g_grch37_indels=one_k_g_grch37_indels,
                            log=log_out, out=csv_out)
        java_heap_mem = mem - 2
        command = java_command('$GATK_HOME/GenomeAnalysisTK.jar', java_heap_mem, gatk_args)
        run_stage(self.state, 'base_recalibration_gatk', command)


    # XXX I'm not sure that --num_cpu_threads_per_data_thread has any benefit here
    def print_reads_gatk(self, inputs, bam_out):
        '''Print reads using GATK'''
        [csv_in, _log], [bam_in, reference_in] = inputs
        mem = int(self.state.config.get_stage_option('print_reads_gatk', 'mem'))
        gatk_args = "-T PrintReads -R {reference} -I {bam} --BQSR {recal_csv} " \
                    "-o {out} --num_cpu_threads_per_data_thread 4".format(reference=reference_in,
                            bam=bam_in, recal_csv=csv_in, out=bam_out)
        java_heap_mem = mem - 2
        command = java_command('$GATK_HOME/GenomeAnalysisTK.jar', java_heap_mem, gatk_args)
        run_stage(self.state, 'print_reads_gatk', command)


    def call_variants_gatk(self, inputs, vcf_out):
        '''Call variants using GATK'''
        bam_in, [reference_in] = inputs
        mem = int(self.state.config.get_stage_option('call_variants_gatk', 'mem'))
        interval_grch37 = self.state.config.get_option('interval_grch37')
        gatk_args = "-T HaplotypeCaller -R {reference} --min_base_quality_score 20 " \
                    "--variant_index_parameter 128000 --emitRefConfidence GVCF " \
                    "--standard_min_confidence_threshold_for_calling 30.0 " \
                    "--num_cpu_threads_per_data_thread 8 " \
                    "--variant_index_type LINEAR " \
                    "--standard_min_confidence_threshold_for_emitting 30.0 " \
                    "-I {bam} -L {interval_list} -o {out}".format(reference=reference_in,
                            bam=bam_in, interval_list=interval_grch37, out=vcf_out)
        java_heap_mem = mem - 2
        command = java_command('$GATK_HOME/GenomeAnalysisTK.jar', java_heap_mem, gatk_args)
        run_stage(self.state, 'call_variants_gatk', command)
