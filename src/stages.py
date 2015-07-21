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

PICARD_JAR = '$PICARD_HOME/lib/picard.jar'
GATK_JAR = '$GATK_HOME/GenomeAnalysisTK.jar'

def java_command(jar_path, mem_in_gb, command_args):
    '''Build a string for running a java command'''
    # Bit of room between Java's max heap memory and what was requested.
    # Allows for other Java memory usage, such as stack.
    java_mem = mem_in_gb - 2
    return 'java -Xmx{mem}g -jar {jar_path} {command_args}'.format(
        jar_path=jar_path, mem=java_mem, command_args=command_args)


def run_java(state, stage, jar_path, mem, args):
    command = java_command(jar_path, mem, args)
    run_stage(state, stage, command)


class Stages(object):
    def __init__(self, state):
        self.state = state
        self.reference = self.get_options('ref_grch37')
        self.mills_grch37 = self.get_options('mills_grch37')
        self.one_k_g_grch37_indels = self.get_options('one_k_g_grch37_indels')
        self.interval_grch37 = self.get_options('interval_grch37')
        #self.CEU_mergeGvcf = self.get_options('CEU_mergeGvcf')
        #self.GBR_mergeGvcf = self.get_options('GBR_mergeGvcf')
        #self.FIN_mergeGvcf = self.get_options('FIN_mergeGvcf')

    def run_picard(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, PICARD_JAR, mem, args)

    def run_gatk(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, GATK_JAR, mem, args)


    def get_stage_options(self, stage, *options):
        return self.state.config.get_stage_options(stage, *options)

    def get_options(self, *options):
        return self.state.config.get_options(*options)


    # def original_reference(self, output):
    #     '''Original reference file'''
    #     pass


    def original_fastqs(self, output):
        '''Original fastq files'''
        pass


    #def index_reference_bwa(self, reference_in, index_file_out):
    #    '''Index the reference genome using BWA'''
    #    command = "bwa index -a bwtsw {ref}".format(ref=reference_in)
    #    run_stage(self.state, 'index_reference_bwa', command)
    #

    #def index_reference_samtools(self, reference_in, index_file_out):
    #    '''Index the reference genome using samtools'''
    #    command = 'samtools faidx {ref}'.format(ref=reference_in)
    #    run_stage(self.state, 'index_reference_samtools', command)


    def align_bwa(self, inputs, bam_out, sample_id):
        '''Align the paired end fastq files to the reference genome using bwa'''
        fastq_read1_in, fastq_read2_in = inputs
        cores = self.get_stage_options('align_bwa', 'cores')
        read_group = '"@RG\tID:{sample}\tSM:{sample}\tPL:Illumina"'.format(sample=sample_id)
        command = 'bwa mem -t {cores} -R {read_group} {reference} {fastq_read1} {fastq_read2} ' \
                  '| samtools view -b -h -o {bam} -' \
                  .format(cores=cores,
                      read_group=read_group,
                      fastq_read1=fastq_read1_in,
                      fastq_read2=fastq_read2_in,
                      reference=self.reference,
                      bam=bam_out)
        run_stage(self.state, 'align_bwa', command)


    def sort_bam_picard(self, bam_in, sorted_bam_out):
        '''Sort the BAM file using Picard'''
        picard_args = 'SortSam INPUT={bam_in} OUTPUT={sorted_bam_out} ' \
                      'VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate ' \
                      'MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=True'.format(
                          bam_in=bam_in, sorted_bam_out=sorted_bam_out)
        self.run_picard('sort_bam_picard', picard_args)


    def mark_duplicates_picard(self, bam_in, outputs):
        '''Mark duplicate reads using Picard'''
        dedup_bam_out, metrics_out = outputs
        picard_args = 'MarkDuplicates INPUT={bam_in} OUTPUT={dedup_bam_out} ' \
                      'METRICS_FILE={metrics_out} VALIDATION_STRINGENCY=LENIENT ' \
                      'MAX_RECORDS_IN_RAM=5000000 ASSUME_SORTED=True ' \
                      'CREATE_INDEX=True'.format(bam_in=bam_in, dedup_bam_out=dedup_bam_out,
                          metrics_out=metrics_out)
        self.run_picard('mark_duplicates_picard', picard_args)


    def chrom_intervals_gatk(self, inputs, intervals_out):
        '''Generate chromosome intervals using GATK'''
        bam_in, _metrics_dup = inputs
        cores = self.get_stage_options('chrom_intervals_gatk', 'cores')
        gatk_args = '-T RealignerTargetCreator -R {reference} -I {bam} ' \
                    '--num_threads {threads} --known {mills_grch37} ' \
                    '--known {one_k_g_grch37_indels} -L {interval_grch37} ' \
                    '-o {out}'.format(reference=self.reference, bam=bam_in,
                            threads=cores, mills_grch37=self.mills_grch37,
                            one_k_g_grch37_indels=self.one_k_g_grch37_indels,
                            interval_grch37=self.interval_grch37,
                            out=intervals_out)
        self.run_gatk('chrom_intervals_gatk', gatk_args)


    def local_realignment_gatk(self, inputs, bam_out):
        '''Local realign reads using GATK'''
        target_intervals_in, bam_in = inputs
        gatk_args = "-T IndelRealigner -R {reference} -I {bam} -L {interval_grch37} " \
                    "-targetIntervals {target_intervals} -known {mills_grch37} " \
                    "-known {one_k_g_grch37_indels} " \
                    "-o {out}".format(reference=self.reference, bam=bam_in,
                            mills_grch37=self.mills_grch37,
                            one_k_g_grch37_indels=self.one_k_g_grch37_indels,
                            interval_grch37=self.interval_grch37,
                            target_intervals=target_intervals_in,
                            out=bam_out)
        run_gatk('local_realignment_gatk', gatk_args)


    # XXX I'm not sure that --num_cpu_threads_per_data_thread has any benefit here
    def base_recalibration_gatk(self, bam_in, outputs):
        '''Base recalibration using GATK'''
        csv_out, log_out = outputs
        gatk_args = "-T BaseRecalibrator -R {reference} -I {bam} " \
                    "--num_cpu_threads_per_data_thread 4 --knownSites {dbsnp_grch37} " \
                    "--knownSites {mills_grch37} --knownSites {one_k_g_grch37_indels} " \
                    "-log {log} -o {out}".format(reference=self.reference, bam=bam_in,
                            mills_grch37=self.mills_grch37, dbsnp_grch37=self.dbsnp_grch37,
                            one_k_g_grch37_indels=self.one_k_g_grch37_indels,
                            log=log_out, out=csv_out)
        run_gatk('base_recalibration_gatk', gatk_args)


    # XXX I'm not sure that --num_cpu_threads_per_data_thread has any benefit here
    def print_reads_gatk(self, inputs, bam_out):
        '''Print reads using GATK'''
        [csv_in, _log], bam_in = inputs
        gatk_args = "-T PrintReads -R {reference} -I {bam} --BQSR {recal_csv} " \
                    "-o {out} --num_cpu_threads_per_data_thread 4".format(reference=self.reference,
                            bam=bam_in, recal_csv=csv_in, out=bam_out)
        run_gatk('print_reads_gatk', gatk_args)


    def call_variants_gatk(self, bam_in, vcf_out):
        '''Call variants using GATK'''
        gatk_args = "-T HaplotypeCaller -R {reference} --min_base_quality_score 20 " \
                    "--variant_index_parameter 128000 --emitRefConfidence GVCF " \
                    "--standard_min_confidence_threshold_for_calling 30.0 " \
                    "--num_cpu_threads_per_data_thread 8 " \
                    "--variant_index_type LINEAR " \
                    "--standard_min_confidence_threshold_for_emitting 30.0 " \
                    "-I {bam} -L {interval_list} -o {out}".format(reference=self.reference,
                            bam=bam_in, interval_list=self.interval_grch37, out=vcf_out)
        run_gatk('call_variants_gatk', gatk_args)


    def combine_gvcf_gatk(self, vcf_files_in, vcf_out):
        '''Combine G.VCF files for all samples using GATK'''
        g_vcf_files = ' '.join(['--variant ' + vcf for vcf in vcf_files_in])
        gatk_args = "-T CombineGVCFs -R {reference} " \
                    "--disable_auto_index_creation_and_locking_when_reading_rods " \
                    "{g_vcf_files} -o {vcf_out}".format(reference=self.reference,
                        g_vcf_files=g_vcf_files, vcf_out=vcf_out)
        run_gatk('combine_gvcf_gatk', gatk_args)


    def genotype_gvcf_gatk(self, merged_vcf_in, vcf_out):
        '''Genotype G.VCF files using GATK'''
        cores = self.get_stage_option('genotype_gvcf_gatk', 'cores')
        gatk_args = "-T GenotypeGVCFs -R {reference} " \
                    "--disable_auto_index_creation_and_locking_when_reading_rods " \
                    "--num_threads {cores} --variant {merged_vcf} --out {vcf_out} " \
                    "--variant {CEU_mergeGvcf} --variant {GBR_mergeGvcf} " \
                    "--variant {FIN_mergeGvcf}".format(reference=self.reference,
                            cores=cores, merged_vcf=merged_vcf_in, vcf_out=vcf_out,
                            CEU_mergeGvcf=self.CEU_mergeGvcf, GBR_mergeGvcf=self.GBR_mergeGvcf,
                            FIN_mergeGvcf=self.FIN_mergeGvcf)
        run_gatk('genotype_gvcf_gatk', gatk_args)
