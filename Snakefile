import snakemake.utils

snakemake.utils.min_version('6.5.0')

configfile: 'snakemake_config.yaml'

onsuccess:
    print('workflow success')

onerror:
    print('workflow error')

DEFAULT_MEM_MB=4 * 1024  # 4 GB
DEFAULT_TIME_HOURS=12

# Specifying this as an input to a rule will disable that rule.
# This can be used in combination with "ruleorder:" to determine what
# rule should be used to create a particular output file.
UNSATISFIABLE_INPUT='unsatisfiable_input_file_path'


def all_input(wildcards):
    inputs = dict()
    run_all_modules = bool(config.get('run_all_modules'))
    run_core_modules = bool(config.get('run_core_modules'))
    should_run_sjc = bool(config.get('should_run_sjc_steps'))
    if run_core_modules or run_all_modules:
        # core modules
        inputs.update(iris_epitope_post_out_files())
        inputs['visualization'] = os.path.join(result_dir(), 'visualization',
                                               'summary.png')

    if should_run_sjc:
        if has_tier_1():
            inputs['sjc_tier1'] = iris_append_sjc_out_file_name_for_tier('tier1')
        if has_tier_3():
            inputs['sjc_tier2tier3'] = iris_append_sjc_out_file_name_for_tier('tier2tier3')

    return inputs


localrules: all
rule all:
    input:
        unpack(all_input),


def result_dir():
    return os.path.join('results', config['run_name'])


def iris_db_path():
    return os.path.join(config['iris_data'], 'db')


def iris_db_sjc_path():
    return os.path.join(config['iris_data'], 'db_sjc')


def iris_exp_matrix_out_matrix():
    run_name = config['run_name']
    basename = 'exp.merged_matrix.{}.txt'.format(run_name)
    return os.path.join(result_dir(), 'exp_matrix', basename)


def gene_exp_matrix_path_for_run():
    from_config = config.get('gene_exp_matrix')
    if from_config:
        return from_config

    if config.get('run_all_modules'):
        return iris_exp_matrix_out_matrix()

    return None


def hla_types_list_for_run():
    from_config = config.get('mhc_list')
    if from_config:
        return from_config

    return os.path.join(result_dir(), 'hla_typing', 'hla_types.list')


def hla_from_patients_for_run():
    from_config = config.get('mhc_by_sample')
    if from_config:
        return from_config

    return os.path.join(result_dir(), 'hla_typing', 'hla_patient.tsv')

def splicing_matrix_path_for_run():
    db_path = iris_db_path()

    return os.path.join(db_path, config['run_name'], 'splicing_matrix')


def sjc_count_path_for_run():
    db_path = iris_db_sjc_path()

    return os.path.join(db_path, config['run_name'], 'sjc_matrix')


def splicing_matrix_txt_path_for_run():
    matrix_path = splicing_matrix_path_for_run()
    file_name = ('splicing_matrix.{}.cov10.{}.txt'
                 .format(config['splice_event_type'], config['run_name']))
    return os.path.join(matrix_path, file_name)


def splicing_matrix_idx_path_for_run():
    return '{}.idx'.format(splicing_matrix_txt_path_for_run())


def sjc_count_txt_path_for_run():
    matrix_path = sjc_count_path_for_run()
    file_name = 'SJ_count.{}.txt'.format(config['run_name'])
    return os.path.join(matrix_path, file_name)


def sjc_count_idx_path_for_run():
    return '{}.idx'.format(sjc_count_txt_path_for_run())


def format_ref_names(config_key):
    configured = config.get(config_key, '')
    # if no ref names -> provide a quoted empty string on the command line
    if not configured.strip():
        return "''"

    return configured


# must have either tier 1 or tier 3
def has_tier_1():
    return len(tier_1_group_names()) > 0


def has_tier_3():
    return len(tier_3_group_names()) > 0


def tier_1_group_names():
    return group_names_from_config_key('tissue_matched_normal_reference_group_names')


def tier_3_group_names():
    return group_names_from_config_key('normal_reference_group_names')


def group_names_from_config_key(key):
    names_str = config.get(key)
    split = names_str.split(',')
    return [x.strip() for x in split if x]


def reference_file_wildcard_constraints():
    reference_files = config.get('reference_files')
    if reference_files:
        file_names = '|'.join([re.escape(file_name)
                               for file_name in reference_files])
        without_gz = '|'.join([re.escape(file_name[:-3])
                               for file_name in reference_files
                               if file_name.endswith('.gz')])
    else:
        no_match = '^$'  # only matches empty string
        file_names = no_match
        without_gz = no_match

    return {'file_names': file_names, 'without_gz': without_gz}


def get_url_for_download_reference_file(wildcards):
    file_name = wildcards.file_name
    return config['reference_files'][file_name]['url']


rule download_reference_file:
    output:
        ref_file=os.path.join('references', '{file_name}'),
    log:
        out=os.path.join('references',
                         'download_reference_file_{file_name}_log.out'),
        err=os.path.join('references',
                         'download_reference_file_{file_name}_log.err'),
    wildcard_constraints:
        file_name=reference_file_wildcard_constraints()['file_names']
    params:
        url=get_url_for_download_reference_file,
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        'curl -L \'{params.url}\''
        ' -o {output.ref_file}'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule unzip_reference_file:
    input:
        gz=os.path.join('references', '{file_name}.gz'),
    output:
        un_gz=os.path.join('references', '{file_name}'),
    log:
        out=os.path.join('references',
                         'unzip_reference_file_{file_name}_log.out'),
        err=os.path.join('references',
                         'unzip_reference_file_{file_name}_log.err'),
    wildcard_constraints:
        file_name=reference_file_wildcard_constraints()['without_gz']
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        ' gunzip -c {input.gz}'
        ' 1> {output.un_gz}'
        ' 2> {log.err}'


def write_param_file_blocklist_param():
    value = config.get('blocklist')
    if value:
        return '--blocklist-file {}'.format(value)

    return ''


def write_param_file_bigwig_param():
    value = config.get('mapability_bigwig')
    if value:
        return '--mapability-bigwig {}'.format(value)

    return ''


def write_param_file_genome_param():
    value = config.get('fasta_name')
    if value:
        reference_path = os.path.join('references', value)
        return '--reference-genome {}'.format(reference_path)

    return ''


def write_param_file_input(wildcards):
    inputs = dict()
    fasta = config.get('fasta_name')
    if fasta:
        inputs['fasta'] = os.path.join('references', fasta)

    return inputs


rule write_param_file:
    input:
        unpack(write_param_file_input),
    output:
        param_file=os.path.join(result_dir(), 'screen.para'),
    log:
        out=os.path.join(result_dir(), 'write_param_file_log.out'),
        err=os.path.join(result_dir(), 'write_param_file_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_3=config['conda_env_3'],
        script=os.path.join('scripts', 'write_param_file.py'),
        group_name=config['run_name'],
        iris_db=iris_db_path(),
        matched_psi_cut=config.get('tissue_matched_normal_psi_p_value_cutoff', ''),
        matched_sjc_cut=config.get('tissue_matched_normal_sjc_p_value_cutoff', ''),
        matched_delta_psi_cut=config.get('tissue_matched_normal_delta_psi_p_value_cutoff', ''),
        matched_fc_cut=config.get('tissue_matched_normal_fold_change_cutoff', ''),
        matched_group_cut=config.get('tissue_matched_normal_group_count_cutoff', ''),
        matched_ref_names=format_ref_names('tissue_matched_normal_reference_group_names'),
        tumor_psi_cut=config.get('tumor_psi_p_value_cutoff', ''),
        tumor_sjc_cut=config.get('tumor_sjc_p_value_cutoff', ''),
        tumor_delta_psi_cut=config.get('tumor_delta_psi_p_value_cutoff', ''),
        tumor_fc_cut=config.get('tumor_fold_change_cutoff', ''),
        tumor_group_cut=config.get('tumor_group_count_cutoff', ''),
        tumor_ref_names=format_ref_names('tumor_reference_group_names'),
        normal_psi_cut=config.get('normal_psi_p_value_cutoff', ''),
        normal_sjc_cut=config.get('normal_sjc_p_value_cutoff', ''),
        normal_delta_psi_cut=config.get('normal_delta_psi_p_value_cutoff', ''),
        normal_fc_cut=config.get('normal_fold_change_cutoff', ''),
        normal_group_cut=config.get('normal_group_count_cutoff', ''),
        normal_ref_names=format_ref_names('normal_reference_group_names'),
        comparison_mode=config['comparison_mode'],
        stat_test_type=config['stat_test_type'],
        use_ratio='--use-ratio' if config.get('use_ratio') else '',
        blocklist=write_param_file_blocklist_param(),
        bigwig=write_param_file_bigwig_param(),
        genome=write_param_file_genome_param(),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} {params.conda_env_3} python {params.script}'
        ' --out-path {output.param_file}'
        ' --group-name {params.group_name}'
        ' --iris-db {params.iris_db}'
        ' --psi-p-value-cutoffs'
        ' {params.matched_psi_cut},{params.tumor_psi_cut},{params.normal_psi_cut}'
        ' --sjc-p-value-cutoffs'
        ' {params.matched_sjc_cut},{params.tumor_sjc_cut},{params.normal_sjc_cut}'
        ' --delta-psi-cutoffs'
        ' {params.matched_delta_psi_cut},{params.tumor_delta_psi_cut},{params.normal_delta_psi_cut}'
        ' --fold-change-cutoffs'
        ' {params.matched_fc_cut},{params.tumor_fc_cut},{params.normal_fc_cut}'
        ' --group-count-cutoffs'
        ' {params.matched_group_cut},{params.tumor_group_cut},{params.normal_group_cut}'
        ' --reference-names-tissue-matched-normal {params.matched_ref_names}'
        ' --reference-names-tumor {params.tumor_ref_names}'
        ' --reference-names-normal {params.normal_ref_names}'
        ' --comparison-mode {params.comparison_mode}'
        ' --statistical-test-type {params.stat_test_type}'
        ' {params.use_ratio}'
        ' {params.blocklist}'
        ' {params.bigwig}'
        ' {params.genome}'
        ' 1> {log.out}'
        ' 2> {log.err}'


# if the necessary files are specified in the config, then
# use them rather than run IRIS format
def copy_splice_matrix_files_input(wildcards):
    inputs = dict()
    inputs['splice_txt'] = config.get('splice_matrix_txt', UNSATISFIABLE_INPUT)
    inputs['splice_idx'] = config.get('splice_matrix_idx', UNSATISFIABLE_INPUT)
    if config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs

ruleorder: copy_splice_matrix_files > iris_format
localrules: copy_splice_matrix_files
rule copy_splice_matrix_files:
    input:
        unpack(copy_splice_matrix_files_input),
    output:
        splice_txt=splicing_matrix_txt_path_for_run(),
        splice_idx=splicing_matrix_idx_path_for_run(),
    shell:
        'cp {input.splice_txt} {output.splice_txt}'
        ' && cp {input.splice_idx} {output.splice_idx}'


def copy_sjc_count_files_input(wildcards):
    inputs = dict()
    inputs['count_txt'] = config.get('sjc_count_txt', UNSATISFIABLE_INPUT)
    inputs['count_idx'] = config.get('sjc_count_idx', UNSATISFIABLE_INPUT)
    if config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs

ruleorder: copy_sjc_count_files > iris_sjc_matrix
localrules: copy_sjc_count_files
rule copy_sjc_count_files:
    input:
        unpack(copy_sjc_count_files_input),
    output:
        count_txt=sjc_count_txt_path_for_run(),
        count_idx=sjc_count_idx_path_for_run(),
    shell:
        'cp {input.count_txt} {output.count_txt}'
        ' && cp {input.count_idx} {output.count_idx}'


def create_star_index_out_dir_param(wildcards, output):
    return os.path.dirname(output.index)


def create_star_index_input(wildcards):
    inputs = dict()
    inputs['gtf'] = os.path.join('references', config['gtf_name'])
    inputs['fasta'] = os.path.join('references', config['fasta_name'])
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule create_star_index:
    input:
        unpack(create_star_index_input),
    output:
        index=os.path.join('references', 'star_index', 'SA'),
    log:
        out=os.path.join('references', 'create_star_index_log.out'),
        err=os.path.join('references', 'create_star_index_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        out_dir=create_star_index_out_dir_param,
        overhang=config['star_sjdb_overhang'],
    threads: config['create_star_index_threads']
    resources:
        mem_mb=config['create_star_index_mem_gb'] * 1024,
        time_hours=config['create_star_index_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} STAR'
        ' --runMode genomeGenerate'
        ' --runThreadN {threads}'
        ' --genomeDir {params.out_dir}'
        ' --genomeFastaFiles {input.fasta}'
        ' --sjdbGTFfile {input.gtf}'
        ' --sjdbOverhang {params.overhang}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def organize_fastqs_sample_details():
    details = dict()
    fastq_dict = config.get('sample_fastqs')
    if not fastq_dict:
        return details

    sample_names = list()
    all_fastqs = list()
    for name, fastqs in fastq_dict.items():
        for fastq in fastqs:
            sample_names.append(name)
            all_fastqs.append(fastq)

    details['sample_names'] = sample_names
    details['fastqs'] = all_fastqs
    return details


def unique_sample_names():
    fastq_dict = config.get('sample_fastqs')
    if not fastq_dict:
        return list()

    return list(fastq_dict.keys())


def organize_fastqs_input(wildcards):
    details = organize_fastqs_sample_details()
    if not details:
        return {'unsatisfiable': UNSATISFIABLE_INPUT}

    return {'fastqs': details['fastqs']}


def organize_fastqs_sample_names_param():
    sample_names = organize_fastqs_sample_details().get('sample_names', list())
    return sample_names


localrules: organize_fastqs
rule organize_fastqs:
    input:
        unpack(organize_fastqs_input),
    output:
        done=touch(os.path.join(result_dir(), 'fastq_dir', 'organize_fastqs.done')),
    params:
        sample_names=organize_fastqs_sample_names_param(),
        out_dir=os.path.join(result_dir(), 'fastq_dir'),
    run:
        import os
        import os.path

        out_dir = params.out_dir
        if os.path.isdir(out_dir):
            files = os.listdir(out_dir)
            if files:
                raise Exception('organize_fastqs: {} already contains files'
                                .format(out_dir))

        for i, sample_name in enumerate(params.sample_names):
            sample_dir = os.path.join(out_dir, sample_name)
            orig_fastq_path = input.fastqs[i]
            fastq_basename = os.path.basename(orig_fastq_path)
            new_fastq_path = os.path.join(sample_dir, fastq_basename)
            os.makedirs(sample_dir, exist_ok=True)
            os.symlink(orig_fastq_path, new_fastq_path)


def iris_makesubsh_mapping_task_out_file_names():
    task_dir = os.path.join(result_dir(), 'mapping_tasks')
    sample_names = unique_sample_names()
    star_tasks = list()
    cuff_tasks = list()
    for sample_name in sample_names:
        star_name = 'STARmap.{}.sh'.format(sample_name)
        cuff_name = 'Cuffquant.{}.sh'.format(sample_name)
        star_tasks.append(os.path.join(task_dir, star_name))
        cuff_tasks.append(os.path.join(task_dir, cuff_name))

    return {'star_tasks': star_tasks, 'cuff_tasks': cuff_tasks}


def iris_makesubsh_mapping_star_done_file_names():
    out_dir = os.path.join(result_dir(), 'process_rnaseq')
    sample_names = unique_sample_names()
    final_bams = list()
    for sample in sample_names:
        align_dir = os.path.join(out_dir, '{}.aln'.format(sample))
        final_bam = os.path.join(align_dir, 'Aligned.sortedByCoord.out.bam')
        final_bams.append(final_bam)

    return final_bams


def iris_makesubsh_mapping_cuff_done_file_names():
    cuff_tasks = iris_makesubsh_mapping_task_out_file_names()['cuff_tasks']
    done_names = list()
    for task in cuff_tasks:
        done_names.append('{}.done'.format(task))

    return done_names


def iris_makesubsh_mapping_star_dir_param(wildcards):
    input = iris_makesubsh_mapping_input(wildcards)
    return os.path.dirname(input['index'])


def iris_makesubsh_mapping_task_dir_param(wildcards, output):
    return os.path.dirname(output.star_tasks[0])


def label_string_param():
    # IRIS uses this value to tell which files are for read 1 or read 2.
    # Specifically it looks for '1{label_string}f' and '1{label_string}f'
    return '.'


def iris_makesubsh_mapping_input(wildcards):
    inputs = dict()
    inputs['organize_fastqs_done'] = os.path.join(result_dir(), 'fastq_dir',
                                                  'organize_fastqs.done')
    inputs['index'] = os.path.join('references', 'star_index', 'SA')
    inputs['gtf'] = os.path.join('references', config['gtf_name'])
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_makesubsh_mapping:
    input:
        unpack(iris_makesubsh_mapping_input),
    output:
        **iris_makesubsh_mapping_task_out_file_names()
    log:
        out=os.path.join(result_dir(), 'iris_makesubsh_mapping_log.out'),
        err=os.path.join(result_dir(), 'iris_makesubsh_mapping_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        fastq_dir=os.path.join(result_dir(), 'fastq_dir'),
        star_dir=iris_makesubsh_mapping_star_dir_param,
        run_name=config['run_name'],
        out_dir=os.path.join(result_dir(), 'process_rnaseq'),
        label_string=label_string_param(),
        task_dir=iris_makesubsh_mapping_task_dir_param,
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS makesubsh_mapping'
        ' --fastq-folder-dir {params.fastq_dir}'
        ' --starGenomeDir {params.star_dir}'
        ' --gtf {input.gtf}'
        ' --data-name {params.run_name}'
        ' --outdir {params.out_dir}'
        ' --label-string {params.label_string}'
        ' --task-dir {params.task_dir}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_star_task_input(wildcards):
    inputs = dict()
    inputs['star_task'] = os.path.join(
        result_dir(), 'mapping_tasks',
        'STARmap.{}.sh'.format(wildcards.sample))
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_star_task:
    input:
        unpack(iris_star_task_input),
    output:
        unsorted_bam=os.path.join(result_dir(), 'process_rnaseq',
                                  '{sample}.aln', 'Aligned.out.bam'),
        sorted_bam=os.path.join(result_dir(), 'process_rnaseq', '{sample}.aln',
                                'Aligned.sortedByCoord.out.bam'),
    log:
        out=os.path.join(result_dir(), 'mapping_tasks',
                         'iris_star_task_{sample}_log.out'),
        err=os.path.join(result_dir(), 'mapping_tasks',
                         'iris_star_task_{sample}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
    threads: config['iris_star_task_threads']
    resources:
        mem_mb=config['iris_star_task_mem_gb'] * 1024,
        time_hours=config['iris_star_task_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} bash'
        ' {input.star_task}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_cuff_task_input(wildcards):
    inputs = dict()
    inputs['cuff_task'] = os.path.join(
        result_dir(), 'mapping_tasks',
        'Cuffquant.{}.sh'.format(wildcards.sample))
    inputs['star_task_done'] = os.path.join(
        result_dir(), 'process_rnaseq',
        '{}.aln'.format(wildcards.sample),
        'Aligned.sortedByCoord.out.bam')
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_cuff_task:
    input:
        unpack(iris_cuff_task_input),
    output:
        cuff_task_done=touch(os.path.join(result_dir(), 'process_rnaseq',
                                          '{sample}.aln', 'cufflinks',
                                          'genes.fpkm_tracking')),
    log:
        out=os.path.join(result_dir(), 'mapping_tasks',
                         'iris_cuff_task_{sample}_log.out'),
        err=os.path.join(result_dir(), 'mapping_tasks',
                         'iris_cuff_task_{sample}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
    threads: config['iris_cuff_task_threads']
    resources:
        mem_mb=config['iris_cuff_task_mem_gb'] * 1024,
        time_hours=config['iris_cuff_task_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} bash'
        ' {input.cuff_task}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_makesubsh_hla_task_out_file_names():
    task_dir = os.path.join(result_dir(), 'hla_tasks')
    sample_names = unique_sample_names()
    hla_tasks = list()
    for sample_name in sample_names:
        hla_name = 'seq2hla.{}.sh'.format(sample_name)
        hla_tasks.append(os.path.join(task_dir, hla_name))

    return {'hla_tasks': hla_tasks}


def iris_hla_task_done_file_names():
    sample_names = unique_sample_names()
    done_file_names = list()
    hla_dir = os.path.join(result_dir(), 'hla_typing')
    for sample in sample_names:
        out_dir = os.path.join(hla_dir, sample)
        expression = os.path.join(out_dir,
                                  '{}-ClassI.expression'.format(sample))
        genotype = os.path.join(out_dir,
                                '{}-ClassI.HLAgenotype4digits'.format(sample))
        done_file_names.append(expression)
        done_file_names.append(genotype)

    return done_file_names


def iris_makesubsh_hla_task_dir_param(wildcards, output):
    return os.path.dirname(output.hla_tasks[0])


def iris_makesubsh_hla_input(wildcards):
    inputs = dict()
    inputs['organize_fastqs_done'] = os.path.join(result_dir(), 'fastq_dir',
                                                  'organize_fastqs.done')
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_makesubsh_hla:
    input:
        unpack(iris_makesubsh_hla_input),

    output:
        **iris_makesubsh_hla_task_out_file_names()
    log:
        out=os.path.join(result_dir(), 'iris_makesubsh_hla_log.out'),
        err=os.path.join(result_dir(), 'iris_makesubsh_hla_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        fastq_dir=os.path.join(result_dir(), 'fastq_dir'),
        run_name=config['run_name'],
        out_dir=os.path.join(result_dir(), 'hla_typing'),
        label_string=label_string_param(),
        task_dir=iris_makesubsh_hla_task_dir_param,
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS makesubsh_hla'
        ' --fastq-folder-dir {params.fastq_dir}'
        ' --data-name {params.run_name}'
        ' --outdir {params.out_dir}'
        ' --label-string {params.label_string}'
        ' --task-dir {params.task_dir}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_hla_task_input(wildcards):
    inputs = dict()
    inputs['hla_task'] = os.path.join(
        result_dir(), 'hla_tasks',
        'seq2hla.{}.sh'.format(wildcards.sample))
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_hla_task:
    input:
        unpack(iris_hla_task_input),
    output:
        expression=os.path.join(result_dir(), 'hla_typing', '{sample}',
                                '{sample}-ClassI.expression'),
        genotype=os.path.join(result_dir(), 'hla_typing', '{sample}',
                              '{sample}-ClassI.HLAgenotype4digits'),
    log:
        out=os.path.join(result_dir(), 'hla_tasks',
                         'iris_hla_task_{sample}_log.out'),
        err=os.path.join(result_dir(), 'hla_tasks',
                         'iris_hla_task_{sample}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
    threads: config['iris_hla_task_threads']
    resources:
        mem_mb=config['iris_hla_task_mem_gb'] * 1024,
        time_hours=config['iris_hla_task_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} bash'
        ' {input.hla_task}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_parse_hla_input(wildcards):
    inputs = dict()
    inputs['hla_tasks_done'] = iris_hla_task_done_file_names()
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_parse_hla:
    input:
        unpack(iris_parse_hla_input),
    output:
        patient=os.path.join(result_dir(), 'hla_typing', 'hla_patient.tsv'),
        types=os.path.join(result_dir(), 'hla_typing', 'hla_types.list'),
        exp=os.path.join(result_dir(), 'hla_typing', 'hla_exp.list'),
    log:
        out=os.path.join(result_dir(), 'iris_parse_hla_log.out'),
        err=os.path.join(result_dir(), 'iris_parse_hla_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        out_dir=os.path.join(result_dir(), 'hla_typing'),
    resources:
        mem_mb=config['iris_parse_hla_mem_gb'] * 1024,
        time_hours=config['iris_parse_hla_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS parse_hla'
        ' --outdir {params.out_dir}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_makesubsh_rmats_task_out_file_names():
    task_dir = os.path.join(result_dir(), 'rmats_tasks')
    sample_names = unique_sample_names()
    rmats_tasks = list()
    for sample_name in sample_names:
        rmats_name = 'rMATS_prep.{}.sh'.format(sample_name)
        rmats_tasks.append(os.path.join(task_dir, rmats_name))

    return {'rmats_tasks': rmats_tasks}


def iris_makesubsh_rmats_done_file_names():
    rmats_tasks = iris_makesubsh_rmats_task_out_file_names()['rmats_tasks']
    done_names = list()
    for task in rmats_tasks:
        done_names.append('{}.done'.format(task))

    return done_names


def iris_makesubsh_rmats_task_dir_param(wildcards, output):
    return os.path.dirname(output.rmats_tasks[0])


def iris_makesubsh_rmats_input(wildcards):
    inputs = dict()
    inputs['star_done'] = iris_makesubsh_mapping_star_done_file_names()
    inputs['gtf'] = os.path.join('references', config['gtf_name'])
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_makesubsh_rmats:
    input:
        unpack(iris_makesubsh_rmats_input),
    output:
        **iris_makesubsh_rmats_task_out_file_names()
    log:
        out=os.path.join(result_dir(), 'iris_makesubsh_rmats_log.out'),
        err=os.path.join(result_dir(), 'iris_makesubsh_rmats_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        rmats_path=config['rmats_path'],
        bam_dir=os.path.join(result_dir(), 'process_rnaseq'),
        run_name=config['run_name'],
        task_dir=iris_makesubsh_rmats_task_dir_param,
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS makesubsh_rmats'
        ' --rMATS-path {params.rmats_path}'
        ' --bam-dir {params.bam_dir}'
        ' --gtf {input.gtf}'
        ' --data-name {params.run_name}'
        ' --task-dir {params.task_dir}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_rmats_task_input(wildcards):
    inputs = dict()
    inputs['rmats_task'] = os.path.join(
        result_dir(), 'rmats_tasks',
        'rMATS_prep.{}.sh'.format(wildcards.sample))
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_rmats_task:
    input:
        unpack(iris_rmats_task_input),
    output:
        # The output files have the format:
        # result_dir()/process_rnaseq/{run}.RL{readLength}/{sample}.tmp/{datetime}_{n}.rmats
        # Just using a .done file instead.
        rmats_task_done=touch(os.path.join(result_dir(), 'rmats_tasks',
                                           'rMATS_prep.{sample}.sh.done')),
    log:
        out=os.path.join(result_dir(), 'rmats_tasks',
                         'iris_rmats_task_{sample}_log.out'),
        err=os.path.join(result_dir(), 'rmats_tasks',
                         'iris_rmats_task_{sample}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
    threads: config['iris_rmats_task_threads']
    resources:
        mem_mb=config['iris_rmats_task_mem_gb'] * 1024,
        time_hours=config['iris_rmats_task_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} bash'
        ' {input.rmats_task}'
        ' 1> {log.out}'
        ' 2> {log.err}'

checkpoint check_read_lengths:
    input:
        rmats_done=iris_makesubsh_rmats_done_file_names(),
    output:
        read_lengths=os.path.join(result_dir(), 'process_rnaseq',
                                  'read_lengths.txt'),
    log:
        out=os.path.join(result_dir(), 'process_rnaseq', 'check_read_lengths_log.out'),
        err=os.path.join(result_dir(), 'process_rnaseq', 'check_read_lengths_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_3=config['conda_env_3'],
        script=os.path.join('scripts', 'check_read_lengths.py'),
        parent_dir=os.path.join(result_dir(), 'process_rnaseq'),
        run_name=config['run_name'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} {params.conda_env_3} python {params.script}'
        ' --parent-dir {params.parent_dir}'
        ' --run-name {params.run_name}'
        ' --out {output.read_lengths}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_makesubsh_rmatspost_input(wildcards):
    inputs = dict()
    inputs['rmats_done'] = iris_makesubsh_rmats_done_file_names()
    inputs['gtf'] = os.path.join('references', config['gtf_name'])
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_makesubsh_rmatspost:
    input:
        rmats_done=iris_makesubsh_rmats_done_file_names(),
        gtf=os.path.join('references', config['gtf_name']),
    output:
        makesubsh_done=touch(os.path.join(result_dir(), 'iris_makesubsh_rmats_post.done')),
    log:
        out=os.path.join(result_dir(), 'iris_makesubsh_rmatspost_log.out'),
        err=os.path.join(result_dir(), 'iris_makesubsh_rmatspost_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        rmats_path=config['rmats_path'],
        bam_dir=os.path.join(result_dir(), 'process_rnaseq'),
        run_name=config['run_name'],
        task_dir=os.path.join(result_dir(), 'rmats_post_tasks'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS makesubsh_rmatspost'
        ' --rMATS-path {params.rmats_path}'
        ' --bam-dir {params.bam_dir}'
        ' --gtf {input.gtf}'
        ' --data-name {params.run_name}'
        ' --task-dir {params.task_dir}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_rmatspost_task_input(wildcards):
    inputs = dict()
    inputs['makesubsh_done'] = os.path.join(
        result_dir(), 'iris_makesubsh_rmats_post.done')
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_rmatspost_task:
    input:
        unpack(iris_rmatspost_task_input),
    output:
        summary=os.path.join(result_dir(), 'process_rnaseq',
                             '{run_name}.RL{read_length}',
                             '{run_name}_RL{read_length}.matrix',
                             'summary.txt'),
    log:
        out=os.path.join(result_dir(), 'rmats_post_tasks',
                         'iris_rmatspost_task_{run_name}_{read_length}_log.out'),
        err=os.path.join(result_dir(), 'rmats_post_tasks',
                         'iris_rmatspost_task_{run_name}_{read_length}_log.err'),
    wildcard_constraints:
        run_name=config['run_name'],
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        # post_task is generated by iris_makesubsh_rmatspost, but since the number of
        # read_lengths is not known until the check_read_lengths checkpoint,
        # that .sh file is not used in the "input" or "output" sections of the snakemake
        post_task=os.path.join(result_dir(), 'rmats_post_tasks',
                               'rMATS_post.{run_name}_RL{read_length}.sh'),
    threads: config['iris_rmatspost_task_threads']
    resources:
        mem_mb=config['iris_rmatspost_task_mem_gb'] * 1024,
        time_hours=config['iris_rmatspost_task_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} bash'
        ' {params.post_task}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def prepare_iris_format_input(wildcards):
    read_lengths_file = checkpoints.check_read_lengths.get().output[0]
    summaries = list()
    run_name = config['run_name']
    input_prefix = os.path.join(result_dir(), 'process_rnaseq')
    with open(read_lengths_file, 'rt') as in_handle:
        for line in in_handle:
            read_length = line.strip()
            run_with_read_length_dot = '{}.RL{}'.format(run_name, read_length)
            run_with_read_length_underscore = '{}_RL{}'.format(run_name, read_length)
            summary = os.path.join(
                input_prefix, run_with_read_length_dot,
                '{}.matrix'.format(run_with_read_length_underscore),
                'summary.txt')
            summaries.append(summary)

    return {'summaries': summaries}


rule prepare_iris_format:
    input:
        unpack(prepare_iris_format_input),
    output:
        matrix=os.path.join(result_dir(), 'process_rnaseq', 'matrix_list.txt'),
        sample=os.path.join(result_dir(), 'process_rnaseq', 'sample_list.txt'),
    log:
        out=os.path.join(result_dir(), 'process_rnaseq',
                         'prepare_iris_format_log.out'),
        err=os.path.join(result_dir(), 'process_rnaseq',
                         'prepare_iris_format_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_3=config['conda_env_3'],
        script=os.path.join('scripts', 'prepare_iris_format.py'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} {params.conda_env_3} python {params.script}'
        ' --matrix-out {output.matrix}'
        ' --sample-out {output.sample}'
        ' --summaries {input.summaries}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_format_input(wildcards):
    inputs = dict()
    inputs['matrix_list'] = os.path.join(result_dir(), 'process_rnaseq',
                                         'matrix_list.txt')
    inputs['sample_list'] = os.path.join(result_dir(), 'process_rnaseq',
                                         'sample_list.txt')
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


# after IRIS format, IRIS index does not need to be run
rule iris_format:
    input:
        unpack(iris_format_input),
    output:
        matrix=splicing_matrix_txt_path_for_run(),
        idx=splicing_matrix_idx_path_for_run(),
    log:
        out=os.path.join(result_dir(), 'iris_format_log.out'),
        err=os.path.join(result_dir(), 'iris_format_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        splice_type=config['splice_event_type'],
        run_name=config['run_name'],
        iris_db=iris_db_path(),
        sample_name_field='2',
    resources:
        mem_mb=config['iris_format_mem_gb'] * 1024,
        time_hours=config['iris_format_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS format'
        ' {input.matrix_list}'
        ' {input.sample_list}'
        ' --splicing-event-type {params.splice_type}'
        ' --data-name {params.run_name}'
        ' --sample-name-field {params.sample_name_field}'
        ' --sample-based-filter'
        ' --iris-db-path {params.iris_db}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def prepare_iris_exp_matrix_input(wildcards):
    fpkm_files = list()
    sample_names = unique_sample_names()
    for name in sample_names:
        fpkm_path = os.path.join(
            result_dir(), 'process_rnaseq', '{}.aln'.format(name), 'cufflinks',
            'genes.fpkm_tracking')
        fpkm_files.append(fpkm_path)

    return {'fpkm_files': fpkm_files}


rule prepare_iris_exp_matrix:
    input:
        unpack(prepare_iris_exp_matrix_input),
    output:
        manifest=os.path.join(result_dir(), 'process_rnaseq',
                              'cufflinks_manifest.txt'),
    log:
        out=os.path.join(result_dir(), 'process_rnaseq',
                         'prepare_iris_exp_matrix_log.out'),
        err=os.path.join(result_dir(), 'process_rnaseq',
                         'prepare_iris_exp_matrix_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_3=config['conda_env_3'],
        script=os.path.join('scripts', 'prepare_iris_exp_matrix.py'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} {params.conda_env_3} python {params.script}'
        ' --out-manifest {output.manifest}'
        ' --fpkm-files {input.fpkm_files}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_exp_matrix_out_dir_param(wildcards, output):
    return os.path.dirname(output.matrix)


def iris_exp_matrix_input(wildcards):
    inputs = dict()
    inputs['manifest'] = os.path.join(result_dir(), 'process_rnaseq',
                                      'cufflinks_manifest.txt')
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_exp_matrix:
    input:
        unpack(iris_exp_matrix_input),
    output:
        matrix=iris_exp_matrix_out_matrix(),
    log:
        out=os.path.join(result_dir(), 'iris_exp_matrix_log.out'),
        err=os.path.join(result_dir(), 'iris_exp_matrix_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        run_name=config['run_name'],
        out_dir=iris_exp_matrix_out_dir_param,
    resources:
        mem_mb=config['iris_exp_matrix_mem_gb'] * 1024,
        time_hours=config['iris_exp_matrix_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS exp_matrix'
        ' --outdir {params.out_dir}'
        ' --data-name {params.run_name}'
        ' {input.manifest}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_makesubsh_extract_sjc_bam_dir_param(wildcards):
    input = iris_makesubsh_extract_sjc_input(wildcards)
    return os.path.dirname(input['bam'])


def iris_makesubsh_extract_sjc_task_dir_param(wildcards, output):
    return os.path.dirname(output.extract_task)


def iris_makesubsh_extract_sjc_input(wildcards):
    inputs = dict()
    inputs['bam'] = os.path.join(
        result_dir(), 'process_rnaseq', '{}.aln'.format(wildcards.sample),
        'Aligned.sortedByCoord.out.bam')
    inputs['gtf'] = os.path.join('references', config['gtf_name'])
    inputs['fasta'] = os.path.join('references', config['fasta_name'])
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_makesubsh_extract_sjc:
    input:
        unpack(iris_makesubsh_extract_sjc_input),
    output:
        extract_task=os.path.join(result_dir(), 'extract_sjc_tasks',
                                  'cmdlist.extract_sjc.{sample}'),
        bam_list=os.path.join(result_dir(), 'extract_sjc_tasks',
                              'bam_folder_list_{sample}.txt'),
    log:
        out=os.path.join(result_dir(), 'extract_sjc_tasks',
                         'iris_makesubsh_extract_sjc_{sample}_log.out'),
        err=os.path.join(result_dir(), 'extract_sjc_tasks',
                         'iris_makesubsh_extract_sjc_{sample}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        bam_dir=iris_makesubsh_extract_sjc_bam_dir_param,
        task_name='{sample}',
        bam_prefix='Aligned.sortedByCoord.out',
        task_dir=iris_makesubsh_extract_sjc_task_dir_param,
    resources:
    shell:
        'echo {params.bam_dir} > {output.bam_list}'
        ' && {params.conda_wrapper} {params.conda_env_2} IRIS'
        ' makesubsh_extract_sjc'
        ' --bam-folder-list {output.bam_list}'
        ' --task-name {params.task_name}'
        ' --gtf {input.gtf}'
        ' --genome-fasta {input.fasta}'
        ' --BAM-prefix {params.bam_prefix}'
        ' --task-dir {params.task_dir}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_extract_sjc_task_input(wildcards):
    inputs = dict()
    inputs['extract_task'] = os.path.join(
        result_dir(), 'extract_sjc_tasks',
        'cmdlist.extract_sjc.{}'.format(wildcards.sample))
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_extract_sjc_task:
    input:
        unpack(iris_extract_sjc_task_input),
    output:
        sj_count=os.path.join(result_dir(), 'process_rnaseq', '{sample}.aln',
                              'SJcount.txt'),
    log:
        out=os.path.join(result_dir(), 'extract_sjc_tasks',
                         'iris_extract_sjc_task_{sample}_log.out'),
        err=os.path.join(result_dir(), 'extract_sjc_tasks',
                         'iris_extract_sjc_task_{sample}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
    resources:
        mem_mb=config['iris_extract_sjc_task_mem_gb'] * 1024,
        time_hours=config['iris_extract_sjc_task_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} bash'
        ' {input.extract_task}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def prepare_iris_sjc_matrix_input(wildcards):
    sample_names = unique_sample_names()
    sj_files = list()
    for name in sample_names:
        sample_dir = '{}.aln'.format(name)
        sj_file = os.path.join(result_dir(), 'process_rnaseq', sample_dir,
                               'SJcount.txt')
        sj_files.append(sj_file)

    return {'sj_files': sj_files}


rule prepare_iris_sjc_matrix:
    input:
        unpack(prepare_iris_sjc_matrix_input),
    output:
        sj_list=os.path.join(result_dir(), 'process_rnaseq', 'sjc_file_list.txt'),
    log:
        out=os.path.join(result_dir(), 'process_rnaseq',
                         'prepare_iris_sjc_matrix_log.out'),
        err=os.path.join(result_dir(), 'process_rnaseq',
                         'prepare_iris_sjc_matrix_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_3=config['conda_env_3'],
        script=os.path.join('scripts', 'prepare_iris_sjc_matrix.py'),
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} {params.conda_env_3} python {params.script}'
        ' --sj-out {output.sj_list}'
        ' --sj-files {input.sj_files}'
        ' 1> {log.out}'
        ' 2> {log.err}'

def iris_sjc_matrix_input(wildcards):
    inputs = dict()
    inputs['sj_list'] = os.path.join(result_dir(), 'process_rnaseq',
                                     'sjc_file_list.txt')
    if not config['run_all_modules']:
        inputs['run_all_modules'] = UNSATISFIABLE_INPUT

    return inputs


rule iris_sjc_matrix:
    input:
        unpack(iris_sjc_matrix_input),
    output:
        count_txt=sjc_count_txt_path_for_run(),
        count_idx=sjc_count_idx_path_for_run(),
    log:
        out=os.path.join(result_dir(), 'iris_sjc_matrix_log.out'),
        err=os.path.join(result_dir(), 'iris_sjc_matrix_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        run_name=config['run_name'],
        sample_name_field='2',
        db_sjc=iris_db_sjc_path(),
    resources:
        mem_mb=config['iris_sjc_matrix_mem_gb'] * 1024,
        time_hours=config['iris_sjc_matrix_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS sjc_matrix'
        ' --file-list-input {input.sj_list}'
        ' --data-name {params.run_name}'
        ' --sample-name-field {params.sample_name_field}'
        ' --iris-db-path {params.db_sjc}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_screen_out_dir_param(wildcards, output):
    return os.path.dirname(output.guided)


def iris_screen_out_files():
    out_files = dict()
    out_dir = os.path.join(result_dir(), 'screen')
    out_prefix = '{}.{}'.format(config['run_name'], config['splice_event_type'])
    out_files['guided'] = os.path.join(
        out_dir, '{}.test.all_guided.txt'.format(out_prefix))
    out_files['voted'] = os.path.join(
        out_dir, '{}.test.all_voted.txt'.format(out_prefix))
    out_files['notest'] = os.path.join(
        out_dir, '{}.notest.txt'.format(out_prefix))
    out_files['tier1'] = os.path.join(
        out_dir, '{}.tier1.txt'.format(out_prefix))
    out_files['tier2tier3'] = os.path.join(
        out_dir, '{}.tier2tier3.txt'.format(out_prefix))

    return out_files


rule iris_screen:
    input:
        parameter_file=os.path.join(result_dir(), 'screen.para'),
        gtf=os.path.join('references', config['gtf_name']),
        splice_txt=splicing_matrix_txt_path_for_run(),
        splice_idx=splicing_matrix_idx_path_for_run(),
    output:
        **iris_screen_out_files()
    log:
        out=os.path.join(result_dir(), 'iris_screen_log.out'),
        err=os.path.join(result_dir(), 'iris_screen_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        splice_event_type=config['splice_event_type'],
        out_dir=iris_screen_out_dir_param,
    resources:
        mem_mb=config['iris_screen_mem_gb'] * 1024,
        time_hours=config['iris_screen_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS screen'
        ' --parameter-fin {input.parameter_file}'
        ' --splicing-event-type {params.splice_event_type}'
        ' --outdir {params.out_dir}'
        ' --translating'  # runs IRIS translate within IRIS screen
        ' --gtf {input.gtf}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_predict_out_dir_param(wildcards, output):
    return os.path.dirname(output.predict_out[0])


def iris_predict_out_file_names():
    tier_names = list()
    if has_tier_1():
        tier_names.append('tier1')
    if has_tier_3():
        tier_names.append('tier2tier3')

    names = list()
    for tier_name in tier_names:
        basename = '{}.{}.{}.txt.ExtraCellularAS.txt'.format(
            config['run_name'], config['splice_event_type'], tier_name)
        name = os.path.join(result_dir(), 'screen', basename)
        names.append(name)

    return names


def iris_predict_input(wildcards):
    inputs = dict()
    inputs['parameter_file'] = os.path.join(result_dir(), 'screen.para')
    inputs['screen_out'] = iris_screen_out_files()['guided']
    inputs['mhc_list'] = hla_types_list_for_run()
    gene_path = gene_exp_matrix_path_for_run()
    if gene_path:
        inputs['gene_exp_matrix'] = gene_path

    return inputs


def iris_predict_gene_exp_param():
    gene_exp_path = gene_exp_matrix_path_for_run()
    if not gene_exp_path:
        return ''

    return ' --gene-exp-matrix {}'.format(gene_exp_path)


rule iris_predict:
    input:
        unpack(iris_predict_input),
    output:
        predict_out=iris_predict_out_file_names(),
    log:
        out=os.path.join(result_dir(), 'iris_predict_log.out'),
        err=os.path.join(result_dir(), 'iris_predict_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        out_dir=iris_predict_out_dir_param,
        task_dir=os.path.join(result_dir(), 'predict_tasks'),
        splice_event_type=config['splice_event_type'],
        iedb_path=config['iedb_path'],
        task_wildcard_string=os.path.join(result_dir(), 'predict_tasks',
                                          'pep2epitope_{}.tier*.*.sh'.format(
                                              config['splice_event_type'])),
        gene_exp=iris_predict_gene_exp_param(),
    resources:
        mem_mb=config['iris_predict_mem_gb'] * 1024,
        time_hours=config['iris_predict_time_hr'],
    shell:
        # Remove any existing task scripts.
        # Usually snakemake automatically removes output before running a job, but
        # in this case the number of output files is not known in advance.
        'if [[ -n "$(ls {params.task_wildcard_string})" ]];'
        ' then rm {params.task_wildcard_string};'
        ' fi;'
        ' {params.conda_wrapper} {params.conda_env_2} IRIS predict'
        ' {params.out_dir}'
        ' --task-dir {params.task_dir}'
        ' --parameter-fin {input.parameter_file}'
        ' --splicing-event-type {params.splice_event_type}'
        ' --iedb-local {params.iedb_path}'
        ' --mhc-list {input.mhc_list}'
        ' {params.gene_exp}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def count_iris_predict_tasks_task_dir_param(wildcards, output):
    return os.path.dirname(output.predict_task_list)


checkpoint count_iris_predict_tasks:
    input:
        predict_out=iris_predict_out_file_names(),
    output:
        predict_task_list=os.path.join(result_dir(), 'predict_tasks',
                                       'predict_tasks_list.txt'),
    log:
        out=os.path.join(result_dir(), 'predict_tasks',
                         'count_iris_predict_tasks_log.out'),
        err=os.path.join(result_dir(), 'predict_tasks',
                         'count_iris_predict_tasks_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_3=config['conda_env_3'],
        script=os.path.join('scripts', 'count_iris_predict_tasks.py'),
        task_dir=count_iris_predict_tasks_task_dir_param,
        splice_type=config['splice_event_type'],
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        '{params.conda_wrapper} {params.conda_env_3} python {params.script}'
        '  --out-list {output.predict_task_list}'
        ' --task-dir {params.task_dir}'
        ' --splice-type {params.splice_type}'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule iris_predict_task:
    input:
        predict_task=os.path.join(result_dir(), 'predict_tasks', '{task_name}.sh'),
    output:
        predict_task_done=touch(os.path.join(result_dir(), 'predict_tasks',
                                             '{task_name}.sh.done')),
    log:
        out=os.path.join(result_dir(), 'predict_tasks', '{task_name}_log.out'),
        err=os.path.join(result_dir(), 'predict_tasks', '{task_name}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
    resources:
        mem_mb=config['iris_predict_task_mem_gb'] * 1024,
        time_hours=config['iris_predict_task_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} bash'
        ' {input.predict_task}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_epitope_post_out_dir_param():
    return os.path.join(result_dir(), 'screen')


def iris_epitope_post_out_files():
    files = dict()
    splice_type = config['splice_event_type']
    tier1_dir = os.path.join(result_dir(), 'screen', '{}.tier1'.format(splice_type))
    if has_tier_1():
        files['tier1_junction'] = os.path.join(tier1_dir, 'epitope_summary.junction-based.txt')
        files['tier1_peptide'] = os.path.join(tier1_dir, 'epitope_summary.peptide-based.txt')
        files['tier1_filtered'] = os.path.join(tier1_dir, 'pred_filtered.score500.txt')

    tier2tier3_dir = os.path.join(result_dir(), 'screen',
                                  '{}.tier2tier3'.format(splice_type))
    if has_tier_3():
        files['tier2tier3_junction'] = os.path.join(
            tier2tier3_dir, 'epitope_summary.junction-based.txt')
        files['tier2tier3_peptide'] = os.path.join(
            tier2tier3_dir, 'epitope_summary.peptide-based.txt')
        files['tier2tier3_filtered'] = os.path.join(
            tier2tier3_dir, 'pred_filtered.score500.txt')

    return files


def iris_predict_task_done_file_names():
    tasks_list_file = checkpoints.count_iris_predict_tasks.get().output[0]
    predict_tasks_done = list()
    with open(tasks_list_file, 'rt') as handle:
        for line in handle:
            task_file = line.strip()
            task_done_file = '{}.done'.format(task_file)
            predict_tasks_done.append(task_done_file)

    return predict_tasks_done


def iris_epitope_post_input(wildcards):
    inputs = dict()
    predict_tasks_done = iris_predict_task_done_file_names()
    inputs['predict_tasks_done'] = predict_tasks_done
    inputs['parameter_file'] = os.path.join(result_dir(), 'screen.para')
    inputs['mhc_by_sample'] = hla_from_patients_for_run()
    gene_exp = gene_exp_matrix_path_for_run()
    if gene_exp:
        inputs['gene_exp_matrix'] = gene_exp

    return inputs


def iris_epitope_post_gene_exp_param():
    gene_exp_path = gene_exp_matrix_path_for_run()
    if not gene_exp_path:
        return ''

    return ' --gene-exp-matrix {}'.format(gene_exp_path)


rule iris_epitope_post:
    input:
        unpack(iris_epitope_post_input),
    output:
        **iris_epitope_post_out_files()
    log:
        out=os.path.join(result_dir(), 'iris_epitope_post_log.out'),
        err=os.path.join(result_dir(), 'iris_epitope_post_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        out_dir=iris_epitope_post_out_dir_param(),
        splice_event_type=config['splice_event_type'],
        gene_exp=iris_epitope_post_gene_exp_param(),
    resources:
        mem_mb=config['iris_epitope_post_mem_gb'] * 1024,
        time_hours=config['iris_epitope_post_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS epitope_post'
        ' --parameter-fin {input.parameter_file}'
        ' --outdir {params.out_dir}'
        ' --splicing-event-type {params.splice_event_type}'
        ' --mhc-by-sample {input.mhc_by_sample}'
        ' {params.gene_exp}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_screen_sjc_out_dir_param(wildcards, output):
    return os.path.dirname(output.screen_sjc_out)


def iris_screen_sjc_out_file_name():
    run_name = config['run_name']
    splice_type = config['splice_event_type']
    name = 'SJ.{}.{}.summary_by_sig_event.txt'.format(run_name, splice_type)
    return os.path.join(result_dir(), 'screen_sjc', name)


rule iris_screen_sjc:
    input:
        parameter_file=os.path.join(result_dir(), 'screen.para'),
        splice_txt=splicing_matrix_txt_path_for_run(),
        splice_idx=splicing_matrix_idx_path_for_run(),
        sjc_count_txt=sjc_count_txt_path_for_run(),
        sjc_count_idx=sjc_count_idx_path_for_run(),
    output:
        screen_sjc_out=iris_screen_sjc_out_file_name(),
    log:
        out=os.path.join(result_dir(), 'iris_screen_sjc_log.out'),
        err=os.path.join(result_dir(), 'iris_screen_sjc_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        splice_event_type=config['splice_event_type'],
        out_dir=iris_screen_sjc_out_dir_param,
    resources:
        mem_mb=config['iris_screen_sjc_mem_gb'] * 1024,
        time_hours=config['iris_screen_sjc_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS screen_sjc'
        ' --parameter-file {input.parameter_file}'
        ' --splicing-event-type {params.splice_event_type}'
        ' --event-list-file {input.splice_txt}'
        ' --outdir {params.out_dir}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_append_sjc_out_dir_param(wildcards):
    input = iris_append_sjc_input(wildcards)
    return os.path.dirname(input['event_list'])


def iris_append_sjc_event_list_for_tier(tier):
    screen_out_files = iris_screen_out_files()
    if tier not in ['tier1', 'tier2tier3']:
        raise Exception('iris_append_sjc_event_list_for_tier({}): unexpected tier'
                        .format(tier))

    return screen_out_files.get(tier)


def iris_append_sjc_out_file_name_for_tier(tier):
    run_name = config['run_name']
    splice_type = config['splice_event_type']
    name = '{}.{}.{}.txt.ijc_info.txt'.format(run_name, splice_type, tier)
    return os.path.join(result_dir(), 'screen', name)


def iris_append_sjc_out_file_name_with_tier_wildcard():
    run_name = config['run_name']
    splice_type = config['splice_event_type']
    name = '{}.{}.{{tier}}.txt.ijc_info.txt'.format(run_name, splice_type)
    return os.path.join(result_dir(), 'screen', name)


def iris_append_sjc_input(wildcards):
    inputs = dict()
    predict_tasks_done = iris_predict_task_done_file_names()
    inputs['predict_tasks_done'] = predict_tasks_done
    inputs['parameter_file'] = os.path.join(result_dir(), 'screen.para')
    inputs['screen_sjc_out'] = iris_screen_sjc_out_file_name()
    inputs['event_list'] = iris_append_sjc_event_list_for_tier(wildcards.tier)

    epitope_post_out_files = iris_epitope_post_out_files()
    junction_key = '{}_junction'.format(wildcards.tier)
    inputs['epitope_post_junction'] = epitope_post_out_files[junction_key]
    return inputs


rule iris_append_sjc:
    input:
        unpack(iris_append_sjc_input),
    output:
        append_sjc_out=iris_append_sjc_out_file_name_with_tier_wildcard(),
    log:
        out=os.path.join(result_dir(), 'iris_append_sjc_{tier}_log.out'),
        err=os.path.join(result_dir(), 'iris_append_sjc_{tier}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        splice_event_type=config['splice_event_type'],
        out_dir=iris_append_sjc_out_dir_param,
    resources:
        mem_mb=config['iris_append_sjc_mem_gb'] * 1024,
        time_hours=config['iris_append_sjc_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS append_sjc'
        ' --sjc-summary {input.screen_sjc_out}'
        ' --splicing-event-type {params.splice_event_type}'
        ' --outdir {params.out_dir}'
        ' --add-ijc-info'  # runs IRIS annotate_ijc within IRIS append_sjc
        ' --parameter-file {input.parameter_file}'
        ' --screening-result-event-list {input.event_list}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def iris_visual_summary_input(wildcards):
    inputs = dict()
    inputs['parameter_file'] = os.path.join(result_dir(), 'screen.para')
    splice_type = config['splice_event_type']
    tier1_dir = os.path.join(result_dir(), 'screen',
                             '{}.tier1'.format(splice_type))
    if has_tier_1():
        inputs['tier1_peptide'] = os.path.join(
            tier1_dir, 'epitope_summary.peptide-based.txt')

    tier2tier3_dir = os.path.join(result_dir(), 'screen',
                                  '{}.tier2tier3'.format(splice_type))
    if has_tier_3():
        inputs['tier2tier3_peptide'] = os.path.join(
            tier2tier3_dir, 'epitope_summary.peptide-based.txt')

    return inputs


rule iris_visual_summary:
    input:
        unpack(iris_visual_summary_input),
    output:
        summary=os.path.join(result_dir(), 'visualization', 'summary.png'),
    log:
        out=os.path.join(result_dir(), 'visualization',
                         'iris_visual_summary_log.out'),
        err=os.path.join(result_dir(), 'visualization',
                         'iris_visual_summary_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        conda_env_2=config['conda_env_2'],
        splice_event_type=config['splice_event_type'],
        screen_dir=os.path.join(result_dir(), 'screen'),
    resources:
        mem_mb=config['iris_visual_summary_mem_gb'] * 1024,
        time_hours=config['iris_visual_summary_time_hr'],
    shell:
        '{params.conda_wrapper} {params.conda_env_2} IRIS visual_summary'
        ' --parameter-fin {input.parameter_file}'
        ' --screening-out-dir {params.screen_dir}'
        ' --out-file-name {output.summary}'
        ' --splicing-event-type {params.splice_event_type}'
        ' 1> {log.out}'
        ' 2> {log.err}'
