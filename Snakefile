import os
import shutil
from datetime import datetime
from glob import glob

from snakemake.utils import validate


configfile: 'config.yaml'
validate(config, schema='schema/config.schema.yaml')

timestamp = datetime.utcnow().strftime('%Y%m%dT%H%M%s')

#fetch parameters that are specifc to this minion run
INPUTDIR = config['INPUTDIR']
OUTPUTDIR = config['OUTPUTDIR']
TMPDIR = config['TMPDIR']
#TODO: add timestamp to tmp dir
# TMPDIR = os.path.join(TMPDIR, f'{timestamp}-miniontoolkit')
TARRED_READS = config['TARRED_READS']

FLOWCELL = config['FLOWCELL']
KIT = config['KIT']
ALBACORE_VERSION = config['ALBACORE_VERSION']
ALBACORE_THREADS = config['ALBACORE_THREADS']
ALBACORE_OUTPUT_FORMAT = config['ALBACORE_OUTPUT_FORMAT']

CLEANUP = config['CLEANUP']
DELETE_FAILED = config['DELETE_FAILED']


#directory path to Miniknow fast5 files are located
FAST5_DIR_PATHS = [os.path.abspath(f) for f in glob(os.path.join(INPUTDIR, '*.tar'))] \
    if TARRED_READS else \
    [os.path.abspath(f) for f in glob(os.path.join(INPUTDIR, '**')) if os.path.isdir(f)]

FAST5_DIRS = {os.path.basename(x).replace('.tar', '') if TARRED_READS else os.path.basename(x): x for x in FAST5_DIR_PATHS}

rule all:
    input:
        expand(os.path.join(OUTPUTDIR, 'progress', '{fast5_dir}.done'), fast5_dir=FAST5_DIRS.keys()),
        os.path.join(OUTPUTDIR, 'results', 'sequencing_summary.txt')


rule aggregate_results:
    input:
        expand(os.path.join(OUTPUTDIR, 'albacore_results', 'results_{fast5_dir}'), fast5_dir=FAST5_DIRS.keys())
    output:
        ss=os.path.join(OUTPUTDIR, 'results', 'sequencing_summary.txt'), #combine all albacore runs into a SINGLE sequencing_summary.txt
        health=os.path.join(OUTPUTDIR, 'results', 'run_health.txt'), # statistics on how well the run went.
    run:
        from scripts.gather_albacore_results import get_albacore_results

        get_albacore_results(input_dirs=input,
                             sequencing_summary=output.ss,
                             run_health=output.health,
                             deletefailed=DELETE_FAILED)

if TARRED_READS:
    rule albacore_on_tarred_reads:
        input:
            lambda wildcards: FAST5_DIRS[wildcards.fast5_dir]
        output:
            dir=directory(os.path.join(OUTPUTDIR, 'albacore_results', 'results_{fast5_dir}')),
            done=touch(os.path.join(OUTPUTDIR, 'progress', '{fast5_dir}.done'))
        threads: ALBACORE_THREADS
        benchmark: os.path.join(OUTPUTDIR, 'benchmarks', '{fast5_dir}-albacore.tsv')
        log:
            os.path.join(OUTPUTDIR, 'logs', '{fast5_dir}-albacore.log')
        conda:
            "envs/albacore-%s.yml" % (ALBACORE_VERSION)
        version:
            ALBACORE_VERSION
        shell:
            """
            mkdir -p {TMPDIR}/{wildcards.fast5_dir}
            tar -C {TMPDIR}/{wildcards.fast5_dir} -xf {input}
            read_fast5_basecaller.py \
              --input {TMPDIR}/{wildcards.fast5_dir} \
              --flowcell {FLOWCELL} --kit {KIT} --barcoding --recursive \
              --output_format {ALBACORE_OUTPUT_FORMAT} \
              --save_path {output.dir} \
              --disable_pings -q 999999999 --worker_threads {threads} > {log} 2>&1
            rm -r {TMPDIR}/{wildcards.fast5_dir}
            """
else:
    rule albacore:
        input:
            lambda wildcards: FAST5_DIRS[wildcards.fast5_dir]
        output:
            dir=directory(os.path.join(OUTPUTDIR, 'albacore_results', 'results_{fast5_dir}')),
            done=touch(os.path.join(OUTPUTDIR, 'progress', '{fast5_dir}.done'))
        threads: ALBACORE_THREADS
        benchmark: os.path.join(OUTPUTDIR, 'benchmarks', '{fast5_dir}-albacore.tsv')
        log:
            os.path.join(OUTPUTDIR, 'logs', '{fast5_dir}-albacore.log')
        conda:
            "envs/albacore-%s.yml" % (ALBACORE_VERSION)
        version:
            ALBACORE_VERSION
        shell:
            """
            read_fast5_basecaller.py \
              --input {input} \
              --flowcell {FLOWCELL} --kit {KIT} --barcoding --recursive \
              --output_format {ALBACORE_OUTPUT_FORMAT} \
              --save_path {output.dir} \
              --disable_pings -q 999999999 --worker_threads {threads} > {log} 2>&1
            """

onsuccess:
    if CLEANUP:
        #check to see if any results are there to be removed
        albacore_output_dir = os.path.join(OUTPUTDIR, 'albacore_results')
        if os.path.isdir(albacore_output_dir):
            print("Cleaning up individual albacore results as indicate by 'CLEANUP' flag.")
            shutil.rmtree(albacore_output_dir)
