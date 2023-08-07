import os
import sys
import logging
import multiprocessing
import subprocess
import click

from snakemake.io import load_configfile

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context

def cli(obj):
    """
    RePaRe - workflow for iterative reconstruction of pangenomes.

    For updates or issues, refer to:
    """

cli.add_command(run_init)
cli.add_command(run_init_sra)

def get_snakefile_1(file="./snakefile_optimized_part1"):
    sf1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf1):
        sys.exit("Unable to locate 'snakefile_optimized_part1'; tried %s" % sf1)
    return sf1

def get_snakefile_2(file="./snakefile_optimzized_part2"):
    sf2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf2):
        sys.exit("Unable to locate 'snakefile_optimized_part2'; tried %s" % sf2)
    return sf2

@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True),
    short_help="run RePaRe workflow",
)
@click.argument(
    "workflow",
    type=click.Choice(
        [
            "parallel",
            "iterative",
            "None",
            "all",
        ]
    ),
)
@click.option(
    "-w",
    "--working-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="location to run RePaRe",
    default=".",
)
@click.option(
    "-c",
    "--config-file",
    type=click.Path(exists=True, resolve_path=True),
    help="config file to use",
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    default=multiprocessing.cpu_count(),
    show_default=True,
    help="maximum number of jobs to run in parallel",
)
@click.option(
    "-@",
    "--max-mem",
    type=float,
    default=None,
    help="Maximum memory to use",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)

def run_workflow_1(
    workflow, working_dir, config_file, jobs, max_mem, dryrun, snakemake_args
):
    """
    Runs RePaRe pipeline

    By default runs both parts of the pipeline (parallel + iterative), but one of the two can be specified.
    """

    if config_file is None:
        config_file = os.path.join("config.yaml") #insert path to default config file
    
    if not os.path.exists(config_file):
        logger(
            f"config_file not found: {config_file}\n" "please check the directory"
        )
        exit(1)
    
    validate_config(config_file, workflow)

    conf = load_configfile(config_file)

    samples_file = conf["samples_file"]

    cmd = (
        "snakemake --snakefile {snakefile} --directory {working_dir} "
        " --rerun-triggers mtime {jobs} --rerun-incomplete "
        "--configfile '{config_file}' "
        " {dryrun} --scheduler greedy "
        " {target_rule} "
        " {args} "
    ).format(
        snakefile = get_snakefile_1(),
        working_dir = working_dir,
        jobs = "--jobs {}".format(jobs) if jobs is not None else "",
        config_file = config_file,
        dryrun = "--dryrun" if dryrun else "",
        args = " ".join(snakemake_args),
        target_rule = workflow if workflow != "None" else "", 
    )
    logger.debug("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e)
        exit(1)

def run_workflow_2(
    workflow, working_dir, config_file, jobs, max_mem, dryrun, snakemake_args
):

    if config_file is None:
        config_file = os.path.join("config.yaml") #insert path to default config file
    
    if not os.path.exists(config_file):
        logger(
            f"config_file not found: {config_file}\n" "please check the directory"
        )
        exit(1)
    
    validate_config(config_file, workflow)

    conf = load_configfile(config_file)

    samples_file = conf["samples_file"]

    cmd = (
        "snakemake --snakefile {snakefile} --directory {working_dir} "
        " --rerun-triggers mtime {jobs} --rerun-incomplete "
        "--configfile '{config_file}' "
        " {dryrun} --scheduler greedy "
        " {target_rule} "
        " {args} "
    ).format(
        snakefile = get_snakefile_2(),
        working_dir = working_dir,
        jobs = "--jobs {}".format(jobs) if jobs is not None else "",
        config_file = config_file,
        dryrun = "--dryrun" if dryrun else "",
        args = " ".join(snakemake_args),
        target_rule = workflow if workflow != "None" else "", 
    )
    logger.debug("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e)
        exit(1)

if __name__ == "__main__":
    cli()