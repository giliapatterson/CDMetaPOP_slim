import argparse
import datetime,time,os,sys
import numpy as np
import subprocess
import polars as pl
from multiprocessing import Pool, Process, Queue, Lock
import multiprocessing as mp
import platform


def run_slim(run_df, sim_info_q, outdir, args, script_dir):
    finished = 0
    already_run = 0
    failed = 0
    # Helper: wrap strings in single quotes for SLiM and normalize separators on Windows
    def q(v):
        return f"'{str(v).replace('\\', '/')}'"
    for row in run_df.iter_rows(named=True):
        rep_start = time.perf_counter()
        rep_output_folder = os.path.join(outdir, f"run{row['run']-1}batch0mc{row['rep']}species0")
        simulation_finished = os.path.join(rep_output_folder, f"finished_{row['seed']}.txt")
        if os.path.exists(simulation_finished) and not args.rerun:
            already_run += 1
            print(f"Run {row['run']} rep {row['rep']} with seed {row['seed']} already finished, skipping. Results: {rep_output_folder}")
        else:
            if not os.path.exists(rep_output_folder):
                os.makedirs(rep_output_folder, exist_ok=True)
            if platform.system() != 'Windows':
                command = ["slim",
                        "-d", f"SEED={row['seed']}",
                        "-d", f"RUNVARS_FILE='{row['runvars']}'",
                        "-d", "PARAM_FOLDER='/'",
                        "-d", f"IND_OUT_FOLDER='{rep_output_folder}'",
                        "-d", f"ALLPOPS_OUT='{os.path.join(rep_output_folder, 'summary_popAllTime.csv')}'",
                        "-d", f"BYCLASS='{os.path.join(rep_output_folder, 'summary_classAllTime.csv')}'",
                        "-d", f"QTL_OUT='{os.path.join(rep_output_folder, 'QTL_overall.csv')}'",
                        "-d", f"QTL_SUBPOPS_OUT='{os.path.join(rep_output_folder, 'QTL_subpops.csv')}'",
                        "-d", f"FINISHED='{simulation_finished}'",
                        os.path.join(script_dir, "scripts", "cdmetapop_slim.slim")]
            else:
                    command = ["slim",
                        "-d", f"SEED={row['seed']}",
                        "-d", f"RUNVARS_FILE={q(row['runvars'])}",
                        "-d", f"PARAM_FOLDER={q('')}",
                        "-d", f"IND_OUT_FOLDER={q(rep_output_folder)}",
                        "-d", f"ALLPOPS_OUT={q(os.path.join(rep_output_folder, 'summary_popAllTime.csv'))}",
                        "-d", f"BYCLASS={q(os.path.join(rep_output_folder, 'summary_classAllTime.csv'))}",
                        "-d", f"QTL_OUT={q(os.path.join(rep_output_folder, 'QTL_overall.csv'))}",
                        "-d", f"QTL_SUBPOPS_OUT={q(os.path.join(rep_output_folder, 'QTL_subpops.csv'))}",
                        "-d", f"FINISHED={q(simulation_finished)}",
                        os.path.join(script_dir, "scripts", "cdmetapop_slim.slim")]
            stdout_file = os.path.join(rep_output_folder, "log.txt")
            with open(stdout_file, 'w') as f:
                output = subprocess.run(command, stdout=f, stderr=subprocess.STDOUT)
            if output.returncode != 0:
                failed += 1
                print(f"Run {row['run']} rep {row['rep']} failed in {time.perf_counter() - rep_start:.2g} seconds\nLog: {stdout_file}\nCommand: {' '.join(command)}")
            else:
                finished += 1
                print(f"Run {row['run']} rep {row['rep']} finished in {time.perf_counter() - rep_start:.2g} seconds. Results: {rep_output_folder}")
    sim_info = {'finished': finished, 'failed': failed, 'already_run': already_run}
    sim_info_q.put(sim_info)
    return (finished, failed, already_run)


def main():
    start_time = datetime.datetime.now()
    foldertime = int(time.time())

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--datadir", help="Data directory: directory that contains the CDMetaPOP RunVars file.")
    parser.add_argument("-i", "--inputfilename", help="RunVars file name. Should be specified relative to datadir.")
    parser.add_argument("-o", "--outputdir", default="cdmetapop_slim_output", help="Directory for simulation outputs. Should be specified relative to datadir.")
    parser.add_argument("-s", "--seed", help="Random seed")
    parser.add_argument('--filetime', action='store_true')
    parser.add_argument('--no-filetime', dest='filetime', action='store_false')
    parser.set_defaults(filetime=True)
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument('--no-rerun', dest='rerun', action='store_false')
    parser.set_defaults(rerun=False)
    parser.add_argument('-c', '--cores', type=int, default=1)
    args = parser.parse_args()

    if args.datadir and args.inputfilename and args.outputdir:
        datadir = args.datadir
        fileans = os.path.join(datadir, args.inputfilename)
        if (args.filetime):
            outdir = os.path.join(datadir, args.outputdir + str(foldertime))
        else:
            outdir = os.path.join(datadir, args.outputdir)
        slim_params = os.path.join(outdir, "slim_parameters")
    else:
        print("User must specify data directory and input file name, e.g., at command line type CDmetaPOP_slim.py -d ../CDmetaPOP_data/ -i RunVars.csv")
        sys.exit(-1)

    # If .ip file does not exist
    if not os.path.exists(fileans):
        print(("Cannot find or open runtime inputs file(%s)" % (fileans)))
        sys.exit(-1)

    # Create output file directory - will automatically put in the data directory
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # Directory of current file
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # ----------------------------------------
    # Convert CDMetaPOP input files to SLiM input files
    # ----------------------------------------

    print("---------------------------------------------------")
    print("CONVERTING PARAMETER FILES TO SLIM FORMAT")

    print(f"Processing {fileans} to SLiM format and writing to the directory {slim_params}")

    time1 = time.perf_counter()

    # First check if the parameters have already been processed
    if not os.path.exists(slim_params):
        os.makedirs(slim_params, exist_ok=True)
    command = ["Rscript", os.path.join(script_dir, "scripts", "make_slim_params.R"),
               "--parameter_directory", datadir,
               "--runvars_file_name", fileans,
               "--output_directory", slim_params]
    try:
        retcode = subprocess.call(command)
        if retcode < 0:
            print("Child was terminated by signal", -retcode, file=sys.stderr)
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)
        print(f"Return code: {retcode}")
    if retcode != 0:
        print(f"Processing parameter files did not finish successfully.")
        sys.exit(-1)

    time2 = time.perf_counter()

    print(f'Processing parameter files finished in {time2 - time1:.2g} seconds.')
    print("---------------------------------------------------")

    # ----------------------------------------
    # Run SLiM simulations
    # ----------------------------------------
    print(f"RUNNING SLIM SIMULATIONS")

    time3 = time.perf_counter()

    # old runvars
    old_runvars = pl.read_csv(fileans)

    if (args.seed):
        rng = np.random.default_rng(int(args.seed))
    else:
        rng = np.random.default_rng()

    # Make a dataframe with the parameters for each mcrun of each run
    # Call each mcrun a rep
    nreps = sum(old_runvars[:, 'mcruns'])
    nruns = len(old_runvars)
    rep_df = pl.DataFrame({'rep': np.array([np.arange(a) for a in old_runvars[:, 'mcruns']]).flatten(),
                           'run': np.repeat(np.arange(len(old_runvars)) + 1, old_runvars[:, 'mcruns'], axis=0),
                           'seed': rng.integers(low=1, high=30000, size=nreps)})
    rep_df = rep_df.with_columns(param_folder=slim_params + os.sep + "run" + pl.col('run').cast(pl.String))
    rep_df = rep_df.with_columns(runvars=pl.col('param_folder') + os.sep + "RunVars_slim.csv")
    rep_df.write_csv(os.path.join(outdir, "simulation_info.csv"))

    # Determine number of cores to use
    # This is the minimum of args.cores, number of replicates to run, and number of cores available
    cores_available = mp.cpu_count()
    cores = min(args.cores, cores_available, nreps)
    print(f"Using {cores} cores out of {cores_available} available.")

    # Chop dataframe into batches and run each batch
    start_time = time.perf_counter()
    runs = []
    breakpoints = np.linspace(0, rep_df.shape[0], cores + 1).astype(int)

    sim_info_queue = Queue()
    for start, end in zip(breakpoints[:-1], breakpoints[1:]):
        subset_df = rep_df[start:end]
        runs.append(Process(target=run_slim, args=(subset_df, sim_info_queue, outdir, args, script_dir)))
    # Now Start Processes
    for i in range(len(runs)):
        runs[i].start()
    # Now Join Processes
    for i in range(len(runs)):
        runs[i].join()
    finished = 0
    failed = 0
    already_run = 0
    while not sim_info_queue.empty():
        info = sim_info_queue.get()
        finished += info['finished']
        failed += info['failed']
        already_run += info['already_run']
    sim_info_queue.close()
    print(f"{len(rep_df)} reps complete in {time.perf_counter() - start_time:.2g} seconds.\nReps: Failed {failed}/{len(rep_df)} | Already run {already_run}/{len(rep_df)} | Finished {finished}/{len(rep_df)}")
    print("---------------------------------------------------")


if __name__ == '__main__':
    start_method = 'fork' if platform.system() != 'Windows' else 'spawn'
    mp.set_start_method(start_method)
    main()