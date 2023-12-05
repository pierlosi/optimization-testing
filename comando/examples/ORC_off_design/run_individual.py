from subprocess import run
import multiprocessing

nb_cores = multiprocessing.cpu_count()

# start_design = 'results_263.15_268.15_273.15_278.15_283.15_288.15_293.15_298.15_303.15_308.15_313.15.pickle'

# with fixed mean design
# start='mean_design/results_289.pickle' fix=True -subdir=mean_design -time=86400


def solve(T, opts):

    TK = T + 273.15
    def write_opt(k, v):
        if isinstance(v, bool):
            return f'-{k}' if v else ''
        else:
            return f'-{k}={v}'
    str_opts = (write_opt(k, v) for k, v in opts.items() if v is not None)
    cmd = f"python case_study.py -sup -rec -od -T={TK} {' '.join(str_opts)}"
    print(f"starting optimization for T = {TK} K ({cmd})\n")

    cp = run(cmd.split(), capture_output=True)
    if cp.stderr:
        print(f"\033[91moptimization for T = {TK} K gave an error:\n{cp.stderr.decode()}\033[0m\n")
    else:
        print(f"optimization for T = {TK} K completed successfully!\n")
        # print(cp.stdout.decode())
    return cp.returncode


if __name__ == '__main__':
    import argparse

    import numpy as np

    ap = argparse.ArgumentParser()
    # ap.add_argument("-tol", "-t", required=False, default=1e-2,
    #                 type=float, help="absolute/relative optimality tolerance")
    ap.add_argument("-maxTime", "-time", required=False, default=10,
                    type=float, help="maximum time")
    # ap.add_argument("-branching_priorities", "-bp", required=False, default=0,
    #                 type=int, help="select strategy for branching priorities")
    # ap.add_argument("-use_superheater", "-sup", required=False, default=False,
    #                 action='store_true', help="whether to use a superheater")
    # ap.add_argument("-use_recuperator", "-rec", required=False, default=False,
    #                 action='store_true', help="whether to use a recuperator")
    # ap.add_argument("-plot", "-p", required=False, default=False,
    #                 action='store_true', help="whether to create a plot of "
    #                 "the subexpression graph of the resulting problem")
    ap.add_argument("-subdir", required=False, default='.',
                    type=str, help="subdirectory in which to store results")
    # ap.add_argument("-obj", required=False, default='TAC',
    #                 type=str, help="objective")
    ap.add_argument("-start_from", "-start", required=False, default=None,
                    type=str, help="Pickle file with start values")
    ap.add_argument("-fix_design", "-fix", required=False, default=False,
                    action='store_true', help="whether to fix design")
    ap.add_argument("-T", required=False, default='-10_41_5',
                    type=str, help="<start>_<stop> or <start>_<stop>_<step> "
                    "temperatures")
    ap.add_argument("-processes", "-n", required=False,
                    default=max(1, nb_cores-1), type=int,
                    help="the number of parallel processes to run")

    opts = vars(ap.parse_args())
    processes = opts.pop('processes')
    TEMPS = np.arange(*map(float, opts.pop('T').split('_')))
    # TEMPS = range(-10, 41, 5)
    # TEMPS = range(-10, 41)
    # TEMPS = [10, 15]
    # TEMPS = range(20, 41, 5)

    print('Running optimizations for temperatures:')
    print('   ', [*TEMPS])
    print(f'using {processes} parallel processes, with the following options:')
    for opt, val in opts.items():
        print(f'    {opt}: {val}')

    import time
    time.sleep(5)

    args = zip(TEMPS, [opts] * len(TEMPS))

    # for Windows, also needs to be within __name__ == '__main__' check
    multiprocessing.freeze_support()

    with multiprocessing.Pool(processes=processes) as pool:
        results = pool.starmap(solve, args)
