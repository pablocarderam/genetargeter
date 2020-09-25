from multiprocessing import cpu_count

from numpy import floor


# from mkl import set_num_threads # mkl is never used outside of this


def configure(
    num_jobs: int = 8,
    test: bool = False,
    subtract: int = 1,
    num_proc: int = None,
    num_thread_per_proc: int = None,
) -> int:

    """
    num_jobs is typically the # of genes we are parallelizing over
    """
    if num_proc is None:
        num_proc = cpu_count() - subtract

    if num_jobs > num_proc:
        num_jobs = num_proc

    if num_thread_per_proc is None:
        num_thread_per_proc = int(floor(num_proc / num_jobs))

    if test:
        num_jobs = 1
        num_thread_per_proc = 1

    # try:
    #     set_num_threads(num_thread_per_proc)
    # except ImportError:
    #     print("MKL not available, so I'm not adjusting the number of threads")

    # print(f"Launching {num_jobs} jobs with {num_thread_per_proc} MKL threads each")

    return num_jobs
