"""A third node for run_all.py (parts 0 and 1 of 2 run on euclid24/25): takes jobs from the END of the list, and only
a job whose folder does not exist yet and which its own part reaches no sooner than two jobs from now (a run takes
~16 min, so the part can not start it while it runs here; when the part reaches it, runs.py sees it finished and
skips it). Stops when no such job is left.
    python3 run_helper.py
"""
import os

import run_all

NPARTS = 2
JOBS = [(c, m) for c in run_all.COSMOLOGIES for m in run_all.MODELS]


def done():
    with open(os.path.join(run_all.HERE, "progress.txt")) as f:
        return {line.split()[1] for line in f if line[0].isdigit()}


def current_position(part, finished):
    """The position (in the list of its part) of the job the part runs now: its first job not finished."""
    mine = JOBS[part::NPARTS]
    return next((i for i, (c, m) in enumerate(mine) if f"{c}/{m}" not in finished), len(mine))


if __name__ == "__main__":
    while True:
        finished = done()
        job = None
        for index in reversed(range(len(JOBS))):
            c, m = JOBS[index]
            part, position = index % NPARTS, index // NPARTS
            if os.path.exists(os.path.join(run_all.HERE, c, m)) or f"{c}/{m}" in finished:
                continue
            if position >= current_position(part, finished) + 2:
                job = (c, m)
            break  # Only the last job not started: never jump over a job, so the parts and this meet once
        if job is None:
            break
        print("helper:", job, flush=True)
        run_all.run(job)
    with open(os.path.join(run_all.HERE, "progress.txt"), "a") as f:
        f.write("ALLDONE helper\n")
