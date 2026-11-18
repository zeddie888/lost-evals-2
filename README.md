# LOST Evaluation v2

This repo is an effort to redesign [lost-evals](https://github.com/UWCubeSat/lost-evals)

## Setup

### Step 1: Compile and symlink LOST:
1. You must have [lost](https://github.com/UWCubeSat/lost) installed on your machine in a `lost` folder
2. Symlink LOST into the top level directory of this project (`lost-evals-2`):
```
ln -s <absolute path of lost folder> lost
```
3. To support Tetra, you will have to be in the right git branch since Tetra is not yet merged to master. Do `git checkout ezhang8-tetra2` 


### Step 2: Download Python packages

`lost-evals-2` shouldn't care about your Python version (I developed this in 3.12.3) but I'd still recommend running `lost-evals-2` in a Python `venv` since you need to install the Python packages in `requirements.txt`.
1. To create a virtual environment, run `python3 -m venv <name of venv>` in the top-level `lost-evals-2` directory, I just do `python3 -m venv venv`
    - You may need to install `venv`. On Ubuntu, do `sudo apt update && sudo apt install python3-venv`
    - Note: creating the venv might take a while (up to a minute or two?)
    - Step 1 should create a folder called `venv` (or whatever you named the venv) in `lost-evals`. You should never need to manually edit any files inside of it
2. To start the venv, run `source <venv name>/bin/activate`, e.g. `source venv/bin/activate`
3. To install required Python packages, run `pip install -r requirements.txt`

## Usage

`main.py` is currently here to give you a basic way of interacting with LOST through `lost_api`
1. Use or add a function to `main.py`, e.g. `test_run_entire_pipeline_random_attitude` is there for you already
    - Note: if you add a function, it must begin with `test_` for now
2. In your terminal (from the `src` folder), invoke the function with `python3 main.py <test function name>`, e.g. `python3 main.py test_run_entire_pipeline_random_attitude`

## Development

If you want to make edits to this repo, I recommend the following setup:
- VSCode
- Use `autopep8` as your Python formatter, install this VSCode extension: [link](https://marketplace.visualstudio.com/items?itemName=ms-python.autopep8)

## Paired Tetra / Pyramid pilot

Build LOST commit `a2303e6080bcbb2763596c26620f267d6c9bcd11` in release mode.
From this repo root, run:

```sh
python3 src/paired.py --lost-dir /absolute/path/to/built/lost --output pilot
python3 src/paired.py --lost-dir /absolute/path/to/built/lost --manifest pilot/cases.json --output replay
python3 -m unittest discover -s tests
```

The runner uses Python 3.10 or later on macOS or Linux. It needs no extra
Python packages. Each run builds both databases from the same catalog and
magnitude cutoff. Choose a new output directory for each run.

By default, the runner tests 100 clean fields, the same 100 fields with five
false detections added to each, and 10 fields containing only false detections.
Both algorithms receive the same ordered star positions and labels. A correct
solve has at least four correct IDs and no wrong IDs. A wrong ID counts as a
wrong solve; fewer than four correct IDs with no wrong IDs counts as no solve.
The runner logs errors separately and excludes pairs it cannot verify.

Open `summary.md` for the comparison, `results.csv` for individual results, and
`run.json` for the database settings and file hashes. Use `cases.json` to replay
the same fields. The star ID timer excludes startup and database loading.
Peak memory includes the database and simulator. The target commit in
`run.json` records which source to build; the binary hash identifies the
executable used for the run.

### Run on a BeagleBone Black

Use a Debian installation with Python 3.10 or newer. Build LOST natively on
the board from the commit above; the macOS executable will not run on it.
For a Debian build, install `g++ make groff xxd libcairo2-dev libeigen3-dev`
and build with `make -j1 release LOST_DISABLE_ASAN=1`. Use one build job to
limit memory use. Check free storage and RAM before building the databases.

Copy this benchmark and the desktop run's `cases.json` to the board. First
run the command above with `--fields 1` and a fresh output directory. Check
that both algorithms pass the known case. Then use `--manifest` with the
desktop case list and use a separate output directory for the board run.
The benchmark needs no camera.

Keep the raw input files from both machines. A shared seed does not guarantee
identical inputs across C++ libraries or architectures. Compare the ordered
coordinates and catalog labels in the per-case `.input` files before treating
the numbers as a direct hardware comparison. If they differ, use stored
centroids with a replay driver before drawing hardware comparisons.
