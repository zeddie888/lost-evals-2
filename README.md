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

The standalone runner needs Python 3.10+ on macOS or Linux and no Python packages.
It invokes LOST directly because the existing wrapper renders plots and does not
expose the deterministic, centroid-only and timing options needed here.
Each run builds and times separate databases from the same catalog and magnitude
cutoff. Use a fresh output directory; nothing is written into the LOST checkout.

The pilot verifies a known smoke case and compares 100 clean fields, their
five-false-detection counterparts, and 10 false-only controls at 12° horizontal
FOV. LOST generates false stars across the sky; the runner selects seeds with
exactly five detections in-frame and checks ordered inputs and labels across
algorithms. Clean and stressed versions use identical camera attitudes. With
no centroid perturbation, comparison uses an exact position threshold. A
correct solve requires four or more correct IDs and no wrong IDs; failures to
execute are separate from algorithm no-solves.

Start with `summary.md`, then refer to `results.csv`, `cases.json`, `run.json`
and the per-case logs under `raw/`. The stage timer excludes startup and
database loading. Peak RSS is measured serially in fresh processes via
`getrusage` and includes the database and simulator allocations. The catalog
and binary hashes identify inputs; the reported target commit is the expected
source revision, not a verification of how the executable was built.
