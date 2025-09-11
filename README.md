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
