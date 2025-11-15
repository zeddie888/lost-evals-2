import inspect
import sys

from lost_api import CompileOptions


def test_compile_lost():
    from lost_api import compile_lost

    options = CompileOptions(rebuild=False)
    compile_lost(options)


def test_clean_lost():
    from lost_api import clean_lost

    clean_lost()


def test_generate_tetra_database():
    from lost_api import generate_tetra_database, TetraDbOptions

    # NOTE: change this to match your image FOV
    options = TetraDbOptions(min_mag=7, tetra_max_angle=20, output="tetra-20.dat")
    generate_tetra_database(options)


def test_run_entire_pipeline_random_attitude():
    from lost_api import parse_attitude_result, PipelineOptions, run_entire_pipeline

    attitude_output_fname = "attitude.txt"
    options = PipelineOptions(
        generate=1,
        generate_random_attitudes=True,
        generate_exposure=1,
        fov=20,
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude=attitude_output_fname,
    )
    run_entire_pipeline(options)
    print(parse_attitude_result(attitude_output_fname))


def test_run_entire_pipeline_given_attitude(
    ra_input: str | float, de_input: str | float, roll_input: str | float
):
    try:
        ra = float(ra_input)
        de = float(de_input)
        roll = float(roll_input)
    except ValueError:
        raise ValueError("ra, de, and roll must be convertible to float")
    from lost_api import parse_attitude_result, PipelineOptions, run_entire_pipeline

    attitude_output_fname = "attitude.txt"
    options = PipelineOptions(
        generate=1,
        generate_ra=ra,
        generate_de=de,
        generate_roll=roll,
        generate_random_attitudes=False,
        generate_exposure=1,
        fov=20,
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude=attitude_output_fname,
    )
    run_entire_pipeline(options)
    print(parse_attitude_result(attitude_output_fname))


def _run_single_pipeline(run_id):
    """Helper function to run a single pipeline iteration."""
    import os

    from lost_api import parse_attitude_result, PipelineOptions, run_entire_pipeline

    # Use unique output file for each parallel run
    # File will be deleted after run is completed
    attitude_output_fname = f"attitude_{os.getpid()}_{run_id}.txt"

    options = PipelineOptions(
        generate=1,
        generate_random_attitudes=True,
        generate_exposure=1,
        fov=20,
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude=attitude_output_fname,
    )
    run_entire_pipeline(options)
    res = parse_attitude_result(attitude_output_fname)

    # Clean up the temporary file from the lost directory
    try:
        lost_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "lost")
        file_path = os.path.join(lost_dir, attitude_output_fname)
        if os.path.exists(file_path):
            os.remove(file_path)
    except:
        pass

    return run_id, res.known


def test_run_pipeline_n_times(n_input: str | int):
    try:
        n = int(n_input)
    except ValueError:
        raise ValueError("n must be convertible to int")

    if n <= 0:
        raise ValueError("n must be a positive integer")

    import os
    from concurrent.futures import as_completed, ProcessPoolExecutor

    # Use all available CPU cores
    max_workers = os.cpu_count()
    print(f"Running {n} iterations in parallel using {max_workers} workers...")

    successful_runs = 0
    completed = 0

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        futures = {executor.submit(_run_single_pipeline, i): i for i in range(n)}

        for future in as_completed(futures):
            run_id, is_successful = future.result()
            completed += 1

            if is_successful:
                successful_runs += 1

            print(f"Run {completed}/{n}: {'Success' if is_successful else 'Failed'}")

    print(f"\n{'='*50}")
    print(f"Results: {successful_runs}/{n} successful runs")
    print(f"Success rate: {(successful_runs/n)*100:.2f}%")
    print(f"{'='*50}")


if __name__ == "__main__":
    current_module = sys.modules[__name__]
    functions = {
        name: func
        for name, func in inspect.getmembers(current_module, inspect.isfunction)
        if name.startswith("test_")
    }

    if len(sys.argv) > 1:
        func_name = sys.argv[1]
        func = functions.get(func_name)
        if func:
            args = sys.argv[2:]
            func(*args)
        else:
            print(f"Function '{func_name}' not found.")
    else:
        print("Select a function to run:")
        for name in functions:
            print(f"- {name}")
