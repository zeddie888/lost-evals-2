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


def test_centroiding_step():
    """Run the pipeline and request that the input centroids be printed to a file,
    then verify the file was created and is non-empty."""
    from lost_api import run_entire_pipeline, PipelineOptions
    import os

    centroid_fname = "input_centroids.txt"
    actual_fname = "actual_centroids.txt"
    options = PipelineOptions(
        generate=1,
        generate_random_attitudes=True,
        generate_exposure=0.6,
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude="attitude.txt",
        print_input_centroids=centroid_fname,
        print_actual_centroids=actual_fname,
    )

    ok = run_entire_pipeline(options)
    if not ok:
        print("Pipeline failed during centroiding step")
        return

    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    lost_dir = os.path.join(project_root, "lost")
    path = os.path.join(lost_dir, centroid_fname)
    if os.path.isfile(path) and os.path.getsize(path) > 0:
        print(f"Centroid output written: {path}")
    else:
        print(f"Centroid output missing or empty: {path}")

    # Try to compute average centroiding accuracy by comparing input vs actual centroids.
    lost_dir = os.path.join(project_root, "lost")
    input_path = os.path.join(lost_dir, centroid_fname)
    actual_path = os.path.join(lost_dir, actual_fname)
    avg_error = None
    try:
        if os.path.isfile(input_path) and os.path.isfile(actual_path):
            import re
            def parse_points(p):
                points = []
                num_re = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
                for line in p:
                    toks = num_re.findall(line)
                    if len(toks) >= 2:
                        x = float(toks[0])
                        y = float(toks[1])
                        points.append((x, y))
                return points

            with open(input_path, "r") as fi, open(actual_path, "r") as fa:
                inp_pts = parse_points(fi)
                act_pts = parse_points(fa)

            n = min(len(inp_pts), len(act_pts))
            if n > 0:
                import math
                errs = [math.hypot(inp_pts[i][0] - act_pts[i][0], inp_pts[i][1] - act_pts[i][1]) for i in range(n)]
                avg_error = sum(errs) / n
                print(f"Average centroid error (pixels) over {n} matches: {avg_error}")
            else:
                print("No matching centroid points found to compute accuracy")
        else:
            print("Input or actual centroid files missing; cannot compute centroid accuracy")
    except Exception as e:
        print(f"Error computing centroid accuracy: {e}")

    return avg_error


def test_star_id_step():
    """Run the pipeline configured to produce star-id comparison output and
    check the comparison file exists and is non-empty."""
    from lost_api import run_entire_pipeline, PipelineOptions
    import os

    compare_fname = "compare_star_ids.txt"
    options = PipelineOptions(
        generate=1,
        generate_random_attitudes=True,
        generate_exposure=0.6,
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude="attitude.txt",
        print_actual_centroids="actual_centroids.txt",
        compare_star_ids=compare_fname,
    )

    ok = run_entire_pipeline(options)
    if not ok:
        print("Pipeline failed during star-id step")
        return None

    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    lost_dir = os.path.join(project_root, "lost")
    path = os.path.join(lost_dir, compare_fname)
    if os.path.isfile(path) and os.path.getsize(path) > 0:
        print(f"Star-id comparison output written: {path}")
    else:
        print(f"Star-id comparison output missing or empty: {path}")

    # Determine number of detected stars by counting lines in actual centroids file if present
    actual_centroids_path = os.path.join(lost_dir, "actual_centroids.txt")
    try:
        if os.path.isfile(actual_centroids_path):
            with open(actual_centroids_path, "r") as f:
                # count non-empty lines that contain a number
                import re
                num_re = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
                count = 0
                for line in f:
                    if num_re.search(line):
                        count += 1
            print(f"Detected stars (from actual_centroids.txt): {count}")
            return count
        else:
            print("actual_centroids.txt not found; cannot count detected stars")
            return None
    except Exception as e:
        print(f"Error counting detected stars: {e}")
        return None


def test_attitude_step():
    """Run the pipeline and parse the attitude output file to verify attitude
    determination works (smoke test)."""
    from lost_api import run_entire_pipeline, parse_attitude_result, PipelineOptions
    import os

    attitude_fname = "attitude_step.txt"
    options = PipelineOptions(
        generate=1,
        generate_random_attitudes=True,
        generate_exposure=0.6,
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude=attitude_fname,
    )

    ok = run_entire_pipeline(options)
    if not ok:
        print("Pipeline failed during attitude step")
        return

    res = parse_attitude_result(attitude_fname)
    if res is None:
        print("Could not parse attitude output")
    elif res.known:
        print(f"Attitude determined: ra={res.ra}, de={res.de}, roll={res.roll}")
    else:
        print("Attitude was not determined for this run")


def test_ablation_study_fovs(runs_per_fov: str | int = 5):
    """Run an ablation study sweeping different FOVs and plot success rate.

    runs_per_fov: number of repeated pipeline runs to average for each FOV.
    The function will write a plot `ablation_fov.png` to the project root.
    """
    try:
        runs = int(runs_per_fov)
    except Exception:
        raise ValueError("runs_per_fov must be an integer or convertible to int")

    from lost_api import run_entire_pipeline, parse_attitude_result, PipelineOptions
    import matplotlib.pyplot as plt
    import os

    # FOV values to sweep (degrees)
    fovs = [5, 10, 15, 20, 25, 30]
    success_rates = []

    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    for fov in fovs:
        successes = 0
        for i in range(runs):
            attitude_fname = f"attitude_fov{fov}_run{i}.txt"
            options = PipelineOptions(
                generate=1,
                generate_random_attitudes=True,
                generate_exposure=0.6,
                centroid_algo="cog",
                database="tetra-20.dat",
                false_stars=0,
                star_id_algo="tetra",
                attitude_algo="dqm",
                centroid_mag_filter=0,
                print_attitude=attitude_fname,
                fov_deg=float(fov)
            )

            ok = run_entire_pipeline(options)
            if not ok:
                print(f"Pipeline run failed for FOV={fov} (run {i})")
                # treat as failure but continue
                continue

            res = parse_attitude_result(attitude_fname)
            if res is None:
                print(f"Could not parse attitude result for FOV={fov} (run {i})")
            elif res.known:
                successes += 1

        rate = successes / runs if runs > 0 else 0.0
        success_rates.append(rate)
        print(f"FOV={fov} deg: success {successes}/{runs} -> {rate:.2%}")

    # Plot results
    plt.figure(figsize=(8, 5))
    plt.plot(fovs, success_rates, marker='o')
    plt.title('Ablation study: attitude success rate vs FOV')
    plt.xlabel('FOV (degrees)')
    plt.ylabel('Success rate')
    plt.ylim(-0.05, 1.05)
    plt.grid(True)

    out_path = os.path.join(project_root, 'ablation_fov.png')
    plt.savefig(out_path, bbox_inches='tight')
    print(f"Saved ablation plot to: {out_path}")


def test_cross_boresight_accuracy(fov_deg: float = 20.0, x_res: int = 1024, y_res: int = 1024):
    """Compute cross-boresight accuracy using centroid and star-id tests.

    This function calls `test_centroiding_step` and `test_star_id_step` to
    obtain the average centroiding error (in pixels) and the number of
    detected stars, then calls `compute_cross_boresight_accuracy` from
    `lost_api` to compute the metric. Returns the computed metric or None
    if it cannot be computed.
    """
    # Run centroiding step and get average centroid error (pixels)
    print("Running centroiding step to compute average centroid accuracy...")
    avg_centroid_error = test_centroiding_step()
    if avg_centroid_error is None:
        print("Could not obtain average centroid accuracy; aborting cross-boresight computation")
        return None

    # Run star-id step and get detected star count
    print("Running star-id step to count detected stars...")
    detected_stars = test_star_id_step()
    if detected_stars is None or detected_stars <= 0:
        print("Could not obtain detected star count; aborting cross-boresight computation")
        return None

    from lost_api import compute_cross_boresight_accuracy

    num_pixels = int(x_res) * int(y_res)
    metric = compute_cross_boresight_accuracy(
        fov_deg=fov_deg,
        avg_centroid_accuracy=avg_centroid_error,
        num_pixels=num_pixels,
        avg_detected_stars=float(detected_stars),
    )

    print(f"Cross-boresight accuracy: {metric}")
    return metric


def test_around_boresight_accuracy(x_res: int = 1024, y_res: int = 1024):
    """Compute around-boresight (roll) accuracy using centroid and star-id tests.

    This function calls `test_centroiding_step` and `test_star_id_step` to
    obtain the average centroiding error (in pixels) and the number of
    detected stars, then calls `compute_around_boresight_accuracy` from
    `lost_api` to compute the roll error metric. Returns the computed metric
    (in radians) or None if it cannot be computed.
    """
    # Run centroiding step and get average centroid error (pixels)
    print("Running centroiding step to compute average centroid accuracy...")
    avg_centroid_error = test_centroiding_step()
    if avg_centroid_error is None:
        print("Could not obtain average centroid accuracy; aborting around-boresight computation")
        return None

    # Run star-id step and get detected star count
    print("Running star-id step to count detected stars...")
    detected_stars = test_star_id_step()
    if detected_stars is None or detected_stars <= 0:
        print("Could not obtain detected star count; aborting around-boresight computation")
        return None

    from lost_api import compute_around_boresight_accuracy
    import math

    num_pixels = int(x_res) * int(y_res)
    metric_rad = compute_around_boresight_accuracy(
        avg_centroid_accuracy=avg_centroid_error,
        num_pixels=num_pixels,
        avg_detected_stars=float(detected_stars),
    )

    # Convert to degrees and arcseconds for convenience
    metric_deg = metric_rad * 180.0 / math.pi
    metric_arcsec = metric_rad * 206265.0

    print(f"Around-boresight (roll) accuracy:")
    print(f"  {metric_rad:.6f} radians")
    print(f"  {metric_deg:.6f} degrees")
    print(f"  {metric_arcsec:.2f} arcseconds")
    return metric_rad


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
