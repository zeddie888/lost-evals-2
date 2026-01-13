import inspect
import sys

from lost_api import CompileOptions, CentroidResult, StarIDResult


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
        fov=25,
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


# def test_centroiding_step():
#     """Run the pipeline and request that the input centroids be printed to a file,
#     then verify the file was created and is non-empty."""
#     from lost_api import run_entire_pipeline, PipelineOptions
#     import os

#     centroid_fname = "input_centroids.txt"
#     actual_fname = "actual_centroids.txt"
#     options = PipelineOptions(
#         generate=1,
#         generate_random_attitudes=True,
#         generate_exposure=0.6,
#         centroid_algo="cog",
#         database="tetra-20.dat",
#         false_stars=0,
#         star_id_algo="tetra",
#         attitude_algo="dqm",
#         centroid_mag_filter=0,
#         print_attitude="attitude.txt",
#         print_input_centroids=centroid_fname,
#         print_actual_centroids=actual_fname,
#     )

#     ok = run_entire_pipeline(options)
#     if not ok:
#         print("Pipeline failed during centroiding step")
#         return

#     project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#     lost_dir = os.path.join(project_root, "lost")
#     path = os.path.join(lost_dir, centroid_fname)
#     if os.path.isfile(path) and os.path.getsize(path) > 0:
#         print(f"Centroid output written: {path}")
#     else:
#         print(f"Centroid output missing or empty: {path}")

#     # Try to compute average centroiding accuracy by comparing input vs actual centroids.
#     lost_dir = os.path.join(project_root, "lost")
#     input_path = os.path.join(lost_dir, centroid_fname)
#     actual_path = os.path.join(lost_dir, actual_fname)
#     avg_error = None
#     try:
#         if os.path.isfile(input_path) and os.path.isfile(actual_path):
#             import re
#             def parse_points(p):
#                 points = []
#                 num_re = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
#                 for line in p:
#                     toks = num_re.findall(line)
#                     if len(toks) >= 2:
#                         x = float(toks[0])
#                         y = float(toks[1])
#                         points.append((x, y))
#                 return points

#             with open(input_path, "r") as fi, open(actual_path, "r") as fa:
#                 inp_pts = parse_points(fi)
#                 act_pts = parse_points(fa)

#             n = min(len(inp_pts), len(act_pts))
#             if n > 0:
#                 import math
#                 errs = [math.hypot(inp_pts[i][0] - act_pts[i][0], inp_pts[i][1] - act_pts[i][1]) for i in range(n)]
#                 avg_error = sum(errs) / n
#                 print(f"Average centroid error (pixels) over {n} matches: {avg_error}")
#             else:
#                 print("No matching centroid points found to compute accuracy")
#         else:
#             print("Input or actual centroid files missing; cannot compute centroid accuracy")
#     except Exception as e:
#         print(f"Error computing centroid accuracy: {e}")

#     return avg_error


# def test_star_id_step():
#     """Run the pipeline configured to produce star-id comparison output and
#     check the comparison file exists and is non-empty."""
#     from lost_api import run_entire_pipeline, PipelineOptions
#     import os

#     compare_fname = "compare_star_ids.txt"
#     options = PipelineOptions(
#         generate=1,
#         generate_random_attitudes=True,
#         generate_exposure=0.6,
#         centroid_algo="cog",
#         database="tetra-20.dat",
#         false_stars=0,
#         star_id_algo="tetra",
#         attitude_algo="dqm",
#         centroid_mag_filter=0,
#         print_attitude="attitude.txt",
#         print_actual_centroids="actual_centroids.txt",
#         compare_star_ids=compare_fname,
#     )

#     ok = run_entire_pipeline(options)
#     if not ok:
#         print("Pipeline failed during star-id step")
#         return None

#     project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#     lost_dir = os.path.join(project_root, "lost")
#     path = os.path.join(lost_dir, compare_fname)
#     if os.path.isfile(path) and os.path.getsize(path) > 0:
#         print(f"Star-id comparison output written: {path}")
#     else:
#         print(f"Star-id comparison output missing or empty: {path}")

#     # Determine number of detected stars by counting lines in actual centroids file if present
#     actual_centroids_path = os.path.join(lost_dir, "actual_centroids.txt")
#     try:
#         if os.path.isfile(actual_centroids_path):
#             with open(actual_centroids_path, "r") as f:
#                 # count non-empty lines that contain a number
#                 import re
#                 num_re = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
#                 count = 0
#                 for line in f:
#                     if num_re.search(line):
#                         count += 1
#             print(f"Detected stars (from actual_centroids.txt): {count}")
#             return count
#         else:
#             print("actual_centroids.txt not found; cannot count detected stars")
#             return None
#     except Exception as e:
#         print(f"Error counting detected stars: {e}")
#         return None


# def test_attitude_step():
#     """Run the pipeline and parse the attitude output file to verify attitude
#     determination works (smoke test)."""
#     from lost_api import run_entire_pipeline, parse_attitude_result, PipelineOptions
#     import os

#     attitude_fname = "attitude_step.txt"
#     options = PipelineOptions(
#         generate=1,
#         generate_random_attitudes=True,
#         generate_exposure=0.6,
#         centroid_algo="cog",
#         database="tetra-20.dat",
#         false_stars=0,
#         star_id_algo="tetra",
#         attitude_algo="dqm",
#         centroid_mag_filter=0,
#         print_attitude=attitude_fname,
#     )

#     ok = run_entire_pipeline(options)
#     if not ok:
#         print("Pipeline failed during attitude step")
#         return

#     res = parse_attitude_result(attitude_fname)
#     if res is None:
#         print("Could not parse attitude output")
#     elif res.known:
#         print(f"Attitude determined: ra={res.ra}, de={res.de}, roll={res.roll}")
#     else:
#         print("Attitude was not determined for this run")

def test_ablation_study_resolutions(runs_per_res: str | int = 50, fov: float = 25.0):
    """Run an ablation study sweeping different image resolutions and plot success rate.

    runs_per_res: number of repeated pipeline runs to average for each resolution.
    fov: field-of-view in degrees (kept constant across all resolutions).
    The function will write plots to the project root.
    """
    try:
        runs = int(runs_per_res)
    except Exception:
        raise ValueError("runs_per_res must be an integer or convertible to int")

    from lost_api import run_entire_pipeline, parse_attitude_result, PipelineOptions
    import matplotlib.pyplot as plt
    import os

    # Resolution values to sweep (width x height)
    resolutions = [256, 512, 1024, 2048]
    success_rates = []
    cross_boresight_accuracies = []
    around_boresight_accuracies = []

    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    for res in resolutions:
        successes = 0
        centroid_errors = []
        detected_stars_list = []
        
        print(f"\n{'='*60}")
        print(f"Testing resolution: {res}x{res}")
        print(f"{'='*60}")
        
        for i in range(runs):
            attitude_fname = f"attitude_res{res}_run{i}.txt"
            centroid_fname = f"input_centroids_res{res}_run{i}.txt"
            actual_centroid_fname = f"actual_centroids_res{res}_run{i}.txt"
            compare_star_ids_fname = f"compare_star_ids_res{res}_run{i}.txt"
            
            options = PipelineOptions(
                generate=1,
                generate_x_res=res,
                generate_y_res=res,
                generate_random_attitudes=True,
                generate_exposure=0.6,
                centroid_algo="cog",
                database="tetra-20.dat",
                false_stars=0,
                star_id_algo="tetra",
                attitude_algo="dqm",
                centroid_mag_filter=0,
                print_attitude=attitude_fname,
                print_input_centroids=centroid_fname,
                print_actual_centroids=actual_centroid_fname,
                compare_star_ids=compare_star_ids_fname,
                fov_deg=fov
            )

            ok = run_entire_pipeline(options)
            if not ok:
                print(f"Pipeline run failed for res={res} (run {i})")
                continue

            res_attitude = parse_attitude_result(attitude_fname)
            if res_attitude is None:
                print(f"Could not parse attitude result for res={res} (run {i})")
            elif res_attitude.known:
                successes += 1

            # Compute centroid error for this run
            lost_dir = os.path.join(project_root, "lost")
            input_path = os.path.join(lost_dir, centroid_fname)
            actual_path = os.path.join(lost_dir, actual_centroid_fname)
            avg_error = None
            try:
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
                if os.path.isfile(input_path) and os.path.isfile(actual_path):
                    with open(input_path, "r") as fi, open(actual_path, "r") as fa:
                        inp_pts = parse_points(fi)
                        act_pts = parse_points(fa)
                    n = min(len(inp_pts), len(act_pts))
                    if n > 0:
                        import math
                        errs = [math.hypot(inp_pts[j][0] - act_pts[j][0], inp_pts[j][1] - act_pts[j][1]) for j in range(n)]
                        avg_error = sum(errs) / n
                if avg_error is not None:
                    centroid_errors.append(avg_error)
                    if i % 10 == 0:  # Print every 10th run to reduce clutter
                        print(f"  Run {i}: centroid error = {avg_error:.4f} pixels")
            except Exception as e:
                print(f"Error computing centroid accuracy for res={res} run {i}: {e}")

            # Count detected stars for this run
            actual_centroids_path = actual_path
            try:
                if os.path.isfile(actual_centroids_path):
                    with open(actual_centroids_path, "r") as f:
                        import re
                        num_re = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
                        count = 0
                        for line in f:
                            if num_re.search(line):
                                count += 1
                    detected_stars_list.append(count)
                    if i % 10 == 0:  # Print every 10th run to reduce clutter
                        print(f"  Run {i}: detected stars = {count}")
            except Exception as e:
                print(f"Error counting detected stars for res={res} run {i}: {e}")

        rate = successes / runs if runs > 0 else 0.0
        success_rates.append(rate)
        print(f"\nResolution={res}x{res}: success {successes}/{runs} -> {rate:.2%}")

        # Compute average centroid error and detected stars for this resolution
        avg_centroid_error = sum(centroid_errors) / len(centroid_errors) if centroid_errors else None
        avg_detected_stars = sum(detected_stars_list) / len(detected_stars_list) if detected_stars_list else None
        
        if avg_centroid_error is not None:
            print(f"  Average centroid error: {avg_centroid_error:.4f} pixels")
        if avg_detected_stars is not None:
            print(f"  Average detected stars: {avg_detected_stars:.2f}")
        
        # Compute cross-boresight and around-boresight accuracy
        cross_acc_arcsec = None
        around_acc_arcsec = None
        if avg_centroid_error is not None and avg_detected_stars is not None and avg_detected_stars > 0:
            from lost_api import compute_cross_boresight_accuracy, compute_around_boresight_accuracy
            try:
                cross_acc = compute_cross_boresight_accuracy(
                    fov_deg=fov,
                    avg_centroid_accuracy=avg_centroid_error,
                    num_pixels=res,  # Use resolution, not total pixels
                    avg_detected_stars=avg_detected_stars,
                )
                # Convert cross-boresight accuracy from degrees to arcseconds
                cross_acc_arcsec = cross_acc * 3600.0
                
                around_acc = compute_around_boresight_accuracy(
                    avg_centroid_accuracy=avg_centroid_error,
                    num_pixels=res,  # Use resolution, not total pixels
                    avg_detected_stars=avg_detected_stars,
                )
                # Convert around-boresight accuracy from radians to arcseconds
                around_acc_arcsec = around_acc * 206265.0
                
                print(f"  Cross-boresight accuracy: {cross_acc_arcsec:.2f} arcseconds")
                print(f"  Around-boresight accuracy: {around_acc_arcsec:.2f} arcseconds")
            except Exception as e:
                print(f"Error computing accuracy metrics for res={res}: {e}")
        
        cross_boresight_accuracies.append(cross_acc_arcsec)
        around_boresight_accuracies.append(around_acc_arcsec)

    # Plot success rates
    plt.figure(figsize=(8, 5))
    plt.plot(resolutions, success_rates, marker='o', label='Success Rate', linewidth=2)
    plt.title(f'Ablation study: attitude success rate vs Resolution (FOV={fov}°)')
    plt.xlabel('Resolution (pixels)')
    plt.ylabel('Success rate')
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.xticks(resolutions, [f'{r}x{r}' for r in resolutions])

    out_path = os.path.join(project_root, 'ablation_resolution.png')
    plt.savefig(out_path, bbox_inches='tight', dpi=150)
    print(f"\nSaved success rate plot to: {out_path}")

    # Plot cross-boresight accuracy (arcseconds)
    plt.figure(figsize=(8, 5))
    plt.plot(resolutions, cross_boresight_accuracies, marker='o', 
             label='Cross-Boresight Accuracy', linewidth=2, color='tab:orange')
    plt.title(f'Ablation study: cross-boresight accuracy vs Resolution (FOV={fov}°)')
    plt.xlabel('Resolution (pixels)')
    plt.ylabel('Cross-boresight accuracy (arcseconds)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.xticks(resolutions, [f'{r}x{r}' for r in resolutions])
    
    out_path_cross = os.path.join(project_root, 'ablation_resolution_cross_boresight.png')
    plt.savefig(out_path_cross, bbox_inches='tight', dpi=150)
    print(f"Saved cross-boresight accuracy plot to: {out_path_cross}")

    # Plot around-boresight accuracy (arcseconds)
    plt.figure(figsize=(8, 5))
    plt.plot(resolutions, around_boresight_accuracies, marker='o', 
             label='Around-Boresight Accuracy', linewidth=2, color='tab:green')
    plt.title(f'Ablation study: around-boresight (roll) accuracy vs Resolution (FOV={fov}°)')
    plt.xlabel('Resolution (pixels)')
    plt.ylabel('Around-boresight accuracy (arcseconds)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.xticks(resolutions, [f'{r}x{r}' for r in resolutions])
    
    out_path_around = os.path.join(project_root, 'ablation_resolution_around_boresight.png')
    plt.savefig(out_path_around, bbox_inches='tight', dpi=150)
    print(f"Saved around-boresight accuracy plot to: {out_path_around}")
    
    # Print summary table
    print(f"\n{'='*80}")
    print(f"SUMMARY: Resolution Ablation Study (FOV={fov}°, {runs} runs per resolution)")
    print(f"{'='*80}")
    print(f"{'Resolution':<15} {'Success Rate':<15} {'Cross-Bore (\")':<20} {'Around-Bore (\")':<20}")
    print(f"{'-'*80}")
    for i, res_val in enumerate(resolutions):
        success = f"{success_rates[i]:.1%}" if success_rates[i] is not None else "N/A"
        cross = f"{cross_boresight_accuracies[i]:.2f}" if cross_boresight_accuracies[i] is not None else "N/A"
        around = f"{around_boresight_accuracies[i]:.2f}" if around_boresight_accuracies[i] is not None else "N/A"
        print(f"{res_val}x{res_val:<11} {success:<15} {cross:<20} {around:<20}")
    print(f"{'='*80}")
    
    return {
        'resolutions': resolutions,
        'success_rates': success_rates,
        'cross_boresight_accuracies': cross_boresight_accuracies,
        'around_boresight_accuracies': around_boresight_accuracies
    }

def test_ablation_study_fovs(runs_per_fov: str | int = 50):
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
    cross_boresight_accuracies = []
    around_boresight_accuracies = []

    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    for fov in fovs:
        successes = 0
        centroid_errors = []
        detected_stars_list = []
        for i in range(runs):
            attitude_fname = f"attitude_fov{fov}_run{i}.txt"
            centroid_fname = f"input_centroids_fov{fov}_run{i}.txt"
            actual_centroid_fname = f"actual_centroids_fov{fov}_run{i}.txt"
            compare_star_ids_fname = f"compare_star_ids_fov{fov}_run{i}.txt"
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
                print_input_centroids=centroid_fname,
                print_actual_centroids=actual_centroid_fname,
                compare_star_ids=compare_star_ids_fname,
                fov_deg=float(fov)
            )

            ok = run_entire_pipeline(options)
            if not ok:
                print(f"Pipeline run failed for FOV={fov} (run {i})")
                continue

            res = parse_attitude_result(attitude_fname)
            if res is None:
                print(f"Could not parse attitude result for FOV={fov} (run {i})")
            elif res.known:
                successes += 1

            # Compute centroid error for this run
            lost_dir = os.path.join(project_root, "lost")
            input_path = os.path.join(lost_dir, centroid_fname)
            actual_path = os.path.join(lost_dir, actual_centroid_fname)
            avg_error = None
            try:
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
                if os.path.isfile(input_path) and os.path.isfile(actual_path):
                    with open(input_path, "r") as fi, open(actual_path, "r") as fa:
                        inp_pts = parse_points(fi)
                        act_pts = parse_points(fa)
                    n = min(len(inp_pts), len(act_pts))
                    if n > 0:
                        import math
                        errs = [math.hypot(inp_pts[j][0] - act_pts[j][0], inp_pts[j][1] - act_pts[j][1]) for j in range(n)]
                        avg_error = sum(errs) / n
                if avg_error is not None:
                    centroid_errors.append(avg_error)
                    print("average centroiding error", avg_error)
            except Exception as e:
                print(f"Error computing centroid accuracy for FOV={fov} run {i}: {e}")

            # Count detected stars for this run
            actual_centroids_path = actual_path
            try:
                if os.path.isfile(actual_centroids_path):
                    with open(actual_centroids_path, "r") as f:
                        import re
                        num_re = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
                        count = 0
                        for line in f:
                            if num_re.search(line):
                                count += 1
                    detected_stars_list.append(count)
                    print("detected stars: ", count)
            except Exception as e:
                print(f"Error counting detected stars for FOV={fov} run {i}: {e}")

            # Clean up files
            # for fname in [attitude_fname, centroid_fname, actual_centroid_fname, compare_star_ids_fname]:
            #     try:
            #         file_path = os.path.join(lost_dir, fname)
            #         if os.path.exists(file_path):
            #             os.remove(file_path)
            #     except Exception:
            #         pass

        rate = successes / runs if runs > 0 else 0.0
        success_rates.append(rate)
        print(f"FOV={fov} deg: success {successes}/{runs} -> {rate:.2%}")

        # Compute average centroid error and detected stars for this FOV
        avg_centroid_error = sum(centroid_errors) / len(centroid_errors) if centroid_errors else None
        avg_detected_stars = sum(detected_stars_list) / len(detected_stars_list) if detected_stars_list else None
        # Compute cross-boresight and around-boresight accuracy
        cross_acc = None
        around_acc = None
        if avg_centroid_error is not None and avg_detected_stars is not None and avg_detected_stars > 0:
            from lost_api import compute_cross_boresight_accuracy, compute_around_boresight_accuracy
            try:
                cross_acc = compute_cross_boresight_accuracy(
                    fov_deg=float(fov),
                    avg_centroid_accuracy=avg_centroid_error,
                    num_pixels=1024,
                    avg_detected_stars=avg_detected_stars,
                )
                # Convert cross-boresight accuracy from degrees to arcseconds
                cross_acc_arcsec = cross_acc * 3600.0 if cross_acc is not None else None
                around_acc = compute_around_boresight_accuracy(
                    avg_centroid_accuracy=avg_centroid_error,
                    num_pixels=1024,
                    avg_detected_stars=avg_detected_stars,
                )
                # Convert around-boresight accuracy from degrees to arcseconds
                around_acc_arcsec = around_acc * 3600.0 if around_acc is not None else None
            except Exception as e:
                print(f"Error computing accuracy metrics for FOV={fov}: {e}")
                cross_acc_arcsec = None
                around_acc_arcsec = None
        else:
            cross_acc_arcsec = None
            around_acc_arcsec = None
        cross_boresight_accuracies.append(cross_acc_arcsec)
        around_boresight_accuracies.append(around_acc_arcsec)

    # Plot results
    plt.figure(figsize=(8, 5))
    plt.plot(fovs, success_rates, marker='o', label='Success Rate')
    plt.title('Ablation study: attitude success rate vs FOV')
    plt.xlabel('FOV (degrees)')
    plt.ylabel('Success rate')
    plt.ylim(-0.05, 1.05)
    plt.grid(True)
    plt.legend()

    out_path = os.path.join(project_root, 'ablation_fov.png')
    plt.savefig(out_path, bbox_inches='tight')
    print(f"Saved ablation plot to: {out_path}")

    # Plot cross-boresight accuracy (arcseconds)
    plt.figure(figsize=(8, 5))
    plt.plot(fovs, cross_boresight_accuracies, marker='o', label='Cross-Boresight Accuracy (arcsec)')
    plt.title('Ablation study: cross-boresight accuracy vs FOV')
    plt.xlabel('FOV (degrees)')
    plt.ylabel('Cross-boresight accuracy (arcseconds)')
    plt.grid(True)
    plt.legend()
    out_path_cross = os.path.join(project_root, 'ablation_cross_boresight.png')
    plt.savefig(out_path_cross, bbox_inches='tight')
    print(f"Saved cross-boresight accuracy plot to: {out_path_cross}")

    # Plot around-boresight accuracy (arcseconds)
    plt.figure(figsize=(8, 5))
    plt.plot(fovs, around_boresight_accuracies, marker='o', label='Around-Boresight Accuracy (arcsec)')
    plt.title('Ablation study: around-boresight (roll) accuracy vs FOV')
    plt.xlabel('FOV (degrees)')
    plt.ylabel('Around-boresight accuracy (arcseconds)')
    plt.grid(True)
    plt.legend()
    out_path_around = os.path.join(project_root, 'ablation_around_boresight.png')
    plt.savefig(out_path_around, bbox_inches='tight')
    print(f"Saved around-boresight accuracy plot to: {out_path_around}")


def test_average_accuracy_across_trials(
    n_trials: int = 100,
    ra: float | None = None,
    de: float | None = None,
    roll: float | None = None,
    x_res: int = 1024,
    y_res: int = 1024,
    fov_deg: float = 25.0
):
    """Compute average centroiding accuracy and detected stars across multiple trials,
    then calculate cross-boresight and around-boresight accuracy metrics.
    
    Args:
        n_trials: Number of trials to run
        ra: Right ascension (if None, uses random attitudes)
        de: Declination (if None, uses random attitudes)
        roll: Roll angle (if None, uses random attitudes)
        x_res: Image width in pixels
        y_res: Image height in pixels
        fov_deg: Field of view in degrees
        
    Returns:
        dict with keys:
            - avg_centroid_error: average centroiding error in pixels
            - avg_detected_stars: average number of detected stars
            - cross_boresight_arcsec: cross-boresight accuracy in arcseconds
            - around_boresight_arcsec: around-boresight accuracy in arcseconds
            - around_boresight_deg: around-boresight accuracy in degrees
            - around_boresight_rad: around-boresight accuracy in radians
    """
    import math
    from lost_api import compute_cross_boresight_accuracy, compute_around_boresight_accuracy
    
    # Determine if using random or given attitude
    use_random = ra is None or de is None or roll is None
    
    centroid_errors = []
    detected_stars_list = []
    successful_trials = 0
    
    print(f"Running {n_trials} trials with {'random' if use_random else 'given'} attitude...")
    
    for i in range(n_trials):
        print(f"\nTrial {i+1}/{n_trials}")
        
        # Run centroiding test
        if use_random:
            centroid_result = test_run_centroiding_random_attitude()
        else:
            centroid_result = test_run_centroiding_given_attitude(ra, de, roll)
        
        # Run star ID test to get detected stars count
        if use_random:
            starid_result = test_run_starID_random_attitude()
        else:
            starid_result = test_run_starID_given_attitude(ra, de, roll)
        
        # Check if we got valid results
        if (centroid_result and 
            centroid_result.centroids_mean_error is not None and
            starid_result and 
            starid_result.starid_num_total is not None and
            starid_result.starid_num_total > 0):
            
            centroid_errors.append(centroid_result.centroids_mean_error)
            detected_stars_list.append(starid_result.starid_num_total)
            successful_trials += 1
            
            print(f"  Centroid error: {centroid_result.centroids_mean_error:.4f} pixels")
            print(f"  Detected stars: {starid_result.starid_num_total}")
        else:
            print(f"  Trial {i+1} failed - no valid results")
    
    if successful_trials == 0:
        print("\nNo successful trials - cannot compute accuracy metrics")
        return None
    
    # Compute averages
    avg_centroid_error = sum(centroid_errors) / len(centroid_errors)
    avg_detected_stars = sum(detected_stars_list) / len(detected_stars_list)
    
    print(f"\n{'='*60}")
    print(f"Results from {successful_trials}/{n_trials} successful trials:")
    print(f"  Average centroid error: {avg_centroid_error:.4f} pixels")
    print(f"  Average detected stars: {avg_detected_stars:.2f}")
    
    # Compute accuracy metrics
    # num_pixels = x_res
    
    try:
        cross_boresight = compute_cross_boresight_accuracy(
            fov_deg=fov_deg,
            avg_centroid_accuracy=avg_centroid_error,
            num_pixels=x_res,
            avg_detected_stars=avg_detected_stars
        )
        cross_boresight_arcsec = cross_boresight * 3600.0  # degrees to arcseconds
        
        around_boresight_rad = compute_around_boresight_accuracy(
            avg_centroid_accuracy=avg_centroid_error,
            num_pixels=x_res,
            avg_detected_stars=avg_detected_stars
        )
        around_boresight_deg = around_boresight_rad * 180.0 / math.pi
        around_boresight_arcsec = around_boresight_rad * 206265.0  # radians to arcseconds
        
        print(f"\nAccuracy Metrics:")
        print(f"  Cross-boresight accuracy: {cross_boresight_arcsec:.2f} arcseconds")
        print(f"  Around-boresight (roll) accuracy:")
        print(f"    {around_boresight_rad:.6f} radians")
        print(f"    {around_boresight_deg:.6f} degrees")
        print(f"    {around_boresight_arcsec:.2f} arcseconds")
        print(f"{'='*60}")
        
        return {
            'avg_centroid_error': avg_centroid_error,
            'avg_detected_stars': avg_detected_stars,
            'cross_boresight_arcsec': cross_boresight_arcsec,
            'around_boresight_rad': around_boresight_rad,
            'around_boresight_deg': around_boresight_deg,
            'around_boresight_arcsec': around_boresight_arcsec,
            'successful_trials': successful_trials,
            'total_trials': n_trials
        }
        
    except Exception as e:
        print(f"\nError computing accuracy metrics: {e}")
        return None

# def test_cross_boresight_accuracy(fov_deg: float = 20.0, x_res: int = 1024, y_res: int = 1024):
#     """Compute cross-boresight accuracy using centroid and star-id tests.

#     This function calls `test_centroiding_step` and `test_star_id_step` to
#     obtain the average centroiding error (in pixels) and the number of
#     detected stars, then calls `compute_cross_boresight_accuracy` from
#     `lost_api` to compute the metric. Returns the computed metric or None
#     if it cannot be computed.
#     """
#     # Run centroiding step and get average centroid error (pixels)
#     print("Running centroiding step to compute average centroid accuracy...")
#     avg_centroid_error = test_centroiding_step()
#     if avg_centroid_error is None:
#         print("Could not obtain average centroid accuracy; aborting cross-boresight computation")
#         return None

#     # Run star-id step and get detected star count
#     print("Running star-id step to count detected stars...")
#     detected_stars = test_star_id_step()
#     if detected_stars is None or detected_stars <= 0:
#         print("Could not obtain detected star count; aborting cross-boresight computation")
#         return None

#     from lost_api import compute_cross_boresight_accuracy

#     num_pixels = int(x_res) * int(y_res)
#     metric = compute_cross_boresight_accuracy(
#         fov_deg=fov_deg,
#         avg_centroid_accuracy=avg_centroid_error,
#         num_pixels=num_pixels,
#         avg_detected_stars=float(detected_stars),
#     )

#     print(f"Cross-boresight accuracy: {metric}")
#     return metric


# def test_around_boresight_accuracy(x_res: int = 1024, y_res: int = 1024):
#     """Compute around-boresight (roll) accuracy using centroid and star-id tests.

#     This function calls `test_centroiding_step` and `test_star_id_step` to
#     obtain the average centroiding error (in pixels) and the number of
#     detected stars, then calls `compute_around_boresight_accuracy` from
#     `lost_api` to compute the roll error metric. Returns the computed metric
#     (in radians) or None if it cannot be computed.
#     """
#     # Run centroiding step and get average centroid error (pixels)
#     print("Running centroiding step to compute average centroid accuracy...")
#     avg_centroid_error = test_centroiding_step()
#     if avg_centroid_error is None:
#         print("Could not obtain average centroid accuracy; aborting around-boresight computation")
#         return None

#     # Run star-id step and get detected star count
#     print("Running star-id step to count detected stars...")
#     detected_stars = test_star_id_step()
#     if detected_stars is None or detected_stars <= 0:
#         print("Could not obtain detected star count; aborting around-boresight computation")
#         return None

#     from lost_api import compute_around_boresight_accuracy
#     import math

#     num_pixels = int(x_res) * int(y_res)
#     metric_rad = compute_around_boresight_accuracy(
#         avg_centroid_accuracy=avg_centroid_error,
#         num_pixels=num_pixels,
#         avg_detected_stars=float(detected_stars),
#     )

#     # Convert to degrees and arcseconds for convenience
#     metric_deg = metric_rad * 180.0 / math.pi
#     metric_arcsec = metric_rad * 206265.0

#     print(f"Around-boresight (roll) accuracy:")
#     print(f"  {metric_rad:.6f} radians")
#     print(f"  {metric_deg:.6f} degrees")
#     print(f"  {metric_arcsec:.2f} arcseconds")
#     return metric_rad


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
        fov=25,
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

def test_run_centroiding_random_attitude() -> CentroidResult:
    """Run the pipeline with random attitude and return centroid results."""
    from lost_api import run_entire_pipeline_C, parse_attitude_result, PipelineOptions

    attitude_output_fname = "attitude.txt"
    options = PipelineOptions(
        generate=1,
        generate_x_res = 1024,
        generate_y_res = 1024,
        generate_random_attitudes=True,
        generate_exposure=0.6,
        fov_deg = 25.0,
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude=attitude_output_fname,
        centroid_compare_threshold=2,
    )
    ret = run_entire_pipeline_C(options)
    print(ret)
    print(parse_attitude_result(attitude_output_fname))
    return ret


def test_run_centroiding_given_attitude(
    ra_input: str | float, de_input: str | float, roll_input: str | float
) -> CentroidResult:
    """Run the pipeline with given attitude and return centroid results."""
    try:
        ra = float(ra_input)
        de = float(de_input)
        roll = float(roll_input)
    except ValueError:
        raise ValueError("ra, de, and roll must be convertible to float")
    from lost_api import run_entire_pipeline_C, parse_attitude_result, PipelineOptions

    attitude_output_fname = "attitude.txt"
    options = PipelineOptions(
        generate=1,
        generate_ra=ra,
        generate_de=de,
        generate_roll=roll,
        generate_random_attitudes=False,
        generate_exposure=0.6,
        fov_deg = 25.0,
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude=attitude_output_fname,
        centroid_compare_threshold=2,
    )
    ret = run_entire_pipeline_C(options)
    print(ret)
    print(parse_attitude_result(attitude_output_fname))
    return ret


def test_run_starID_random_attitude() -> StarIDResult:
    """Run the pipeline with random attitude and return star ID results."""
    from lost_api import run_entire_pipeline_S, parse_attitude_result, PipelineOptions

    attitude_output_fname = "attitude.txt"
    options = PipelineOptions(
        generate=1,
        generate_x_res = 1024,
        generate_y_res = 1024,
        generate_random_attitudes=True,
        generate_exposure=0.6,
        fov_deg = 25.0,
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude=attitude_output_fname,
        centroid_compare_threshold=2,
    )
    ret = run_entire_pipeline_S(options)
    print(ret)
    print(parse_attitude_result(attitude_output_fname))
    return ret


def test_run_starID_given_attitude(
    ra_input: str | float, de_input: str | float, roll_input: str | float
) -> StarIDResult:
    """Run the pipeline with given attitude and return star ID results."""
    try:
        ra = float(ra_input)
        de = float(de_input)
        roll = float(roll_input)
    except ValueError:
        raise ValueError("ra, de, and roll must be convertible to float")
    from lost_api import run_entire_pipeline_S, parse_attitude_result, PipelineOptions

    attitude_output_fname = "attitude.txt"
    options = PipelineOptions(
        generate=1,
        generate_ra=ra,
        generate_de=de,
        generate_roll=roll,
        generate_random_attitudes=False,
        generate_exposure=0.6,
        fov_deg = 25.0,
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude=attitude_output_fname,
    )
    ret = run_entire_pipeline_S(options)
    print(ret)
    print(parse_attitude_result(attitude_output_fname))
    return ret


def test_test():
    """Test both star ID and centroiding with given attitude."""
    a = test_run_starID_given_attitude(200, 40, 300)
    b = test_run_centroiding_given_attitude(200, 40, 300)
    print(a)
    print(b)


def evaluate_centroid(ax, params, runner, special_paramss, x_vals):
    """Evaluate centroid algorithms and plot on given Axes.

    - ax: matplotlib Axes to plot on
    - params: module/object providing centroid configuration (expects attributes used below)
    - runner: object with a `run_lost(cmd_args_list)` function that returns a dict with
      a 'centroids_mean_error' key
    - special_paramss: iterable of special-parameter lists to append to the base command
    - x_vals: x coordinates for plotting (one per element of special_paramss)
    """
    import math

    for algo_name, algo_params in params.centroid_algos:
        y_vals = []
        for special_params in special_paramss:
            cmd = ['--generate', str(params.centroid_num_trials),
                   '--generate-random-attitudes=true',
                   '--compare-centroids=-']
            # Allow params to provide base args and algo-specific args as lists
            if hasattr(params, 'centroid_base_args') and params.centroid_base_args:
                cmd += list(params.centroid_base_args)
            if algo_params:
                cmd += list(algo_params)
            cmd += list(special_params)

            ran = runner.run_lost(cmd)
            if isinstance(ran, dict) and 'centroids_mean_error' in ran:
                y_vals.append(ran['centroids_mean_error'])
            else:
                print('WARNING! No centroids detected or unexpected runner output!')
                y_vals.append(math.nan)
        ax.plot(x_vals, y_vals, label=algo_name, marker='.')

    ax.set_ylabel('Centroid Error (pixels, average)')
    ax.legend()


# def test_centroid(special_paramss, x_vals, params_module: str = 'params', runner_module: str = 'runner', out_fname: str = 'centroid_eval.png'):
#     """Helper to run evaluate_centroid using modules by name and save a plot.

#     special_paramss: list of lists of extra args to append for each x value
#     x_vals: list of x coordinates for plotting
#     params_module / runner_module: importable module names providing `params` and `runner`
#     """
#     try:
#         params = __import__(params_module)
#     except Exception as e:
#         print(f"Could not import params module '{params_module}': {e}")
#         return
#     try:
#         runner = __import__(runner_module)
#     except Exception as e:
#         print(f"Could not import runner module '{runner_module}': {e}")
#         return

#     import matplotlib.pyplot as plt

#     fig, ax = plt.subplots(figsize=(8, 5))
#     evaluate_centroid(ax, params, runner, special_paramss, x_vals)
#     fig.savefig(out_fname, bbox_inches='tight')
#     print(f"Saved centroid evaluation to: {out_fname}")


def _run_single_pipeline(run_id, options=None):
    """Helper function to run a single pipeline iteration with random attitude."""
    import os

    from lost_api import parse_attitude_result, PipelineOptions, run_entire_pipeline

    # Use unique output file for each parallel run
    # File will be deleted after run is completed
    attitude_output_fname = f"attitude_{os.getpid()}_{run_id}.txt"
    if options is None:
        options = PipelineOptions(
            generate=1,
            generate_random_attitudes=True,
            generate_exposure=1,
            fov=25,
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


def test_run_pipeline_n_times(n_input: str | int, options = None):
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
    if options is None:
        options = PipelineOptions(
            generate=1,
            generate_random_attitudes=True,
            generate_exposure=1,
            fov=25,
            centroid_algo="cog",
            database="tetra-20.dat",
            false_stars=0,
            star_id_algo="tetra",
            attitude_algo="dqm",
            centroid_mag_filter=0,
            print_attitude=attitude_output_fname,
        )
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        futures = {executor.submit(_run_single_pipeline(options), i): i for i in range(n)}

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
