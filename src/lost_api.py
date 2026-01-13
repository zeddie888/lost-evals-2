import os
import subprocess
import threading
import time
from dataclasses import dataclass

from tqdm import tqdm

"""
This file provides a Python API for interacting with LOST

Functions that are part of the public API should be wrappers or helpers
on top of LOST-produced outputs. Evaluation logic should go in evaluators.py
"""

# Public data structures and types


@dataclass
class CompileOptions:
    rebuild: bool = False  # will rebuild project with `make clean` first


@dataclass
class TetraDbOptions:
    # minimum magnitude for a catalog star to be used for a Tetra pattern,
    # remember that higher magnitude = dimmer star!
    min_mag: int = 7
    # maximum angle (in degrees) between stars in a Tetra star pattern
    tetra_max_angle: int = 20
    output: str = "tetra.dat"  # filename for db


@dataclass
class PipelineOptions:
    generate: int = 1
    generate_x_res: int = 1024
    generate_y_res: int = 1024
    generate_ra: float | None = None
    generate_de: float | None = None
    generate_roll: float | None = None
    # Exposure time (seconds)
    generate_exposure: float = 0.6
    generate_random_attitudes: bool = True
    fov: float = 20
    # Pipeline parameters
    centroid_algo: str = "cog"
    database: str = "tetra.dat"
    false_stars: int = 0
    star_id_algo: str = "tetra"
    # note: always use dqm, never quest
    attitude_algo: str = "dqm"
    # do not use any catalog stars with magnitude > this value,
    # i.e. dimmer than this value. 0 means no filtering
    centroid_mag_filter: int = 0
    # fname to output attitude estimates to
    print_attitude: str = "attitude.txt"
    # Field-of-view in degrees (optional). If provided, the pipeline will be
    # invoked with a --fov <deg> argument. If None, the pipeline default is used.
    fov_deg: float | None = None
    # Optional outputs for intermediate stages (centroiding, star-id)
    # If set, these will be passed directly to the LOST CLI and written
    # inside the `lost` directory.
    print_input_centroids: str | None = None
    print_actual_centroids: str | None = None
    compare_star_ids: str | None = None
    # centroid compare threshold
    centroid_compare_threshold: int = 2


@dataclass
class AttitudeResult:
    known: bool
    ra: float | None = None
    de: float | None = None
    roll: float | None = None


# Public API


def clean_lost() -> bool:
    try:
        lost_dir = _get_lost_dir()
        result = subprocess.run(
            ["make", "clean"],
            cwd=lost_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
        if result.returncode != 0:
            print(f"'make clean' failed with return code {result.returncode}")
            print(result.stderr.decode().strip())
            return False
        # NOTE: currently lost's make clean does not remove the binary
        return True
    except subprocess.CalledProcessError as e:
        print(f"Clean failed: {e.stderr.decode().strip()}")
        return False
    except Exception as e:
        print(f"Error: {e}")
        return False


def compile_lost(options: CompileOptions | None) -> bool:
    # For me, takes ~45s after rebuild, otherwise ~20s
    options = options or CompileOptions()
    lost_dir = _get_lost_dir()
    lost_binary_path = os.path.join(lost_dir, "lost")

    try:
        # Delete the old lost binary if it exists
        if os.path.isfile(lost_binary_path):
            os.remove(lost_binary_path)

        if options.rebuild and not clean_lost():
            return False

        with tqdm(
            total=1,
            desc="Compiling LOST...",
            bar_format="{desc} {bar} {elapsed}",
            ncols=60,
        ) as pbar:
            # Start a background thread to animate the bar
            stop_spinner = False

            def animate():
                while not stop_spinner:
                    pbar.update(0)  # Just to refresh the bar
                    time.sleep(0.1)

            spinner_thread = threading.Thread(target=animate)
            spinner_thread.start()

            result = subprocess.run(
                ["make", "-j4", "LOST_DISABLE_ASAN=1"],
                cwd=lost_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )

            stop_spinner = True
            spinner_thread.join()
            pbar.update(1)

        if result.returncode != 0:
            print(f"'make' failed with return code {result.returncode}")
            print(result.stderr.decode().strip())
            return False
        return os.path.isfile(lost_binary_path)
    except subprocess.CalledProcessError as e:
        print(f"Compilation failed: {e.stderr.decode().strip()}")
        return False
    except Exception as e:
        print(f"Error: {e}")
        return False


def generate_tetra_database(options: TetraDbOptions) -> bool:
    """Generate database for Tetra star ID algorithm.

    For sanity checking, note that with min_mag=7 and tetra_max_angle=20, the
    database size is about 41.6 MB and contains 2,591,420 patterns

    Detailed output:
        Narrowed catalog has 8773 stars.
        Tetra max angle: 20
        Tetra processed catalog has 3222 stars.
        Number of pattern stars: 1770
        Tetra found 2591420 patterns
        Generated database with 41609592 bytes
    """
    lost_dir = _get_lost_dir()
    output_path = os.path.join(lost_dir, options.output)
    cmd = [
        "./lost",
        "database",
        "--min-mag", str(options.min_mag),
        "--tetra",
        "--tetra-max-angle", str(options.tetra_max_angle),
        "--output",
        options.output,
    ]
    try:
        result = subprocess.run(cmd, cwd=lost_dir, check=True)
        if result.returncode != 0:
            print(f"Database generation failed with return code {result.returncode}")
            return False
        return os.path.isfile(output_path)
    except subprocess.CalledProcessError as e:
        print(f"Database generation failed: {e}")
        return False
    except Exception as e:
        print(f"Error: {e}")
        return False


def run_entire_pipeline(options: PipelineOptions) -> bool:
    lost_dir = _get_lost_dir()

    # Either use random attitude or all of ra, de, roll must be provided
    if not options.generate_random_attitudes and not all(
        v is not None
        for v in [options.generate_ra, options.generate_de, options.generate_roll]
    ):
        print(
            "Error: Must specify RA, DE, and ROLL when generate_random_attitudes is False."
        )
        return False

    cmd = [
        "./lost",
        "pipeline",
        "--generate", str(options.generate),
        "--generate-x-resolution", str(options.generate_x_res),
        "--generate-y-resolution", str(options.generate_y_res),
        "--generate-exposure", str(options.generate_exposure),
        "--fov", str(options.fov),
        "--centroid-algo", options.centroid_algo,
        "--database", options.database,
        "--false-stars", str(options.false_stars),
        "--star-id-algo", options.star_id_algo,
        "--attitude-algo", options.attitude_algo,
        "--centroid-mag-filter", str(options.centroid_mag_filter),
        "--print-attitude", options.print_attitude,
        "--plot-input", "input-foo-test.png",
        "--plot-raw-input", "raw-input-foo-test.png",
    ]

    if options.generate_random_attitudes:
        cmd += ["--generate-random-attitudes", "true"]
    else:
        cmd += [
            "--generate-ra", str(options.generate_ra),
            "--generate-de", str(options.generate_de),
            "--generate-roll", str(options.generate_roll),
        ]
    if options.fov_deg is not None:
        cmd += ["--fov", str(options.fov_deg)]

    # Optional intermediate outputs
    if getattr(options, "print_input_centroids", None):
        cmd += ["--print-input-centroids", options.print_input_centroids]
    if getattr(options, "print_actual_centroids", None):
        cmd += ["--print-actual-centroids", options.print_actual_centroids]
    if getattr(options, "compare_star_ids", None):
        cmd += ["--compare-star-ids", options.compare_star_ids]

    try:
        result = subprocess.run(
            cmd,
            cwd=lost_dir,
        )
        if result.returncode != 0:
            print(f"Pipeline failed with return code {result.returncode}")
            print(result.stderr.decode().strip())
            return False
        return True
    except subprocess.CalledProcessError as e:
        print(f"Pipeline failed: {e.stderr.decode().strip()}")
        return False
    except Exception as e:
        print(f"Error: {e}")
        return False


def parse_attitude_result(attitude_fname: str) -> AttitudeResult | None:
    """Parse the attitude output file generated by a full LOST pipeline run (i.e. run_entire_pipeline)."""
    try:
        lost_dir = _get_lost_dir()
        attitude_path = os.path.join(lost_dir, attitude_fname)
        with open(attitude_path, "r") as f:
            lines = [line.strip() for line in f if line.strip()]
        attitude = {}
        for line in lines:
            if " " in line:
                key, value = line.split(" ", 1)
                attitude[key] = value
        if attitude.get("attitude_known", "0") == "1":
            ra_str = attitude.get("attitude_ra")
            de_str = attitude.get("attitude_de")
            roll_str = attitude.get("attitude_roll")
            if ra_str is None or de_str is None or roll_str is None:
                return AttitudeResult(known=False)
            ra = float(ra_str)
            de = float(de_str)
            roll = float(roll_str)
            return AttitudeResult(known=True, ra=ra, de=de, roll=roll)
        else:
            return AttitudeResult(known=False)
    except Exception as e:
        print(f"Error parsing attitude result: {e}")
        return None


# Centroiding test result
@dataclass
class CentroidResult:
    centroids_num_correct: int | None = None
    centroids_num_extra: int | None = None
    centroids_mean_error: float | None = None

def run_entire_pipeline_C(options: PipelineOptions) -> CentroidResult | None:
    """Run the entire pipeline with centroid comparison output.
    
    Returns a CentroidResult with centroid statistics, or None on error.
    """
    lost_dir = _get_lost_dir()

    # Either use random attitude or all of ra, de, roll must be provided
    if not options.generate_random_attitudes and not all(
        v is not None
        for v in [options.generate_ra, options.generate_de, options.generate_roll]
    ):
        print(
            "Error: Must specify RA, DE, and ROLL when generate_random_attitudes is False."
        )
        return None

    cmd = [
        "./lost",
        "pipeline",
        "--generate", str(options.generate),
        "--generate-x-resolution", str(options.generate_x_res),
        "--generate-y-resolution", str(options.generate_y_res),
        "--generate-exposure", str(options.generate_exposure),
        "--fov", str(options.fov),  
        "--centroid-algo", options.centroid_algo,
        "--database", options.database,
        "--false-stars", str(options.false_stars),
        "--star-id-algo", options.star_id_algo,
        "--attitude-algo", options.attitude_algo,
        "--centroid-mag-filter", str(options.centroid_mag_filter),
        "--print-attitude", options.print_attitude,
        "--plot-input", "input-foo-zeddie-test.png",
        "--compare-centroids", "-",
        "--centroid-compare-threshold", str(options.centroid_compare_threshold),
    ]

    # Check if assigned attitude or go random
    if options.generate_random_attitudes:
        cmd += ["--generate-random-attitudes", "true"]
    else:
        cmd += [
            "--generate-ra", str(options.generate_ra),
            "--generate-de", str(options.generate_de),
            "--generate-roll", str(options.generate_roll),
        ]

    # ADD: Support for optional fov_deg override
    if options.fov_deg is not None:
        cmd += ["--fov", str(options.fov_deg)]

    try:
        result = subprocess.run(
            cmd, cwd=lost_dir, stdout=subprocess.PIPE, text=True
        )
        if result.returncode != 0:
            print(f"Pipeline failed with return code {result.returncode}")
            print(result.stderr)
            return None
    except subprocess.CalledProcessError as e:
        print(f"Pipeline failed: {e.stderr}")
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None

    output = result.stdout
    cent_correct = None
    cent_extra = None
    cent_error = None
    # get each line of outputs and grab the needed values
    for raw in output.splitlines():
        print(raw)
        # clean up each line's spaces, continue if empty
        line = raw.strip()
        if not line:
            continue
        if line.startswith("centroids_num_correct"):
            try:
                cent_correct = int(line.split()[-1])
            except Exception:
                print("Failure on centroids_num_correct")
        elif line.startswith("centroids_num_extra"):
            try:
                cent_extra = int(line.split()[-1])
            except Exception:
                print("Failure on centroids_num_extra")
        elif line.startswith("centroids_mean_error"):
            try:
                cent_error = float(line.split()[-1])
            except Exception:
                print("Failure on centroids_mean_error")
    return CentroidResult(
        centroids_num_correct=cent_correct,
        centroids_num_extra=cent_extra,
        centroids_mean_error=cent_error,
    )


#STAR ID
@dataclass
class StarIDResult:
    starid_num_correct: int | None = None
    starid_num_incorrect: int | None = None
    starid_num_total: int | None = None
    
def run_entire_pipeline_S(options: PipelineOptions) -> StarIDResult | None:
    """Run the entire pipeline with star ID comparison output.
    
    Returns a StarIDResult with star ID statistics, or None on error.
    """
    lost_dir = _get_lost_dir()

    # Either use random attitude or all of ra, de, roll must be provided
    if not options.generate_random_attitudes and not all(
        v is not None
        for v in [options.generate_ra, options.generate_de, options.generate_roll]
    ):
        print(
            "Error: Must specify RA, DE, and ROLL when generate_random_attitudes is False."
        )
        return None

    cmd = [
        "./lost",
        "pipeline",
        "--generate", str(options.generate),
        "--generate-x-resolution", str(options.generate_x_res),
        "--generate-y-resolution", str(options.generate_y_res),
        "--generate-exposure", str(options.generate_exposure),
        "--fov", str(options.fov),  
        "--centroid-algo", options.centroid_algo,
        "--database", options.database,
        "--false-stars", str(options.false_stars),
        "--star-id-algo", options.star_id_algo,
        "--attitude-algo", options.attitude_algo,
        "--centroid-mag-filter", str(options.centroid_mag_filter),
        "--print-attitude", options.print_attitude,
        "--plot-input", "input-foo-zeddie-test.png",
        "--compare-star-ids", "-",
    ]

    if options.generate_random_attitudes:
        cmd += ["--generate-random-attitudes", "true"]
    else:
        cmd += [
            "--generate-ra", str(options.generate_ra),
            "--generate-de", str(options.generate_de),
            "--generate-roll", str(options.generate_roll),
        ]

    # ADD: Support for optional fov_deg override
    if options.fov_deg is not None:
        cmd += ["--fov", str(options.fov_deg)]

    try:
        result = subprocess.run(
            cmd, cwd=lost_dir, stdout=subprocess.PIPE, text=True
        )
        if result.returncode != 0:
            print(f"Pipeline failed with return code {result.returncode}")
            print(result.stderr)
            return None
    except subprocess.CalledProcessError as e:
        print(f"Pipeline failed: {e.stderr}")
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None

    output = result.stdout
    sid_correct = None
    sid_incorrect = None
    sid_total = None
    # get each line of outputs and grab the needed values
    for raw in output.splitlines():
        print(raw)
        # clean up each line's spaces, continue if empty
        line = raw.strip()
        if not line:
            continue
        if line.startswith("starid_num_correct"):
            try:
                sid_correct = int(line.split()[-1])
            except Exception:
                print("Failure on starid_num_correct")
        elif line.startswith("starid_num_incorrect"):
            try:
                sid_incorrect = int(line.split()[-1])
            except Exception:
                print("Failure on starid_num_incorrect")
        elif line.startswith("starid_num_total"):
            try:
                sid_total = int(line.split()[-1])
            except Exception:
                print("Failure on starid_num_total")
    return StarIDResult(
        starid_num_correct=sid_correct,
        starid_num_incorrect=sid_incorrect,
        starid_num_total=sid_total,
    )

def compute_cross_boresight_accuracy(
    fov_deg: float,
    avg_centroid_accuracy: float,
    num_pixels: int,
    avg_detected_stars: float,
) -> float:
    """Compute the cross-boresight accuracy metric.

    Formula:
        E_cross_boresight = (A * E_centroid) / (N_pixel * sqrt(N_star))

    Where:
        A = FOV (field-of-view in degrees)
        E_centroid = average centroiding accuracy (pixels)
        N_pixel = # of pixels across the focal plane
        N_star = average number of detected stars

    Args:
        fov_deg: field-of-view in degrees (A).
        avg_centroid_accuracy: average centroiding accuracy in pixels (E_centroid).
        num_pixels: # of pixels across the focal plane
        avg_detected_stars: average number of detected stars (N_star, can be fractional if averaged).

    Returns:
        The cross-boresight accuracy as a float.

    Raises:
        ValueError: if any input is non-positive where a positive value is required.

    Notes:
        - The function assumes numeric inputs. It will raise on zero/negative values for
          `num_pixels` or `avg_detected_stars` because those would make the denominator
          zero or undefined.
    """
    # Validate inputs
    if fov_deg is None:
        raise ValueError("fov_deg must be provided and positive")
    if avg_centroid_accuracy is None:
        raise ValueError("avg_centroid_accuracy must be provided")
    try:
        fov = float(fov_deg)
        avg_cent = float(avg_centroid_accuracy)
        pixels = int(num_pixels)
        avg_stars = float(avg_detected_stars)
    except Exception as e:
        raise ValueError(f"Invalid input types: {e}")

    if pixels <= 0:
        raise ValueError("num_pixels must be > 0")
    if avg_stars <= 0:
        raise ValueError("avg_detected_stars must be > 0")

    # Compute denominator: num_pixels * sqrt(avg_detected_stars)
    import math

    denom = pixels * math.sqrt(avg_stars)
    result = fov * (avg_cent) / denom
    return float(result)


def compute_around_boresight_accuracy(
    avg_centroid_accuracy: float,
    num_pixels: int,
    avg_detected_stars: float,
) -> float:
    """Compute the around-boresight (roll) accuracy metric.

    Formula:
        E_roll = atan( E_centroid / (0.3825 * N_pixel) ) * (1 / sqrt(N_star))

    Where:
        E_centroid: average centroiding accuracy (pixels)
        N_pixel: total number of pixels in the image (width * height)
        N_star: average number of detected stars

    Args:
        avg_centroid_accuracy: average centroiding accuracy in pixels.
        num_pixels: total number of pixels in the image (width * height).
        avg_detected_stars: average number of detected stars (can be fractional if averaged).

    Returns:
        The around-boresight roll error.

    Raises:
        ValueError: if any input is invalid or non-positive where required.

    Notes:
        - The result is in radians. To convert to degrees: result * 180 / math.pi
        - To convert to arcseconds: result * 206265 (radians to arcseconds)
        - The constant 0.3825 is empirically determined for the camera/optics system.
    """
    import math

    # Validate inputs
    if avg_centroid_accuracy is None:
        raise ValueError("avg_centroid_accuracy must be provided")
    try:
        avg_cent = float(avg_centroid_accuracy)
        pixels = int(num_pixels)
        avg_stars = float(avg_detected_stars)
    except Exception as e:
        raise ValueError(f"Invalid input types: {e}")

    if pixels <= 0:
        raise ValueError("num_pixels must be > 0")
    if avg_stars <= 0:
        raise ValueError("avg_detected_stars must be > 0")

    # Compute: atan( E_centroid / (0.3825 * N_pixel) ) * (1 / sqrt(N_star))
    numerator = avg_cent / (0.3825 * pixels)
    atan_term = math.atan(numerator)
    result = atan_term * (1.0 / math.sqrt(avg_stars))

    return float(result)

# Internal helpers, not meant for use outside this module


def _get_lost_dir() -> str:
    """Get the path to the symlinked lost directory."""
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    lost_dir = os.path.join(project_root, "lost")
    if not os.path.isdir(lost_dir):
        raise FileNotFoundError(
            f"LOST directory not found at expected location: {lost_dir}"
        )
    return lost_dir
