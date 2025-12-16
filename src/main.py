import sys
import inspect

from lost_api import CompileOptions


def test_compile_lost():
    from lost_api import compile_lost
    options = CompileOptions(
        rebuild=False
    )
    compile_lost(options)


def test_clean_lost():
    from lost_api import clean_lost
    clean_lost()


def test_generate_tetra_database():
    from lost_api import generate_tetra_database, TetraDbOptions
    options = TetraDbOptions(
        min_mag=7,
        tetra_max_angle=20,
        output="tetra-20.dat"
    )
    generate_tetra_database(options)


def test_run_entire_pipeline_random_attitude():
    from lost_api import run_entire_pipeline, parse_attitude_result, PipelineOptions
    attitude_output_fname = "attitude.txt"
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
        print_attitude=attitude_output_fname
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
    from lost_api import run_entire_pipeline, PipelineOptions
    options = PipelineOptions(
        generate=1,
        generate_ra=ra,
        generate_de=de,
        generate_roll=roll,
        generate_random_attitudes=False,
        generate_exposure=0.6,
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude="attitude.txt"
    )
    run_entire_pipeline(options)

from lost_api import CentroidResult
def test_run_centroiding_random_attitude() -> CentroidResult:
    from lost_api import run_entire_pipeline_C, parse_attitude_result, PipelineOptions
    attitude_output_fname = "attitude.txt"
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
        print_attitude=attitude_output_fname,
        centroid_compare_threshold = 2
    )
    ret = run_entire_pipeline_C(options)
    print(ret)
    print(parse_attitude_result(attitude_output_fname))
    return ret

def test_run_centroiding_given_attitude(ra_input: str | float, de_input: str | float, roll_input: str | float) -> CentroidResult:
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
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude=attitude_output_fname,
        centroid_compare_threshold=2
    )
    ret = run_entire_pipeline_C(options)
    print(ret)
    print(parse_attitude_result(attitude_output_fname))
    return ret

from lost_api import StarIDResult
def test_run_starID_random_attitude() -> StarIDResult:
    from lost_api import run_entire_pipeline_S, parse_attitude_result, PipelineOptions
    attitude_output_fname = "attitude.txt"
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
        print_attitude=attitude_output_fname,
        centroid_compare_threshold=2
    )
    ret = run_entire_pipeline_S(options)
    print(ret)
    print(parse_attitude_result(attitude_output_fname))
    return ret

def test_run_starID_given_attitude(ra_input: str | float, de_input: str | float, roll_input: str | float)-> StarIDResult:
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
        centroid_algo="cog",
        database="tetra-20.dat",
        false_stars=0,
        star_id_algo="tetra",
        attitude_algo="dqm",
        centroid_mag_filter=0,
        print_attitude=attitude_output_fname
    )
    ret = run_entire_pipeline_S(options)
    print(ret)
    print(parse_attitude_result(attitude_output_fname))
    return ret

def test_test():
    a = test_run_starID_given_attitude(200, 40, 300)
    b = test_run_centroiding_given_attitude(200, 40, 300)
    print(a)
    print(b)

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
