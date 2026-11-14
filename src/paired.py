"""Run a serial, paired LOST star-ID pilot (Python standard library only)."""

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import random
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path

ALGORITHMS = ("tetra", "py")
LOST_COMMIT = "a2303e6080bcbb2763596c26620f267d6c9bcd11"
DATABASE_OPTIONS = {
    "tetra": ["--tetra", "--tetra-max-angle", "20"],
    "py": ["--kvector", "--kvector-min-distance", "0.2",
           "--kvector-max-distance", "20", "--kvector-distance-bins", "10000"],
}


def score(correct: int, wrong: int) -> str:
    if wrong:
        return "wrong"
    return "correct" if correct >= 4 else "no_solve"


def input_counts(data: dict[str, str]) -> tuple[int, int]:
    count = int(data["num_input_centroids"])
    true = sum(f"input_centroid_{i}_id" in data for i in range(count))
    return true, count - true


def read_fields(path: Path) -> dict[str, str]:
    return dict(line.split(maxsplit=1) for line in path.read_text().splitlines())


def run(command, cwd, log, *, memory=False):
    with log.open("w") as stream:
        if memory:
            # A fresh Python parent measures only its one LOST child.
            command = [sys.executable, "-c", MEMORY_CHILD, *map(str, command)]
        try:
            result = subprocess.run(command, cwd=cwd, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=120,
                                    env={**os.environ, "LOST_BSC_PATH":
                                         str(cwd / "bright-star-catalog.tsv")})
        except subprocess.TimeoutExpired:
            return "timeout"
    return {0: "ok", 124: "timeout"}.get(
        result.returncode, "crash" if result.returncode < 0 else "failed"
    )


MEMORY_CHILD = """
import os, resource, subprocess, sys
try:
    result = subprocess.run(sys.argv[1:], timeout=110)
except subprocess.TimeoutExpired:
    sys.exit(124)
if result.returncode < 0:
    os.kill(os.getpid(), -result.returncode)
peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
print("peak_rss_bytes", peak if sys.platform == "darwin" else peak * 1024)
sys.exit(result.returncode)
"""


def pipeline(binary, case):
    ra, de, roll = case["attitude"]
    return [binary, "pipeline", "--generate", "1", "--generate-centroids-only",
            "--fov", "12", "--generate-time-based-seed", "false",
            "--generate-seed", str(case["seed"]), "--generate-ra", str(ra),
            "--generate-de", str(de), "--generate-roll", str(roll),
            "--generate-cutoff-mag", "-5" if case["condition"] == "negative" else "100",
            "--generate-false-stars", str(case["injected"])]


def read_input(path):
    data = read_fields(path)
    true, false = input_counts(data)
    for i in range(true + false):
        for axis in ("x", "y"):
            if not math.isfinite(float(data[f"input_centroid_{i}_{axis}"])):
                raise ValueError("non-finite centroid")
    return data


def generate_cases(binary, cwd, out, count):
    rng = random.Random(394859)
    cases = [{"id": "smoke", "condition": "smoke", "seed": 1, "injected": 0,
              "attitude": [218, 12, 302]}]
    probe = out / "probe.txt"
    log = out / "probe.log"
    for i in range(count + 10):
        attitude = [rng.uniform(0, 360), math.degrees(math.asin(rng.uniform(-1, 1))),
                    rng.uniform(0, 360)]
        negative = i >= count
        if not negative:
            cases.append({"id": f"clean-{i}", "condition": "clean", "seed": i + 1,
                          "injected": 0, "attitude": attitude})
        case = {"id": f"{'negative' if negative else 'false'}-{i}",
                "condition": "negative" if negative else "false",
                "injected": 1200, "attitude": attitude}
        for _ in range(150):
            case["seed"] = rng.randrange(1, 2**31)
            if run(pipeline(binary, case) + ["--print-input-centroids", probe],
                   cwd, log) != "ok":
                raise RuntimeError(f"Case generation failed; see {log}")
            true, false = input_counts(read_input(probe))
            if false == 5 and (not negative or true == 0):
                cases.append(case)
                break
        else:
            raise RuntimeError(f"No five-false-centroid case found for {case['id']}")
    return cases


def evaluate(binary, cwd, out, case, algorithm, database):
    prefix = out / f"{case['id']}-{algorithm}"
    inputs, comparison, speed = (Path(f"{prefix}.{suffix}") for suffix in
                                 ("input", "ids", "speed"))
    cmd = pipeline(binary, case) + ["--database", database, "--star-id-algo", algorithm,
                                    "--print-input-centroids", inputs]
    log = Path(f"{prefix}.log")
    status = run(cmd + ["--centroid-compare-threshold", "0",
                        "--compare-star-ids", comparison, "--print-speed", speed], cwd, log)
    row = dict(case=case["id"], condition=case["condition"], algorithm=algorithm,
               status=status, outcome="", correct_ids="", wrong_ids="", true_centroids="",
               false_centroids="", starid_ns="", peak_rss_bytes="", memory_status="")
    if status != "ok":
        return row, None
    try:
        data = read_input(inputs)
        true, false = input_counts(data)
        ids = read_fields(comparison)
        correct = int(ids["starid_num_correct"])
        wrong = int(ids["starid_num_incorrect"])
        timing = int(read_fields(speed).get("starid_average_ns", "0"))
        if (int(ids["starid_num_total"]) != true or timing < 0
                or not (0 <= correct <= true and 0 <= wrong <= true + false - correct)):
            raise ValueError("invalid match counts")
    except (OSError, KeyError, ValueError) as exc:
        with log.open("a") as stream:
            print(f"Invalid output: {exc}", file=stream)
        row["status"] = "invalid_output"
        return row, None
    row.update(outcome=score(correct, wrong), correct_ids=correct, wrong_ids=wrong,
               true_centroids=true, false_centroids=false, starid_ns=timing or "")
    return row, data


def profile_memory(binary, cwd, out, case, row, database):
    log = out / f"{case['id']}-{row['algorithm']}.memory.log"
    status = run(pipeline(binary, case) + ["--database", database,
                 "--star-id-algo", row["algorithm"]], cwd, log, memory=True)
    if status == "ok":
        values = [line.split()[-1] for line in log.read_text().splitlines()
                  if line.startswith("peak_rss_bytes ")]
        if len(values) == 1 and values[0].isdigit() and int(values[0]) > 0:
            row["peak_rss_bytes"] = int(values[0])
        else:
            status = "invalid_output"
    row["memory_status"] = status


def readout(rows, costs):
    lines = ["# Tetra vs. Pyramid pilot", "",
             "Paired star-ID comparison at 12° horizontal FOV. False and failed solves "
             "are shown separately. Unverified pairs are excluded from comparison; RSS includes the database and input setup.", "",
             "| Algorithm | Condition | Cases | Correct | Wrong | No solve | Unverified | "
             "Median / p95 µs (n) | Max RSS MiB (n) |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for algorithm in ALGORITHMS:
        for condition in ("clean", "false", "negative"):
            group = [r for r in rows if r["algorithm"] == algorithm and r["condition"] == condition]
            verified = [r for r in group if r["paired_input_verified"]]
            counts = [sum(r["outcome"] == label for r in verified)
                      for label in ("correct", "wrong", "no_solve")]
            elapsed = sorted(r["starid_ns"] for r in verified if r["starid_ns"] != "")
            memory = [r["peak_rss_bytes"] for r in verified if r["peak_rss_bytes"] != ""]
            timing = (f"{statistics.median(elapsed)/1000:.1f} / "
                      f"{elapsed[math.ceil(.95 * len(elapsed)) - 1]/1000:.1f} ({len(elapsed)})"
                      if elapsed else "—")
            rss = f"{max(memory)/2**20:.1f} ({len(memory)})" if memory else "—"
            lines.append("| " + " | ".join(map(str, (algorithm, condition, len(group),
                         *counts, len(group)-sum(counts), timing, rss))) + " |")
    lines += ["", "On negative controls, any returned ID is a false accept. Ten controls "
              "cannot establish a rare false-accept rate.", "", "Database costs:"]
    for algorithm, cost in costs.items():
        lines.append(f"- {algorithm}: {cost['bytes']/2**20:.2f} MiB, "
                     f"build {cost['build_seconds']:.2f} s.")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lost-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--fields", type=int, default=100)
    args = parser.parse_args()
    if args.fields < 1 or not (args.lost_dir / "lost").is_file():
        parser.error("provide a built LOST directory and a positive field count")
    cwd = args.lost_dir.resolve()
    binary = cwd / "lost"
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    raw = output / "raw"
    raw.mkdir()
    if args.manifest:
        cases = json.loads(args.manifest.read_text())
    else:
        cases = generate_cases(binary, cwd, raw, args.fields)
    if (not isinstance(cases, list) or not cases
            or any(not re.fullmatch(r"[A-Za-z0-9_-]+", c["id"]) for c in cases)
            or len({c["id"] for c in cases}) != len(cases)
            or sum(c["condition"] == "smoke" for c in cases) != 1
            or any(c["condition"] not in ("smoke", "clean", "false", "negative") for c in cases)):
        parser.error("manifest must have unique safe case IDs and exactly one smoke case")
    (output / "cases.json").write_text(json.dumps(cases, indent=2) + "\n")
    costs = {}
    databases = {}
    for algorithm, options in DATABASE_OPTIONS.items():
        path = output / f"{algorithm}.dat"
        flags = ["--min-mag", "7", *options]
        started = time.perf_counter()
        log = raw / f"database-{algorithm}.log"
        if run([binary, "database", *flags, "--output", path], cwd, log) != "ok":
            raise RuntimeError(f"Database generation failed; see {log}")
        costs[algorithm] = dict(bytes=path.stat().st_size,
                                build_seconds=time.perf_counter() - started, flags=flags)
        databases[algorithm] = path
    revision = subprocess.run(["git", "-C", cwd, "rev-parse", "--show-prefix", "HEAD"],
                              capture_output=True, text=True)
    actual_commit = (revision.stdout.strip() if revision.returncode == 0
                     and revision.stdout.startswith("\n") else None)
    metadata = dict(platform=platform.platform(), machine=platform.machine(),
                    lost_target_commit=LOST_COMMIT, lost_checkout_commit=actual_commit,
                    lost_binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    catalog_sha256=hashlib.sha256((cwd / "bright-star-catalog.tsv").read_bytes()).hexdigest(),
                    databases=costs)
    (output / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    rows = []
    for case in cases:
        pair = [evaluate(binary, cwd, raw, case, algo, databases[algo]) for algo in ALGORITHMS]
        a, b = (data for _, data in pair)
        for row, _ in pair:
            row["paired_input_verified"] = a is not None and b is not None
        if a is not None and b is not None and a != b:
            raise RuntimeError(f"Unpaired centroid input: {case['id']}")
        for data in (a, b):
            if data is not None:
                true, false = input_counts(data)
                expected_false = 5 if case["condition"] in ("false", "negative") else 0
                if false != expected_false or (case["condition"] == "negative" and true):
                    raise RuntimeError(f"Unexpected centroid population: {case['id']}")
        if case["condition"] == "smoke" and any(r["outcome"] != "correct" for r, _ in pair):
            raise RuntimeError("Both algorithms must pass the smoke case; inspect raw logs")
        rows.extend(row for row, _ in pair)
    for index, case in enumerate(cases):
        start = index * len(ALGORITHMS)
        for row in rows[start:start + len(ALGORITHMS)]:
            profile_memory(binary, cwd, raw, case, row, databases[row["algorithm"]])
    with (output / "results.csv").open("w", newline="") as table:
        writer = csv.DictWriter(table, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "summary.md").write_text(readout(rows, costs))
    print(output / "summary.md")


if __name__ == "__main__":
    main()
