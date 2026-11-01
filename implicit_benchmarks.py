"""Bounded comparisons and compression diagnostics for implicit auxiliary designs."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys
import time
import datetime

import contraction as con
import implicit_design as im
import modular_lab as modular
import reconstruction_lab as lab
from run_benchmarks import bdd_count


def form_hash(design):
    canonical = [[[[v, [[Fraction(c).numerator, Fraction(c).denominator] for c in form[v]]]
                   for v in sorted(form)] for form in clause] for clause in design.P]
    return hashlib.sha256(json.dumps(canonical, separators=(",", ":")).encode()).hexdigest()


def inspect_design(design, n):
    q = design.M * design.degP // design.Q
    points = []
    for t in range(1, design.Q + 1):
        design.numeric_forms(t)
        point = design.at(t)
        ranks, sizes, distinct = [], [], []
        for group, rows in zip(point.groups, point.rows):
            basis = im.integer_row_basis(rows)
            rank = len(basis)
            ranks.append(rank)
            sizes.append(len(group))
            if rank == len(group):
                distinct.append(1 << len(group))
            elif rank == 0:
                distinct.append(1)
            elif len(group) <= 10:
                vectors = [[row.get(v, 0) for v in group] for row in basis]
                keys = {tuple(sum(c for i, c in enumerate(row) if mask >> i & 1) for row in vectors)
                        for mask in range(1 << len(group))}
                distinct.append(len(keys))
            else:
                distinct.append(None)
        active = len(lab.corrections(design, t, q)[1])
        points.append({"point": t, "group_sizes": sizes, "measurement_ranks": ranks,
                       "raw_measurement_rows": [len(rows) for rows in point.rows],
                       "distinct_measurements": distinct,
                       "unrestricted_tuples": math.prod(distinct) if all(v is not None for v in distinct) else None,
                       "unrestricted_tuple_lower_bound": 1 << sum(ranks),
                       "local_assignments": sum(1 << size for size in sizes),
                       "active_corrections": active,
                       "subset_queries": sum(math.comb(active, j) for j in range(min(q, active) + 1))})
    return {"nominal_auxiliary_coordinates": design.info["D"], "q": q,
            "declared_degree": design.degP,
            "actual_degree": max((len(cs) - 1 for clause in design.P for form in clause
                                   for cs in form.values()), default=0),
            "form_sha256": form_hash(design), "points": points,
            "operator_entries_queried": design.info.get("operator_entries_queried"),
            "operator_rows_built": design.info.get("operator_rows_built"),
            "labels_used": design.info["labels_used"], "copies": design.info["copies"]}


def worker(config):
    n, seed, Q, R = config["case"]
    cnf = lab.ref.random_bounded_occurrence_3cnf(n, n, random.Random(seed))
    method = config["method"]
    start = time.perf_counter()
    disabled_audits = 0
    if method.startswith("eager_no_span_audit"):
        # Research ablation: this function is used only by the reference
        # builder's exhaustive global row-space assertions, not its P values.
        original_solve = lab.ref._solve_in_row_space
        def omit_global_membership_audit(*_):
            nonlocal disabled_audits
            disabled_audits += 1
            return []
        lab.ref._solve_in_row_space = omit_global_membership_audit
        try:
            design = lab.ref.build_aux_design(cnf, 2, Q, R, random.Random(seed))
        finally:
            lab.ref._solve_in_row_space = original_solve
    elif method.startswith("eager"):
        design = lab.ref.build_aux_design(cnf, 2, Q, R, random.Random(seed))
    else:
        design = im.build_implicit_design(cnf, 2, Q, R, random.Random(seed),
            measurements="spectral" if method.startswith("spectral") else "required",
            schedules=config.get("schedules", "reference"),
            operator="sparse" if method.startswith("sparse") else "entries")
    build_seconds = time.perf_counter() - start
    answer, certificate, solve_seconds = None, None, None
    if config.get("solve", True):
        start = time.perf_counter()
        q = design.M * design.degP // Q
        reconstruction = (con.TaylorContraction(design, n) if method.endswith("contraction")
                          else lab.cofactor_reconstruction)
        answer, certificate = modular.exact_modular_count(
            design, n, q, reconstruction=reconstruction, use_disjoint_bound=True)
        solve_seconds = time.perf_counter() - start
        if answer != lab.ref.count_solutions(cnf):
            raise ArithmeticError("Full count disagrees with brute force")
    # Diagnostics/validation are untimed and are not counted as solving.
    return {"build_seconds": build_seconds, "solve_seconds": solve_seconds,
            "total_seconds": build_seconds + (solve_seconds or 0),
            "count": answer, "certificate": certificate, "design": inspect_design(design, n),
            "disabled_global_membership_audits": disabled_audits}


def run_child(config, timeout):
    try:
        result = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                                 "--worker", json.dumps(config)], capture_output=True,
                                text=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "timeout_seconds": timeout}
    if result.returncode:
        raise RuntimeError(result.stderr)
    return {"status": "complete", **json.loads(result.stdout)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker")
    parser.add_argument("--profile", choices=("comparison", "scan", "audit_ablation", "measurement"), default="comparison")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=30)
    parser.add_argument("--output", type=Path, default=Path("local-results/implicit.json"))
    args = parser.parse_args()
    if sys.flags.optimize or args.repeats < 1 or args.timeout <= 0:
        parser.error("Use normal Python, positive repeats, and a positive timeout")
    if args.worker:
        print(json.dumps(worker(json.loads(args.worker))))
        return
    report = {"format": "implicit-auxiliary-design-v1",
              "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "python": platform.python_version(), "platform": platform.platform(),
              "profile": args.profile, "repeats": args.repeats, "timeout_seconds": args.timeout,
              "timing": "Fresh child per trial, imports/input excluded; design and full requested solve included. Diagnostics excluded. No extrapolation of timeouts.",
              "diagnostic": "Unrestricted measurement tuples and correction subset counts are separate statistics; their product is not a runtime estimate. Rank alone does not determine the number of distinct measurements.",
              "memory": "No process-memory benchmark was performed.", "cases": []}
    if args.profile != "scan":
        if args.profile == "comparison":
            configurations = [(case, "reference") for case in
                ((6, 1, 2, 2), (6, 2, 3, 2), (8, 1, 3, 2), (8, 1, 4, 2), (8, 1, 4, 3))]
            methods = ["eager_measurement", "spectral_measurement", "required_measurement",
                       "eager_contraction", "required_contraction"]
        elif args.profile == "audit_ablation":
            configurations = [(case, "reference") for case in ((8, 1, 3, 2), (8, 1, 4, 2), (8, 1, 4, 3))]
            methods = ["eager_no_span_audit_contraction", "sparse_contraction", "required_contraction"]
        else:
            configurations = [((8, seed, 4, 2), "reference") for seed in (71, 72, 73)]
            configurations += [((8, seed, 8, 2), "balanced") for seed in (71, 72, 73)]
            methods = ["spectral_measurement", "required_measurement", "required_contraction"]
        for case, schedule in configurations:
            print(f"Comparing {case}...", flush=True)
            row = {"case": case, "schedules": schedule, "methods": {name: {"runs": []} for name in methods}}
            for repeat in range(args.repeats):
                for name in methods[repeat % len(methods):] + methods[:repeat % len(methods)]:
                    record = row["methods"][name]
                    if record.get("status") == "timeout":
                        continue
                    sample = run_child({"case": case, "method": name, "schedules": schedule}, args.timeout)
                    record["status"] = sample["status"]
                    if sample["status"] == "complete":
                        record["runs"].append({key: sample[key] for key in ("build_seconds", "solve_seconds", "total_seconds")})
                        record.update({key: sample[key] for key in
                                       ("count", "certificate", "design", "disabled_global_membership_audits")})
                    else:
                        record["timeout_seconds"] = args.timeout
            hashes = {value["design"]["form_sha256"] for value in row["methods"].values() if "design" in value}
            if len(hashes) != 1:
                raise ArithmeticError("Comparison methods used different replacement polynomials")
            for record in row["methods"].values():
                if record["status"] == "complete":
                    record["medians"] = {key: statistics.median(run[key] for run in record["runs"])
                                         for key in ("build_seconds", "solve_seconds", "total_seconds")}
            cnf = lab.ref.random_bounded_occurrence_3cnf(case[0], case[0], random.Random(case[1]))
            baselines = {"brute_force": lambda: (lab.ref.count_solutions(cnf), {}),
                         "ordinary_elimination": lambda: con.conventional_count(cnf),
                         "bdd": lambda: bdd_count(cnf)}
            import dd.autoref
            row["baselines"] = {}
            for name, method in baselines.items():
                samples = []
                for _ in range(args.repeats):
                    start = time.perf_counter(); count, _ = method(); samples.append(time.perf_counter() - start)
                row["baselines"][name] = {"seconds": samples, "median_seconds": statistics.median(samples), "count": count}
            if len({record["count"] for record in row["methods"].values() if "count" in record}
                   | {record["count"] for record in row["baselines"].values()}) != 1:
                raise ArithmeticError("Counting baselines disagree")
            report["cases"].append(row)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
            print(json.dumps({"case": case, "seconds": {name: value.get("medians", value["status"])
                              for name, value in row["methods"].items()}}), flush=True)
    else:
        configurations = [((8, 71, Q, R), mode) for Q, R in ((4, 2), (8, 2), (8, 5), (12, 3), (20, 3))
                          for mode in ("reference", "balanced") if Q <= 12 or mode == "balanced"]
        configurations += [((20, 71, Q, R), "balanced") for Q, R in ((4, 2), (8, 3), (12, 3), (20, 3), (12, 7))]
        for case, schedule in configurations:
            for method in ("spectral", "required"):
                print(f"Inspecting {case}, {schedule}, {method}...", flush=True)
                sample = run_child({"case": case, "method": method, "schedules": schedule, "solve": False}, args.timeout)
                report["cases"].append({"case": case, "schedules": schedule, "method": method, **sample})
                args.output.parent.mkdir(parents=True, exist_ok=True)
                args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output.resolve(), flush=True)


if __name__ == "__main__":
    main()
