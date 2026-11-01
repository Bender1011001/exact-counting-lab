"""Run the portable package's independent exact checks."""
import json
import sys
import reconstruction_lab as lab
import modular_lab as modular
import cofactor_checks
import contraction_checks
import bound_checks
import implicit_checks

if __name__ == "__main__":
    if sys.flags.optimize:
        raise SystemExit("Run without -O: the inherited research checks use assertions.")
    report = {"coefficient_checks": lab.verify(), "modular_checks": modular.verify(),
                      "cofactor_checks": cofactor_checks.verify(),
                      "contraction_checks": contraction_checks.verify(),
                      "certified_bound_checks": bound_checks.verify(),
                      "implicit_design_checks": implicit_checks.verify()}
    report["total_checks"] = (report["coefficient_checks"] + report["modular_checks"]
                              + report["cofactor_checks"]["checks"]
                              + sum(report["contraction_checks"].values())
                              + report["certified_bound_checks"]
                              + sum(report["implicit_design_checks"].values()))
    print(json.dumps(report))
