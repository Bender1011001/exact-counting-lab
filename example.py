"""Small reproducible example: python example.py (standard library only)."""
import random
import contraction
import implicit_design
import modular_lab
import reconstruction_lab as lab


def main():
    cnf = lab.ref.random_bounded_occurrence_3cnf(6, 6, random.Random(1))
    design = implicit_design.build_implicit_design(cnf, 2, 2, 2, random.Random(1))
    order = design.M * design.degP // design.Q
    count, certificate = modular_lab.exact_modular_count(
        design, cnf.n, order, reconstruction=contraction.TaylorContraction(design, cnf.n),
        use_disjoint_bound=True)
    brute = lab.ref.count_solutions(cnf)
    ordinary, _ = contraction.conventional_count(cnf)
    if count != brute or count != ordinary or count != 32:
        raise RuntimeError((count, brute, ordinary))
    print("Six-variable synthetic formula, seed 1")
    print(f"Polynomial pipeline: {count}")
    print(f"Brute force: {brute}")
    print(f"Ordinary elimination: {ordinary}")
    print(f"CRT primes used: {len(certificate['primes'])}")
    print("Agreement checks passed. Use the benchmark runners for timings.")


if __name__ == '__main__':
    main()
