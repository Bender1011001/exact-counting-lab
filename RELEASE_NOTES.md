# v0.2.0 — experimental field fibers and Taylor planning

This release adds two opt-in components and improves the research documentation.

- Finite-field measurement preparation (`measurement_mode="field"`) rebuilds the required basis modulo each prime and retains assignment multiplicities. A residue-only control is also available.
- Certified Taylor planning (`use_form_degree=True`) bounds the degree from the actual occurrence forms and requests sufficient Hermite conditions.
- The full Python suite passes 6,430 checks. An independent Mathematica 15.0.1 expansion agrees on 63 coefficients.
- The six-case, three-repeat benchmark includes 216 complete runs, exact-count and form-hash audits, conventional baselines, small-prime controls, and separate measurement-domain diagnostics.
- README now includes a quick start, CI/release/license badges, a pipeline diagram, balanced results, and contribution/citation guidance.

Combined full-pipeline gains on the three balanced-schedule cases were 6.62x, 8.08x, and 2.58x. Two reference-schedule cases became 3.93x and 2.79x slower because lower Taylor orders disabled cofactor pruning. Small-prime runs were slower overall, and no measurement-domain reduction occurred in the measured cases. The new modes therefore remain opt-in.

Conventional counting methods still win all six cases. These are component experiments on small synthetic inputs, including low-degree and constant deformations. No general solver superiority, new asymptotic bound, novelty, or full-paper certification is claimed. See FIELD_FIBERS.md for exact methodology, ablations, and limitations.

Run `python example.py` and `python verify.py`. Benchmarks additionally use the optional pinned BDD requirement. Wolfram is optional; the core remains standard-library Python. The attributed upstream reference remains unchanged.
