# v0.1.0 — experimental research release

Independent exact-counting experiments extending Chelsea A. McElveen's MIT-licensed community reference for the OpenAI ETH preprint.

Included: exact input-specific auxiliary construction, spectral and required measurement bases, cardinality-weighted reconstruction, cofactor pruning, finite-field interpolation and certified CRT recovery, a conventional contraction backend, synthetic benchmarks with established baselines, derivations, raw results, and 4,639 verification checks. Run `python example.py` for a six-variable example with count 32, then `python verify.py`.

Controlled full-pipeline improvements against the eager constructor with its expensive global membership audit disabled were 1.54–10.68x on three small cases. Conventional methods remain faster on those cases. Larger apparent reference speedups include removing exhaustive self-check overhead. Negative results and degeneracies are retained.

This is a research prototype. No general solver superiority, new asymptotic bound, mathematical novelty, full-paper proof certification, or security-analysis advantage is established. It implements experiments on the rational-interpolation route, not the sharper Appendix D algorithm. Cached local assignments can use exponential space. The finite prime pool and explicit construction limits constrain supported inputs; the modular adapter is serial and not thread-safe.

Historical timings were collected locally on October 30, 2026; see raw JSON for samples and environment. Attribution, source commit and SHA-256 are in README.md; the upstream license is retained in third_party/LICENSE.txt. No proprietary fixtures or licensed applications are included.
