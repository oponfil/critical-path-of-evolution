# The Critical Path of Evolution

A theoretical, interdisciplinary preprint that proposes a new resolution of the **Fermi paradox** — not "we are alone" and not "others are far ahead," but **synchronization**: most civilizations emerge at approximately the same cosmic time.

**Author:** Oleg Ponfilenok (independent researcher)

---

## Overview

The work treats the evolution of life — from the thermal gate of the cosmic microwave background to a technological civilization, over a time *T* ≈ 13.8 Gyr — as a single chain of *N* ~ 10¹⁴ small, **high-probability** steps (a *critical path*: the rate-limiting, fastest route of increasing complexity), rather than a sequence of rare **"hard steps."**

Because the step times are independent, the concentration of their sum (additivity of variance) compresses the *relative* spread of "exit times" as 1/√*N*. Different planets therefore reach civilization almost simultaneously. The observed absence of "loud" civilizations (**cosmic silence**) is read as an **upper bound** on the interplanetary dispersion of exit times, σ, and on the number *n* of planets currently on the critical path.

## Core idea

- **Critical path, not hard steps.** Evolution follows the single fastest route through complexity (kinematically unique; selected thermodynamically by maximum entropy production). Apparent "key transitions" (eukaryogenesis, multicellularity) are long chains of simple sub-steps, not single improbable jumps.
- **Synchronization.** Summing many independent steps narrows the spread of emergence times to σ ~ 10³ yr — negligible on cosmological scales — so no civilization is meaningfully "older" than the rest.
- **Two cross-checks meet.** The dispersion inferred from silence (astronomy) and the elementary "tick" inferred from microbiology agree in order of magnitude.

## Key quantitative results

| Quantity | Value | Meaning |
|---|---|---|
| σ = √(μT) | ≈ 800–1400 yr | upper bound on interplanetary dispersion of exit times (from silence) |
| μ = σ²/T | ≈ 25–75 min | mean elementary step time ≈ microbial generation time |
| *n* | ~ 10⁴–10⁵ | civilization-bearing planets in the Galaxy (matches Scherf–Lammer "Earth-like Habitats") |
| *t*₍loud₎ | ~ 10³ yr | time for our civilization to reach Kardashev Type II |

## Testable predictions

- **Eukaryogenesis (decisive test):** it should decompose into a chain of high-probability sub-steps; a single indecomposable, slow link would break synchrony.
- **Planet age ↔ stage of life:** Earth-like planets much older than Earth should host advanced (not simple) life, within σ ~ 10³ yr of our level; simple microbial life is expected only on younger planets.
- **Convergent form:** technological civilizations are expected to converge on a *functional* humanoid plan (manipulator limbs, large social brain, tool use).
- **Kardashev climb:** continued Type I → II growth on a ~10³ yr horizon.
- **Type I detectability:** an instrument able to detect Type-I civilizations interstellarly would expect ~5 already-visible neighbors (a null result does *not* falsify the model, since *n* is only an upper bound).

## Repository contents

| File | Description |
|---|---|
| `critical_path_evolution_en.tex` | English version of the paper (LaTeX source) |
| `critical_path_evolution_ru.tex` | Russian version of the paper (LaTeX source) |
| `critical_path_evolution_en.pdf` | Compiled English PDF |
| `critical_path_evolution_ru.pdf` | Compiled Russian PDF |
| `refs.bib` | Shared bibliography (used by both versions) |

## Building the PDFs

Requires a TeX distribution (TeX Live or MiKTeX). The English version uses Times fonts (`mathptmx`); the Russian version uses Cyrillic support (`babel`/`T2A`) and the `comfortaa` font package.

```sh
# English
pdflatex critical_path_evolution_en
bibtex   critical_path_evolution_en
pdflatex critical_path_evolution_en
pdflatex critical_path_evolution_en

# Russian
pdflatex critical_path_evolution_ru
bibtex   critical_path_evolution_ru
pdflatex critical_path_evolution_ru
pdflatex critical_path_evolution_ru
```

The two passes after `bibtex` resolve citations and cross-references.

## Status

Working preprint — a hypothesis paper spanning cosmology, non-equilibrium thermodynamics, population genetics, astrobiology, and SETI. Feedback and criticism are welcome.

© 2026 Oleg Ponfilenok. All rights reserved.
