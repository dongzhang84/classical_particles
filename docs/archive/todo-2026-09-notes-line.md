> **已归档（2026-11-24）。** 这是 2026-09 笔记那条线（`notes/notes.pdf`，ζ 参数与观测约束）的待办清单。2026-11-18 起主线改为带电小黑洞的一阶运动，笔记线暂停（见提交 `b4ea336`），这份清单不再维护。其中第 2 项（附录 A.3 的 `T[F_S]` 积分独立核对）在 `README.md` 里仍标为"待独立核对"，笔记线恢复时从这里接着做。

# TODO

Last updated: 2026-09-02

## Current Checkpoint

The active deliverable is `notes/notes.pdf`. The LaTeX source and figures live
under `notes/`, but repository policy is to commit only `notes/notes.pdf`.

The current PDF has completed the first pass through Sec. 7. The main result is
that the useful null-test parameter is the radiation-reaction strength
coefficient `zeta`, not a revived alternative equation of motion. In the
classical regime, storage-ring damping data constrain `zeta` most directly. In
the transition regime, experiments measure the product `zeta g(chi)`, so a scan
in `chi` is needed to separate the classical coefficient from quantum
suppression.

Latest committed checkpoint:

- `45d90c4` Correct the classical boundary in Fig. 7, and read zeta off the SLC numbers

## Next Tasks

### 1. Rebuild the Storage-Ring `zeta` Bound

Goal: turn the SLC agreement quoted in Sec. 7.4 into a defensible bound, or
state clearly why the public data are insufficient.

Tasks:

- Locate the original SLC damping-ring paper and any machine-parameter tables
  needed to reconstruct `U0,cl`, damping partition numbers, beam energy, and RF
  calibration uncertainties.
- Separate measured quantities from model inputs.
- Estimate or quote uncertainties for `U0meas/U0calc` and
  `tau_calc/tau_meas`.
- Decide whether Sec. 7 should quote an interval for `zeta` or only a
  percent-level consistency check.

Done when:

- Sec. 7.4 contains either a real uncertainty budget or an explicit statement
  that the paper only supports an order-of-magnitude consistency estimate.
- `notes/notes.pdf` has been rebuilt and committed by itself.

### 2. Independently Check the Appendix A.3 `T[F_S]` Integral

Goal: settle whether the `2/3 dot(a)` term attributed to the singular field in
the 2006 thesis is real, a bookkeeping artifact, or a mistaken identification
of field momentum flow with force.

Tasks:

- Recompute the relevant stress-energy flux in retarded coordinates.
- Compare the result with the Detweiler-Wighting / Harte / Pound statement
  that the singular field renormalizes mass but exerts no force.
- Update `REVIEW_2026.md`, and if needed the historical discussion in
  `notes/notes.tex`, so the status of the lambda-kappa motivation is precise.

Done when:

- The README no longer has an unresolved "待独立核对" note, or it names the
  unresolved subproblem exactly.
- The text distinguishes the phenomenological usefulness of `zeta` from the
  historical claim that it came from the 2006 singular-field calculation.

### 3. Improve the `(B, gamma)` Parameter Map

Goal: make the regime figure self-contained enough that a reader can see where
the classical, transition, and QED regions are without reconstructing the
scalings from the text.

Tasks:

- Check the boundaries currently shown in Fig. 7:
  `tau0 omega_B = 1`, `B = B_S`, `chi = 0.1`, `chi = 1`, and any
  radiation-reaction-limited curves.
- Add representative points or bands for storage rings, laser-electron
  experiments, pulsars, magnetars, and FRB-relevant outer magnetospheres if the
  numbers can be sourced cleanly.
- Ensure the caption says which regions are usable for classical `zeta` fits
  and which are exclusion zones.

Done when:

- Fig. 7 supports the Sec. 7 argument without relying on memory of
  `REVIEW_2026.md`.
- The figure and caption survive a fresh PDF rebuild.

### 4. Decide the Next Paper Shape

Goal: avoid letting the project sprawl back into "solve radiation reaction."

Candidate narrow title:

> Constraints on a phenomenological radiation-reaction strength in the
> classical and quantum-transition regimes

Proposed scope:

- Define `zeta` only as a null-test parameter with predicted value 1.
- Use storage rings as the clean classical bound.
- Use laser-electron experiments only where a `chi` scan can separate
  `zeta` from `g(chi)`.
- Use magnetars primarily to mark the failure of the classical theory, not as a
  classical fitting target.

Done when:

- The abstract and final section of `notes/notes.tex` say the same narrow
  thing.
- No section implies that the lambda-kappa family is a viable replacement for
  LD/LL.

## Workflow Reminders

- After any change to `notes.tex` or `notes/figures/`, rebuild from `notes/`:

```sh
pdflatex notes && bibtex notes && pdflatex notes && pdflatex notes
```

- Commit only the PDF:

```sh
git add notes/notes.pdf
git diff --cached --name-only
```

The second command must print exactly:

```text
notes/notes.pdf
```

## Recent Changelog

- Completed Sec. 7 as the `zeta` observability section.
- Added the reduced-order equation with a free radiation-reaction coefficient.
- Added the sensitivity table showing which observables scale linearly or weakly
  with `zeta`.
- Added the verdict table: storage-ring damping is the tightest clean
  constraint; synchrotron burn-off is astrophysically interesting but weak;
  cooling times and curvature quantities are too degenerate to be useful.
- Added the storage-ring damping figure and the auxiliary-quantity degeneracy
  figure.
- Added Sec. 7.5 explaining why `zeta` and `g(chi)` multiply in the
  quantum-transition regime and why a scan in `chi` is needed.
- Corrected the classical boundary in Fig. 7 and used SLC numbers to read off
  `zeta` at the percent-consistency level.
