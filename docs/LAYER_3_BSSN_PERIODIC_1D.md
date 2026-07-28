# Layer 3.4 — BSSN on a periodic one-dimensional grid

Layer 3.3 delivered the BSSN algebra at a *point*. Every spatial-derivative term
in the evolution equations vanished identically there, because every field was
spatially constant. This increment supplies the missing ingredient: a spatial
grid. It is the first increment in which BSSN is a partial differential equation
rather than an algebraic identity.

It is also deliberately the smallest such increment that can be validated
honestly.

## 1. Why one spatially varying coordinate

The obvious next step after Layer 3.3 is "a 3D grid". That step is not taken
here, for the same reason Layer 3.2 declined to put ADM on a grid: a result that
cannot be checked against something independent is not a result.

A one-dimensional periodic domain buys three things a 3D domain does not:

- **Exact analytic derivative oracles.** On a periodic domain a smooth field's
  derivatives are known in closed form, so the finite-difference operators can be
  validated against `sin(kx)` rather than against themselves.
- **No boundary problem.** Periodicity is an exact boundary condition. An outer
  radiative boundary is an open research problem and a large source of error
  that would contaminate every convergence measurement made here.
- **An exact evolution oracle.** A linearized transverse-traceless wave
  propagating along `x` is an exact solution of *linearized* general relativity
  in this gauge, with a closed-form time dependence to compare against.

Fields vary only along `x`, but every grid point stores complete
three-dimensional tensors — `gammatilde_ij` and `Atilde_ij` are full symmetric
3x3 objects, `Gammatilde^i` is a full 3-vector. This is a **1D3V reduction**:
one spatial dimension, three-component vectors, three-by-three tensors. It is
**not** a general three-dimensional numerical-relativity solver and must not be
described as one.

## 2. How this reuses Layer 3.3 — and the measurement that decided it

Every Layer 3.1 and Layer 3.3 evaluator has the same shape:

```rust
pub fn bssn_evolution_rhs<G: Metric<3>>(
    spatial_metric: &G, extrinsic_curvature: &impl SpatialTensorField,
    coordinates: &[f64; 3], gauge: &BssnGauge, sources: &AdmSources,
    settings: &AdmEvolutionSettings,
) -> Result<(BssnState, BssnEvolutionRhs), BssnError>
```

The fields are **sampled by coordinate**, and all spatial derivatives are taken
internally by nested central differences at `settings.spatial_step` and
`settings.metric_step`. There is no seam through which precomputed grid
derivatives can be injected.

That left two options: refactor `bssn.rs` to accept injected derivatives, or
supply a field that *reads from the grid* and let the existing machinery
difference it. The second duplicates nothing, but risks a real pathology: the
nested difference samples `x ± 2 dx`, and a second difference of spacing `2 dx`
can decouple the even and odd grid points into independent sub-lattices — the
classic checkerboard mode.

This was measured before it was decided. A grid-lookup `Metric<3>` (nearest
periodic index) was fed to `conformal_ricci` with
`spatial_step = metric_step = dx`, against a converged analytic reference:

| N   | dx       | Linf(R_grid − R_exact) | ratio | even points | odd points |
| --- | -------- | ---------------------- | ----- | ----------- | ---------- |
| 32  | 0.031250 | `2.574330e-3`          | —     | `2.574330e-3` | `2.522839e-3` |
| 64  | 0.015625 | `6.462596e-4`          | 3.983 | `6.462596e-4` | `6.430096e-4` |
| 128 | 0.007812 | `1.617200e-4`          | 3.996 | `1.617200e-4` | `1.615198e-4` |
| 256 | 0.003906 | `4.042696e-5`          | 4.000 | `4.042696e-5` | `4.041133e-5` |

Two conclusions about **accuracy**, both load-bearing:

- The error ratio converges to **exactly 4.000** per halving. The grid path is
  second-order accurate.
- Even- and odd-index *errors* agree to **0.04%** at `N = 256`, so the
  truncation error shows no even/odd split on smooth data.

On that basis the grid-provider route was adopted, and `bssn.rs` is **not
modified at all**: no BSSN equation, no Ricci engine, no constraint evaluator,
and no integrator is duplicated by this increment.

**That probe measured accuracy, and accuracy was not the whole question.** A
separate stability measurement, made once evolution was running, found the same
stencil to be unstable — see section 13, which is the most important section of
this document. The two results are not in conflict: accuracy is a statement
about smooth data, stability is a statement about the highest frequency the grid
can represent, and one does not imply the other.

This also answers a question Layer 3.3 could not. There, the conformal Ricci
decomposition agreed with the independently computed physical Ricci to `~1e-6`,
and that was correctly described as a *nested finite-difference floor* — a fixed
number with no error model. On a grid it is not a floor at all: it is truncation
error, and it converges at the discretization order.

**The consequence for callers is a hard requirement**: on the grid,
`spatial_step` and `metric_step` **must both equal `dx`**. Any other value makes
the internal difference sample coordinates that are not grid points, where the
provider can only return its nearest neighbour — first-order garbage. The grid
system sets these itself and does not take them from the caller.

## 3. Grid convention

A **half-open** uniform periodic domain:

```text
x_n = x_min + n * dx,    n = 0, ..., N-1,    dx = (x_max - x_min) / N
```

The upper endpoint `x_max` is **not** stored: it is the periodic image of
`x_min`. Storing both would duplicate a degree of freedom and make the
determinant of the discrete system singular.

Periodic indexing wraps with Euclidean remainder, never with unsigned integer
wraparound:

```text
wrap(i) = i.rem_euclid(N)
```

so the left neighbour of index `0` is `N-1` and the right neighbour of `N-1` is
`0`, with no underflow possible on any offset, positive or negative.

Rejected at construction: fewer points than the stencil requires, non-finite
bounds, non-positive domain length, and a non-finite or non-positive spacing.

## 4. Finite-difference operators

Second-order centred, periodic:

```text
D1 f_i = ( f_{i+1} - f_{i-1} ) / ( 2 dx )
D2 f_i = ( f_{i+1} - 2 f_i + f_{i-1} ) / ( dx^2 )
```

Validated against `f(x) = sin(kx)`, whose derivatives `k cos(kx)` and
`-k^2 sin(kx)` are exact, at several wave numbers below the Nyquist limit and
over at least four resolutions. Second-order convergence is **measured**, not
asserted from the stencil formula.

A fourth-order stencil is deliberately **not** implemented. It is not required
for a validated second-order pipeline and would compete for validation effort
with the parts of this increment that are actually new.

These operators are used for diagnostics and for the manufactured-state oracle.
They are **not** the differencing used inside the BSSN right-hand side — that is
Layer 3.3's nested difference, characterised in section 2. The documentation
says so plainly rather than implying a single unified stencil.

## 5. Field storage and the flat state

`BssnGridState` stores **17 scalar component arrays** in structure-of-arrays
order. Component `c` occupies the contiguous slice `[c*N, (c+1)*N)`:

| slot | component | slot | component |
| ---- | --------- | ---- | --------- |
| 0 | `phi` | 8..13 | `Atilde_xx, xy, xz, yy, yz, zz` |
| 1..6 | `gammatilde_xx, xy, xz, yy, yz, zz` | 14..16 | `Gammatilde^x, ^y, ^z` |
| 7 | `K` | | |

Symmetric rank-2 tensors store only their six independent components, so
symmetry is preserved **structurally** — it cannot be broken by a round trip,
because the redundant entries do not exist.

Structure-of-arrays is chosen because a component derivative reads `f[i-1]`,
`f[i]`, `f[i+1]` at stride 1, and because the flat layout `scirust_sim` requires
is then the storage itself: encoding and decoding are contiguous copies, not
scattered gathers. There is no separate "flattening" step to get wrong.

The flat method-of-lines state is exactly this array, so its length is `17 N`.

## 6. Method of lines

The semidiscrete system is `dY/dt = F(t, Y)` with `Y` the `17 N` array above.
Time integration reuses `scirust_sim::simulate` — the platform's existing
deterministic fixed-step RK4 — through its existing `System` trait, exactly as
Layer 3.2 and `GeodesicSystem` already do. **No integrator is written here.**

Per right-hand-side evaluation the grid system loops over points in ascending
index order (deterministic), builds the local ADM pair by reconstruction, calls
Layer 3.3's `bssn_evolution_rhs`, and writes the 17 components into the output
slice. Reductions accumulate in ascending index order, so they are
bit-reproducible.

## 7. Gauge

Prescribed only: `alpha = 1`, `beta^i = 0`. Gauge is **not evolved**. There is no
1+log slicing, no Gamma-driver shift, no gauge damping, and no evolved gauge
variable. This is a **major limitation**, not an implementation detail: live
gauge is what makes strong-field BSSN evolution work in practice, and none of it
is present.

## 8. Constraint monitoring

Free evolution. Constraints are **monitored and never enforced**. Projection
exists as an explicit opt-in experiment mode, disabled by default, and is never
applied between RK4 stages in the default evolution.

Monitored separately — never blended into a single score:

- BSSN algebraic: unit conformal determinant, trace-free `Atilde`, and
  `Gammatilde^i` consistency.
- ADM physical, from reconstructed fields: the Hamiltonian constraint and the
  three momentum-constraint components.

Each is reduced over the grid to a signed mean, `L1`, discrete `L2`, `Linf`, the
index attaining the maximum (lowest index wins ties, deterministically), and a
non-finite count.

## 9. Validation oracles

- **A — stationary Minkowski.** Identity conformal metric, `phi = 0`, `K = 0`,
  `Atilde = 0`, `Gammatilde^i = 0`, vacuum. Every right-hand side and every
  constraint must stay at machine level for many steps, with no drift.
- **B — manufactured periodic state.** A smooth periodic BSSN field with known
  derivatives. It does **not** satisfy the Einstein equations and is labelled a
  manufactured numerical oracle, not a physical solution. It validates packing,
  derivative extraction, and pointwise assembly, and it measures spatial
  convergence.
- **C — linearized transverse-traceless wave.** `h_yy = A sin(k(x-t))`,
  `h_zz = -A sin(k(x-t))`, with `K_ij = -(1/2) d_t gamma_ij` following the
  repository's sign convention. This is an exact solution of *linearized*
  general relativity in this gauge. The code evolves the **full nonlinear** BSSN
  system, so the measured difference contains both `O(A^2)` nonlinearity and
  `O(dx^2)` truncation; the amplitude is chosen so truncation dominates and the
  convergence measurement is meaningful. This is **not** nonlinear-wave
  validation.
- **D — off-constraint state.** A controlled periodic constraint violation,
  evolved freely, showing the expected initial residual and the absence of
  silent repair.
- **E — algebraic projection.** Injected determinant and trace violations,
  measured, explicitly projected, and compared against the unprojected
  evolution.

## 10. Numerical dissipation

Kreiss–Oliger dissipation is implemented but **disabled by default**
(`sigma = 0`), and is never enabled implicitly:

```text
Q f_i = -(sigma / 16 dx) ( f_{i+2} - 4 f_{i+1} + 6 f_i - 4 f_{i-1} + f_{i-2} )
```

The bracket's Fourier symbol is `4 (1 - cos theta)^2`: zero for a constant
field, maximal (`16`) at the Nyquist mode. For smooth `f` the bracket is
`O(dx^4)`, so `Q f = O(dx^3)` and second-order accuracy is preserved.

It is **numerical dissipation on the evolved variables**, not constraint
damping: it neither targets nor reduces the constraint residuals by
construction.

The undissipated scheme was measured first, as it must be — dissipation added
before measurement hides an instability rather than revealing one. And in this
case dissipation demonstrably **does not cure** the instability; see section 13.

## 11. Determinism

No RNG, no wall clock, no parallelism, no hidden global state. Point iteration
and reduction order are fixed and ascending. Identical inputs produce
byte-identical output, verified by running the experiment twice and comparing
bytes.

## 12. Measured results

All figures from `bssn_periodic_1d_evolution`, byte-identical across runs.

**Oracle A — stationary Minkowski.** At `N = 16, 32, 64, 128` the state change
after 1 unit of coordinate time is `0.000000e0` — *exactly* stationary, not
merely small. The right-hand side of a spatially constant field is exactly zero
because a centred difference subtracts identical values, so RK4 adds exactly
zero and there is no rounding to accumulate. Every constraint is `0.000000e0`.

**Oracle B — manufactured state, conformal Ricci reconstruction.**

| N | dx | mismatch L1 | L∞ | ‖R̃‖ | ‖R^φ‖ | order |
|---|---|---|---|---|---|---|
| 32 | 3.125e-2 | 6.224e-5 | 9.891e-5 | 1.961e-1 | 2.599e-3 | — |
| 64 | 1.563e-2 | 1.601e-5 | 2.521e-5 | 1.981e-1 | 2.624e-3 | **1.97** |
| 128 | 7.813e-3 | 4.029e-6 | 6.333e-6 | 1.986e-1 | 2.630e-3 | **1.99** |
| 256 | 3.906e-3 | 1.009e-6 | 1.585e-6 | 1.987e-1 | 2.632e-3 | **2.00** |

Layer 3.3 could only report this mismatch as a fixed `~1e-6` *floor* with no
error model. On a grid it is truncation error and converges at the
discretisation order. Both parts are genuinely nonzero, so the agreement is not
two zeros matching.

**Oracle C — linearized wave, short-time spatial convergence** (`A = 1e-6`,
`k = 2π`, `t = 0.1`, `C = 0.25`): metric `L∞` error `9.373e-9 → 2.416e-9 →
6.063e-10 → 1.520e-10` at `N = 16..128`, observed order **1.96 / 1.99 / 2.00**.

The wave's *constraint* residual behaves differently, and the difference is
physical rather than numerical. At fixed `A = 1e-6` the Hamiltonian residual is
flat in resolution (`7.74e-11 → 8.98e-11` across `N = 32..256`), but at fixed
`N = 64` it scales as `A^2` — measured exponents **3.96 / 4.01 / 4.00 per
amplitude doubling**. It is the quadratic nonlinearity the linearized oracle
neglects, not discretisation error. The truncation contribution cancels because
the linear-order Ricci *scalar* is proportional to `d^2(h_yy + h_zz)` and
`h_yy + h_zz = 0` holds exactly at every grid point. Asserting `dx`-convergence
of this quantity would have been asserting something false.

**Temporal convergence.** Measured against a finely-stepped run at the *same*
spatial resolution, so the spatial truncation error cancels exactly and what
remains is purely the RK4 error of the semidiscrete ODE system. Difference `L∞`
`7.169e-9 → 4.490e-10 → 2.807e-11 → 1.768e-12`, observed order **4.00 / 4.00 /
3.99**. Measuring against the analytic PDE solution instead would be dominated
by `dx^2` and would report a misleadingly low order.

**Benchmarks** (machine-dependent wall clock; the computation is deterministic):

| operation | N=32 | N=64 | N=128 |
|---|---|---|---|
| periodic first derivative | 73.0 ns | 142.0 ns | 282.6 ns |
| periodic second derivative | 71.6 ns | 143.7 ns | 297.3 ns |
| full BSSN grid RHS | 930 µs | 1.84 ms | 3.66 ms |
| constraint monitor | 1.07 ms | 2.13 ms | 4.23 ms |

Both the derivatives and the right-hand side scale linearly in `N` (RHS ratios
1.98 and 1.99 under doubling), as a one-dimensional grid should. Single-point
costs: conformal connection 3.55 µs, conformal Ricci 25.1 µs. Flat-state round
trip at `N = 128` is 181.6 ns — a memcpy, because storage *is* the flat layout.
Initial data from analytic ADM fields at `N = 64` costs 88.3 µs. Ten RK4 steps
at `N = 32` cost 36.7 ms, i.e. 3.67 ms per step — exactly four right-hand-side
evaluations, confirming no hidden work per stage.

Per grid point the RHS costs ~28.6 µs against Layer 3.3's ~10 µs local figure;
the difference is the `bssn_to_adm` reconstruction the provider performs on
every nested-difference sample. `bssn_grid_rhs` performs **zero heap
allocations** — the input is a borrowed view and the output is written in place,
and the providers are `Copy`. Storage is 17 f64 = **136 bytes per grid point**;
`simulate` adds six full-state buffers during integration.

## 13. The reused stencil is accurate but not stable

This is the increment's most important result, and it is negative.

Composing an outer and an inner difference, each of width `dx`, produces spatial
operators spanning `2 dx`. The Fourier symbol of a `2 dx`-spaced second
difference is `2 cos(2 theta) - 2`, which **vanishes at `theta = pi`**. The
Nyquist mode — the highest frequency the grid can represent — lies in the null
space of the principal part. It is invisible to the term that should control it,
and it grows.

Measured, at `C = 0.25` unless noted:

| N | outcome |
|---|---|
| 16 | survives to `t = 1` |
| 32 | survives to `t = 1` at **every** Courant factor from `0.1` to `2.0` |
| 64 | fails at `t ≈ 0.98`, at **every** Courant factor from `0.1` to `2.0` |
| 128 | fails at `t ≈ 0.50`, at **every** Courant factor from `0.1` to `2.0` |

Two features identify this as a spatial instability rather than a CFL violation:
the onset time is **independent of the timestep** across a factor of twenty in
`dt`, and it roughly **halves as the resolution doubles**. Refining the grid
makes it worse. Tracking the Nyquist projection of `gammatilde_yy` shows it
rising out of rounding noise (`1e-16`) through `5.2e-9` and then diverging.

Minkowski is the control: its right-hand side is exactly zero, nothing seeds the
mode, and it stays stationary indefinitely (`N = 128` to `t = 5`, Nyquist
projection exactly `0.000e0`).

**Dissipation does not cure it.** At `N = 128`, raising `sigma` from `0` to
`0.5` moves the onset only from `t = 0.50` to `t = 0.60`. At `N = 64`,
`sigma = 0.05` avoids the abort but inflates the physical wave amplitude by
`2.3e3`; `sigma = 0.2` by `53`; `sigma = 0.5` by `1.96`. "Did not reach
infinity" is not "stable". Since the fourth-order Kreiss–Oliger operator damps
`theta = pi` hardest, the fact that it fails shows the unstable content is not
confined to the Nyquist mode.

**So no dissipation is enabled by default, and none is used to make any result
in this increment look better than it is.**

The usable envelope is therefore coarse grids and short times, and every
convergence measurement above was taken inside it. The principled fix is to
remove the `2 dx` span by giving Layer 3.3 a derivative-injecting entry point so
the grid can supply true `dx`-spaced stencils. That is a refactor of `bssn.rs`
and is deliberately **not** attempted here; it is the recommended next
increment.

## 14. Known limitations

Stated plainly, because the gap between this and a numerical-relativity code is
large:

- One spatially varying coordinate only. Transverse derivatives are exactly zero
  **by construction**, not by approximation.
- Periodic domain only. No outer boundary, radiative or otherwise.
- Prescribed gauge only. No live gauge of any kind.
- Weak, smooth fields only. No singular spacetimes.
- **The scheme is not stable outside the envelope in section 13.** `N >= 64`
  fails at `t < 1` at every Courant factor tested. This is a working, validated
  *pipeline*, not a usable evolution code.
- **Strong hyperbolicity is not proven by these tests.** Passing a weak-field
  convergence test is evidence, not proof; hyperbolicity is a property of the
  PDE system established analytically, and nothing here establishes it. Nor does
  the measured instability disprove it for BSSN — the failure is in this
  *discretisation*, not in the continuum formulation.
- No general three-dimensional validation.
- No black holes, no punctures, no excision.
- No adaptive mesh refinement.
- No waveform extraction.
- No observational validation.
