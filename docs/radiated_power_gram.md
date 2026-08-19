# Radiated-power Gram (WASM status)

Handoff for [approximate_matched_basis.md](approximate_matched_basis.md). **§5–9 (through \(T\)) are implemented**: isolated fields → Hermitian \(P_H\) (product kernel and \(J_0\)), \(R\), optional \(jX(\Delta x,\Delta y)\), per-port match (real \(z_0\), or complex conjugate when \(X_{nn}\neq 0\)) or a common complex \(z_c\), power-wave \(S\) and \(T\). The visualizer has Isolated/Matched coupling, Per-port/Common Z0 match style, and Geometric/Conjugate steer. \(F^\mathrm{emb}\) is not formed on a quadrature grid (the GUI uses \(w=Ta\)).

Intent: keep the **φ-dependent** product quadrature as the general kernel (needed later for rotated / polarized patterns). Do not auto-switch `compute()` or `runPradJob` to J0. Use `compute_j0` only for planar axisymmetric \(D(\mu)\).

---

## What \(P_H\) is

For identical co-aligned planar elements, the doc’s dual-pol Gram collapses to a scalar steering matrix. Positions \((x_p,y_p)\) are in wavelengths at \(f_0\). Electrical phase uses the same convention as the far-field kernel: \(k = 2\pi\cdot\texttt{frequency_scale}\).

Power-conserving **directivity** \(D\) with \(\int D\,d\Omega = 4\pi\), front hemisphere only:

| `element_kind` | \(D(\mu)\) for \(\mu=\cos\theta\ge 0\) |
| --- | --- |
| `PATTERN_ISOTROPIC` (0) | \(2\) (+3 dBi hemisphere) |
| `PATTERN_COS_N` (1) | \(2(n+1)\,\mu^n\); \(n=0\) matches isotropic |

Isolated field at quadrature sample \(s\):

\[
A_{p,s}=\sqrt{\omega_s D(\hat r_s)}\,\exp(-j\,k\,\hat r_s\cdot\mathbf r_p).
\]

Then, with peak-wave \(P_0=0.5\,\mathrm{W}\),

\[
P_H=\frac{P_0}{4\pi}\,A A^H
\quad\text{(Hermitianized)}.
\]

\(\alpha\), \(\eta_0\), and \(r_\mathrm{ref}\) cancel because \(D\) is already power-normalized. Diagonal \(P_{H,pp}\approx P_0\). For these axisymmetric planar elements \(P_H\) is real (imaginary part is roundoff).

**\(N\)** = number of elements. **\(M\)** = number of quadrature **directions**, \(M=n_\mu n_\phi\), not array size. \(A\) is \(N\times M\). Cost of the Gram is \(O(N^2 M)\).

---

## Quadrature (Gauss-μ × φ)

Front hemisphere, true spherical measure \(d\Omega=d\mu\,d\phi\) (no \(\sin\theta\)):

- **μ:** Gauss–Legendre on \([0,1]\) (map standard nodes \(\xi\in[-1,1]\) by \(\mu=(\xi+1)/2\), \(w_\mu=w_\xi/2\)).
- **φ:** \(n_\phi\) equal nodes on \([0,2\pi)\), trapezoid weight \(\Delta\phi=2\pi/n_\phi\).
- Sample order: **μ outer, φ inner**. Direction \(\hat r=(u,v,\mu)\) with \(u=\sqrt{1-\mu^2}\cos\phi\), \(v=\sqrt{1-\mu^2}\sin\phi\).
- Weight \(\omega_{ij}=w_{\mu,i}\,\Delta\phi\).

Do **not** reuse the display far-field grids (non-standard \((\theta,\phi)\) charts). Rule-of-thumb orders (not auto-applied):

\[
n_\mu \approx \lceil\pi D_\lambda f_\mathrm{scale}\rceil+12,\qquad n_\phi\approx 2 n_\mu\text{ (even)},
\]

where \(D_\lambda\) is the array electrical diameter in wavelengths.

---

## Files

| Path | Role |
| --- | --- |
| [wasm/src/quadrature.rs](../wasm/src/quadrature.rs) | Gauss–Legendre + `HemisphereQuad` (1-D μ kept for J0) |
| [wasm/src/bessel.rs](../wasm/src/bessel.rs) | f32 \(J_0\) (fdlibm `j0f`) |
| [wasm/src/prad.rs](../wasm/src/prad.rs) | Isolated \(A\), blocked Hermitian Gram, `P0`, J0 unique-ρ Gram |
| [wasm/src/match_s.rs](../wasm/src/match_s.rs) | f64 \(Z=R+jX\), per-port or common \(z_0\), Kurokawa \(S\) and \(T\) |
| [wasm/src/lib.rs](../wasm/src/lib.rs) | `RadiatedPowerKernel` wasm-bindgen |
| [js/wasm/farfield-worker.js](../js/wasm/farfield-worker.js) | Same worker: far-field tiles **and** `run_prad` panels. Browser Dedicated Worker or Node `worker_threads` |
| [js/wasm/farfield-pool.js](../js/wasm/farfield-pool.js) | `runPradJob`, `mergeGrams`, `pradWorkerCount`, `stopFarfieldPool`; Node can pass `{Worker, wasmPath, workers}` |
| [tests/prad-gram.test.js](../tests/prad-gram.test.js) | Accuracy + panel/worker merge + J0 vs product |
| [tests/prad-bench.test.js](../tests/prad-bench.test.js) | Main-thread product/J0 and worker timings through \(64\times 64\) |
| [tests/matched-s.test.js](../tests/matched-s.test.js) | N=1 identity, 8×8 \(S_{ii}\)/worst \(|S_{ij}|\), \(T\) GEMV, conjugate phases |
| [js/phasedarray/matched.js](../js/phasedarray/matched.js) | GEMV \(w=Ta\), reflection ratio, conjugate \(F^\mathrm{emb}\) phases, \(n_\mu\) |

Rebuild after Rust changes: `./wasm/build.ps1` (simd + scalar into `js/wasm/`).

The display kernel in [wasm/src/kernel.rs](../wasm/src/kernel.rs) is unchanged.

---

## WASM API (`RadiatedPowerKernel`)

```text
new()
set_quadrature(n_mu, n_phi)
fill_isolated(x, y, frequency_scale, element_kind, element_n)
fill_isolated_range(..., sample0, sample_count)  # sample_count==0 → remainder
form_gram()
compute(...)                    # fill_isolated + form_gram
compute_j0(...)                 # axisymmetric planar fast path (Gauss-μ only)
take_re() / take_im()           # row-major N×N f32 Gram
form_matched_s(z_ref, x, y, x_nn, alpha, beta, aniso, z_common_re, z_common_im)
                                # R, optional jX(Δx,Δy), z0, S (z_ref≤0 → 50 Ω;
                                # Re(z_common)>0 → common z_c, else per-port)
take_z0() / take_z0_im()        # length-N f64; imag is 0 unless X_nn≠0 or Im(z_c)≠0
take_s_re() / take_s_im()       # row-major N×N f64
take_t_re() / take_t_im()       # row-major N×N f64, Kurokawa \(T\)
match_iterations() / match_residual()
n_samples()                     # full quadrature M = n_mu*n_phi
n_elements()
```

`form_gram` uses the **current** panel width of \(A\) (after a range fill, that is \(M/W\), not full \(M\)). Each panel’s \(P\) is already scaled by \(P_0/4\pi\); **summing panels** is the full Gram.

Geometry: `x`, `y` as `Float32Array`, wavelengths. Gram is f32; match/\(S\) promote to f64.

---

## J0 fast path (implemented)

For co-aligned planar elements with \(D=D(\mu)\) only, the φ integral is exact:

\[
(P_H)_{pq}=\frac{P_0}{2}\int_0^1 D(\mu)\,J_0\!\big(k\rho_{pq}\sqrt{1-\mu^2}\big)\,d\mu,
\]

\(\rho_{pq}=\sqrt{(x_p-x_q)^2+(y_p-y_q)^2}\). Same Gauss-μ as the product kernel (`HemisphereQuad.mu1d` / `w_mu`). Result is real symmetric; imag is stored as zeros.

`compute_j0(x, y, frequency_scale, element_kind, element_n)`:

- Requires a prior `set_quadrature`; uses \(n_\mu\) only, ignores \(n_\phi\).
- Does **not** fill \(A\). Same `take_re` / `take_im` buffers as `compute`.
- Unique-ρ cache keyed by \(\rho^2\) f32 bits (`dx*dx + dy*dy`), so lattice reflections \((3,4)\)/\((4,3)\) hit. \(\rho=0\) short-circuits to \(\sum c_i = P_0\) (isotropic Gauss-μ is exact).
- f32 \(J_0\) is fdlibm `j0f` in [wasm/src/bessel.rs](../wasm/src/bessel.rs). μ-sum in f64, store f32.
- HashMap uses a deterministic hasher (wasm32 has no `getrandom`).
- Cost \(O(U n_\mu + N^2)\) with \(U\ll N^2/2\) on rectangular lattices. At large \(N\) the pair loop that writes \(P_H\) dominates; time is almost independent of \(n_\mu\).

J0 is the exact φ integral; the product kernel approximates φ with trapezoid \(n_\phi\). They agree as \(n_\phi\to\infty\) at **fixed** \(n_\mu\). Tests compare J0(\(n_\mu\)) to product(\(n_\mu\), large \(n_\phi\)), not to a coarse φ grid.

Workers stay on the product / sample-panel path. Do not worker-split J0 unless φ-dependent patterns are dropped.

---

## Matched \(S\) (implemented, §6–8)

Runs on whatever Gram is already in the kernel (`compute_j0` or `compute` / panel merge). Does **not** fill \(A\). Positions `x`,`y` are used only to build optional \(X(\Delta x,\Delta y)\).

\[
R = 2\,Z_\mathrm{ref}\,P_H,\qquad
Z = R + jX(\Delta x,\Delta y),
\]

\[
X_{pq}=X_{nn}(d_\min/\rho)^\alpha
\cos\bigl(\beta(\rho/d_\min-1)\bigr)
\bigl(1+A\cos 2\varphi\bigr),
\quad \varphi=\mathrm{atan2}(\Delta y,\Delta x).
\]

\(X_{nn}=0\), empty `x`/`y`, or length mismatch skips \(X\) (real \(z_0\) path). \(\beta=0\), \(A=0\) is the older radial power law. Constants: \(Z_\mathrm{ref}=50\,\Omega\), \(\varepsilon_z=10^{-9}\,\Omega\), \(K_\mathrm{max}=200\), \(\tau=10^{-3}\,\Omega\).

- **\(X_{nn}=0\):** match on \(\Re(R)\); \(S=D^{-1/2}(R-D)\,\mathrm{solve}(R+D,\,D^{1/2})\). Real SPD / Hermitian Cholesky.
- **\(X_{nn}\neq 0\):** conjugate match \(z_0=Z_\mathrm{in}^*\); Kurokawa \(S=G^{-1}(Z-D^*)(Z+D)^{-1}G\) with \(G=\mathrm{diag}\sqrt{\Re(z_0)}\). Complex LU on \(Z+D\). After \(K_\mathrm{max}/4\), under-relax \(\beta=1/2\).

`form_matched_s(z_ref, x, y, x_nn, alpha, beta, aniso, z_common_re, z_common_im)` with `z_ref <= 0` or non-finite uses \(50\,\Omega\). Finite \(\Re(z_c)>0\) skips the per-port solver and sets every \(z_{0,p}=z_c\) (Kurokawa if \(\Im(z_c)\neq 0\) or \(X_{nn}\neq 0\)). Otherwise after a successful per-port match, \(\mathrm{diag}(S)\approx 0\). For J0 / axisymmetric planar elements and \(X_{nn}=0\) with real \(z_0\), \(R\), \(S\), and \(T\) are real.

One-port check: isolated \(P_H=P_0=1/2\) gives \(R=Z_\mathrm{ref}\). Per-port or common \(z_c=50\) both give \(z_0=50\), \(S=0\), \(T=1\).

The visualizer uses J0 + `form_matched_s` on geometry / element / `frequency_scale` / match style / \(z_c\) / \(X_{nn},\alpha,\beta,A\) changes, then \(w=T a\) in `create_farfield_vectors` when Coupling is Matched. Oscillation \(\beta\) is passed as \(\beta_{\mathrm{UI}}\times\texttt{frequency_scale}\). Conjugate steer sets \(a\) from \(\arg(T^T F^\mathrm{iso})\) at the commanded \((\theta,\phi)\). Illumination is applied after \(T\). Full \(F^\mathrm{emb}\) on the plot mesh is not built.

---

## Workers

Same pool as far-field visualization. Split the **sample axis** (like `ax2` rows), not element pairs:

\[
P_H \propto \sum_s a(s)\,a(s)^H
\quad\Rightarrow\quad
P = \sum_w A_{:,S_w} A_{:,S_w}^H.
\]

- `runPradJob({x,y,frequencyScale,elementKind,elementN,nMu,nPhi})` → `{re, im}`.
- Merge is a **sum** of \(N\times N\) tiles (far-field **concatenates** intensity). Tiles are added as they arrive.
- No `SharedArrayBuffer`.
- In-flight partial Grams are capped at ~512 MiB (`pradWorkerCount`): each tile is \(8N^2\) bytes (re+im f32). At \(N=4096\) that forces **4** workers even if 8 were started.
- Node tests pass `Worker` from `node:worker_threads` and `wasmPath` so the worker can `readFile` the `.wasm` (browser workers `fetch` it).

φ-dependent \(D(\mu,\phi)\) would only change `amp[s]`; the split stays valid.

---

## Tests

```powershell
cargo test --manifest-path wasm/Cargo.toml
node --test tests/prad-gram.test.js
node --test tests/prad-bench.test.js
node --test tests/matched-s.test.js
```

Rust covers \(N=1\) power, coincident/far pairs, Hermitian, quadrature convergence, naive Gram, **panel sum = full Gram**, **J0 vs product**, unique-ρ collapse, and **matched \(S\)/\(T\)**: N=1 \(S=0\), \(T=1\), far-pair, residual \(<\tau\), symmetry, pairwise \(jX(\Delta x,\Delta y)\) (including irregular / sunflower-like geometry, oscillation \(\beta\), and location \(A\)). JS prints 8×8 \(S_{ii}\) / worst \(|S_{ij}|\), checks \(w=Ta\), conjugate vs geometric phases, a three-element \(X_{nn}\) case, and \(\beta\)/\(A\) overlays.

---

## Benchmark (measured)

Rectangular grids, \(0.5\lambda\) spacing, isotropic hemisphere, `frequency_scale=1`. \(M=n_\mu n_\phi\). SIMD unless noted. Workers: 8 started; 64×64 uses 4 because of the Gram-memory cap.

Product kernel (fill \(A\) is cheap; **Gram is the cost**):

| array | \(N\) | \(M\) | main SIMD | workers SIMD | speedup |
| --- | --- | --- | --- | --- | --- |
| 8×8 | 64 | 2048 | 2.3 ms | 0.9 ms | — |
| 16×16 | 256 | 2048 | 42 ms | 9.2 ms | ~5× |
| 32×32 | 1024 | 2048 | 0.95 s | 142 ms | ~7× |
| 64×64 | 4096 | 512 | 5.5 s | 2.1 s | ~2.6× |
| 64×64 | 4096 | 2048 | **56 s** | **12.8 s** | ~4× |

Scalar product Gram is ~1.5–2× slower than SIMD. Workers change wall-clock, not \(O(N^2 M)\).

J0 fast path, main thread (SIMD ≈ scalar; \(n_\phi\) unused):

| array | \(N\) | \(n_\mu\) | J0 | vs product \(M=2048\) |
| --- | --- | --- | --- | --- |
| 8×8 | 64 | 16–32 | ~0.1 ms | — |
| 16×16 | 256 | 16–32 | ~1–2 ms | ~20× |
| 32×32 | 1024 | 16–32 | ~63 ms | ~15× |
| 64×64 | 4096 | 16 or 32 | **~2.5 s** | ~22× |

J0 time is almost independent of \(n_\mu\): unique-\(\rho\) hits, then the \(N^2\) write of \(P_H\) dominates.

Suggested default orders vs electrical size were **not** used in the bench; 32/64 is a moderate grid, not sized to a 64×64 aperture (\(D_\lambda\sim 32\) would want larger \(M\)).

---

## Not done (continue here)

1. **\(F^\mathrm{emb}\) on a quadrature grid** — pattern mixing \(F^\mathrm{emb}=T^T F^\mathrm{iso}\) for power-balance checks. The GUI uses \(w=Ta\) into the existing AF kernel instead. `compute_j0` never builds \(A\); use `fill_isolated` if fields on the hemisphere are needed.
2. **φ-dependent / polarized elements** — product grid is the right layout (`N×M`, later \(N\times 2M\)). `fill_amp` currently uses \(D(\mu)\) only. Workers stay on this path.
3. **f64 Gram** if later \(P_\mathrm{loss}\) residuals are too large; hot kernel is f32 to match the far-field crate. Match/\(S\)/\(T\) already run in f64.

---

## Design constraints for a follow-on agent

- Do not compute \(P_H\) on the plot mesh; use the hemisphere product quadrature (or J0, which uses the same Gauss-μ).
- Do not assume the visualizer’s isotropic apply (it does not zero the back hemisphere on plots). The Gram uses hemispheric \(D=2\).
- Preserve the sample-panel / \(A A^H\) path for general / φ-dependent / polarized patterns. J0 is optional and planar-axisymmetric only.
- Do not auto-switch `compute()` or `runPradJob` to J0. Do not drop the worker product path.
- For method §9 mixing, `compute_j0` never builds \(A\); keep using `fill_isolated` if \(F^\mathrm{iso}\) is needed on the quadrature grid.
- Rebuild both wasm artifacts; node tests load `js/wasm/simd` and `js/wasm/scalar`.
