# Radiated-power Gram (WASM status)

Handoff for [approximate_matched_basis.md](approximate_matched_basis.md). **§5–8 are implemented**: isolated fields → Hermitian \(P_H\) (φ-dependent product kernel and axisymmetric \(J_0\) fast path), then \(R = 2 Z_\mathrm{ref} P_H\), simultaneous real \(z_0\) match, and power-wave \(S\). \(T\) and \(F^\mathrm{emb}\) are not started. There is no visualizer UI.

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
| [wasm/src/match_s.rs](../wasm/src/match_s.rs) | f64 \(R\), real \(z_0\) match, power-wave \(S\) (custom Cholesky) |
| [wasm/src/lib.rs](../wasm/src/lib.rs) | `RadiatedPowerKernel` wasm-bindgen |
| [js/wasm/farfield-worker.js](../js/wasm/farfield-worker.js) | Same worker: far-field tiles **and** `run_prad` panels. Browser Dedicated Worker or Node `worker_threads` |
| [js/wasm/farfield-pool.js](../js/wasm/farfield-pool.js) | `runPradJob`, `mergeGrams`, `pradWorkerCount`, `stopFarfieldPool`; Node can pass `{Worker, wasmPath, workers}` |
| [tests/prad-gram.test.js](../tests/prad-gram.test.js) | Accuracy + panel/worker merge + J0 vs product |
| [tests/prad-bench.test.js](../tests/prad-bench.test.js) | Main-thread product/J0 and worker timings through \(64\times 64\) |
| [tests/matched-s.test.js](../tests/matched-s.test.js) | N=1 identity + 8×8 isotropic \(\lambda/2\) J0 \(S\) printout |

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
form_matched_s(z_ref)           # R, z0 match, S on the current Gram (z_ref≤0 → 50 Ω)
take_z0()                       # length-N f64
take_s_re() / take_s_im()       # row-major N×N f64
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

Runs on whatever Gram is already in the kernel (`compute_j0` or `compute` / panel merge). Does **not** fill \(A\).

\[
R = 2\,Z_\mathrm{ref}\,P_H,\qquad
z_0 \text{ from } \Re(R),\qquad
S = D^{-1/2}(R-D)\,\mathrm{solve}(R+D,\,D^{1/2}).
\]

Constants: \(Z_\mathrm{ref}=50\,\Omega\), \(\varepsilon_z=10^{-9}\,\Omega\), \(K_\mathrm{max}=200\), \(\tau=10^{-3}\,\Omega\). Match is the method’s fixed point on \(\Re(R)\); \(S\) uses full Hermitian \(R\). Custom f64 Cholesky (real SPD for the match, Hermitian PD for \(R+D\)); a tiny pivot is floored rather than panicking.

`form_matched_s(z_ref)` with `z_ref <= 0` or non-finite uses \(50\,\Omega\). After a successful match, \(\mathrm{diag}(S)\approx 0\). For J0 / axisymmetric planar elements \(R\) and \(S\) are real.

One-port check: isolated \(P_H=P_0=1/2\) gives \(R=Z_\mathrm{ref}\), \(z_0=50\), \(S=0\).

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

Rust covers \(N=1\) power, coincident/far pairs, Hermitian, quadrature convergence, naive Gram, **panel sum = full Gram**, **J0 vs product** (same \(n_\mu\), large \(n_\phi\)) plus unique-ρ collapse on an 8×8 lattice, and **matched \(S\)**: N=1 \(S=0\), far-pair weak coupling, match residual \(<\tau\), \(S\) symmetry. JS repeats the Gram checks on simd/scalar WASM and prints 8×8 isotropic \(\lambda/2\) \(S_{ii}\) plus worst-case \(|S_{ij}|\) from J0.

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

1. **Visualizer UI** — geometry / element / frequency_scale → `runPradJob` or main-thread kernel, then `form_matched_s`. Pool already starts with far-field WASM init.
2. **Method §9–11** — \(T\) and \(F^\mathrm{emb}=T^T F^\mathrm{iso}\), power-balance checks. Need to **keep** \(A\) or \(F^\mathrm{iso}\) on the quadrature grid for mixing (`compute_j0` never builds \(A\)).
3. **φ-dependent / polarized elements** — product grid is the right layout (`N×M`, later \(N\times 2M\)). `fill_amp` currently uses \(D(\mu)\) only. Workers stay on this path.
4. **f64 Gram** if later \(P_\mathrm{loss}\) residuals are too large; hot kernel is f32 to match the far-field crate. Match/\(S\) already run in f64.

---

## Design constraints for a follow-on agent

- Do not compute \(P_H\) on the plot mesh; use the hemisphere product quadrature (or J0, which uses the same Gauss-μ).
- Do not assume the visualizer’s isotropic apply (it does not zero the back hemisphere on plots). The Gram uses hemispheric \(D=2\).
- Preserve the sample-panel / \(A A^H\) path for general / φ-dependent / polarized patterns. J0 is optional and planar-axisymmetric only.
- Do not auto-switch `compute()` or `runPradJob` to J0. Do not drop the worker product path.
- For method §9 mixing, `compute_j0` never builds \(A\); keep using `fill_isolated` if \(F^\mathrm{iso}\) is needed on the quadrature grid.
- Rebuild both wasm artifacts; node tests load `js/wasm/simd` and `js/wasm/scalar`.
