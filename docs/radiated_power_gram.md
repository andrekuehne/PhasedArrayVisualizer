# Radiated-power Gram (WASM status)

Handoff for [approximate_matched_basis.md](approximate_matched_basis.md). **Only §5 is implemented:** isolated fields from array geometry plus the existing element factor, then the Hermitian radiated-power matrix \(P_H\). Match, \(S\), \(T\), and \(F^\mathrm{emb}\) are not started. There is no visualizer UI.

Intent: keep a **φ-dependent** product quadrature as the general kernel (needed later for rotated / polarized patterns). Axisymmetric \(J_0\) collapse is documented as a fast path, not built.

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
| [wasm/src/quadrature.rs](../wasm/src/quadrature.rs) | Gauss–Legendre + `HemisphereQuad` |
| [wasm/src/prad.rs](../wasm/src/prad.rs) | Isolated \(A\), blocked Hermitian Gram, `P0` |
| [wasm/src/lib.rs](../wasm/src/lib.rs) | `RadiatedPowerKernel` wasm-bindgen |
| [js/wasm/farfield-worker.js](../js/wasm/farfield-worker.js) | Same worker: far-field tiles **and** `run_prad` panels. Browser Dedicated Worker or Node `worker_threads` |
| [js/wasm/farfield-pool.js](../js/wasm/farfield-pool.js) | `runPradJob`, `mergeGrams`, `pradWorkerCount`, `stopFarfieldPool`; Node can pass `{Worker, wasmPath, workers}` |
| [tests/prad-gram.test.js](../tests/prad-gram.test.js) | Accuracy + panel/worker merge |
| [tests/prad-bench.test.js](../tests/prad-bench.test.js) | Main-thread and worker timings through \(64\times 64\) |

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
take_re() / take_im()           # row-major N×N f32
n_samples()                     # full quadrature M = n_mu*n_phi
n_elements()
```

`form_gram` uses the **current** panel width of \(A\) (after a range fill, that is \(M/W\), not full \(M\)). Each panel’s \(P\) is already scaled by \(P_0/4\pi\); **summing panels** is the full Gram.

Geometry: `x`, `y` as `Float32Array`, wavelengths. f32 throughout.

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
```

Rust covers \(N=1\) power, coincident/far pairs, Hermitian, quadrature convergence, naive Gram, **panel sum = full Gram**. JS repeats that on simd/scalar WASM and checks sequential + real worker merge vs main thread.

---

## Benchmark (measured)

Rectangular grids, \(0.5\lambda\) spacing, isotropic hemisphere, `frequency_scale=1`. \(M=n_\mu n_\phi\). SIMD unless noted. Workers: 8 started; 64×64 uses 4 because of the Gram-memory cap.

| array | \(N\) | \(M\) | main SIMD | workers SIMD | speedup |
| --- | --- | --- | --- | --- | --- |
| 8×8 | 64 | 2048 | 2.4 ms | 1.0 ms | — |
| 16×16 | 256 | 2048 | 54 ms | 10 ms | ~5× |
| 32×32 | 1024 | 2048 | 0.95 s | 181 ms | ~5× |
| 64×64 | 4096 | 512 | 5.5 s | 2.0 s | ~2.8× |
| 64×64 | 4096 | 2048 | **63 s** | **15.3 s** | ~4× |

Fill \(A\) is cheap; **Gram is the cost**. Scalar is ~1.5–2× slower than SIMD on the Gram. Workers change wall-clock, not \(O(N^2 M)\).

Suggested default orders vs electrical size were **not** used in the bench; 32/64 is a moderate grid, not sized to a 64×64 aperture (\(D_\lambda\sim 32\) would want larger \(M\)).

---

## Not done (continue here)

1. **Visualizer UI** — geometry / element / frequency_scale → `runPradJob` or main-thread kernel. Pool already starts with far-field WASM init.
2. **Method §6–11** — \(R=2 Z_\mathrm{ref} P_H\), real \(z_0\) match, \(S\) and \(T\), \(F^\mathrm{emb}=T^T F^\mathrm{iso}\), power-balance checks. \(O(N^3)\) after the Gram; need to **keep** \(A\) or \(F^\mathrm{iso}\) on the quadrature grid for mixing (today \(A\) is discarded after `form_gram` except inside the kernel).
3. **φ-dependent / polarized elements** — product grid is the right layout (`N×M`, later \(N\times 2M\)). `fill_amp` currently uses \(D(\mu)\) only.
4. **Optional \(J_0\) fast path** (axisymmetric \(D(\mu)\), planar arrays only):

   \[
   (P_H)_{pq}=\frac{P_0}{2}\int_0^1 D(\mu)\,J_0\!\big(k\rho_{pq}\sqrt{1-\mu^2}\big)\,d\mu.
   \]

   Real, depends only on distance; unique-\(\rho\) cache on rectangular lattices. Keep Gauss-μ, φ as the general kernel and as a check. Do not replace the worker path if φ-dependent patterns remain a requirement.
5. **f64 Gram** if later \(P_\mathrm{loss}\) residuals are too large; hot kernel is f32 to match the far-field crate.

---

## Design constraints for a follow-on agent

- Do not compute \(P_H\) on the plot mesh; use the hemisphere product quadrature.
- Do not assume the visualizer’s isotropic apply (it does not zero the back hemisphere on plots). The Gram uses hemispheric \(D=2\).
- Preserve a sample-panel / \(A A^H\) path for general patterns even if \(J_0\) is added.
- Rebuild both wasm artifacts; node tests load `js/wasm/simd` and `js/wasm/scalar`.
