# Green-function coupling: sequential work packages

Handoff plan for finite-array Green coupling. **No Floquet / infinite-array path.**
\(F^\mathrm{iso}\) and complex \(Z\) come from the same environment Green function;
existing \(S\), \(T\), and array-factor machinery stay. First physics: **horizontal
short dipole over infinite PEC** (finite \(\ell\), not a point Hertzian). Later:
patch over grounded slab, same interface.

**Green PEC (implemented):** Element dropdown option **PEC dipole**. Matching
(Isolated / Per-port / Common Z0 / Propagation) is hidden in that mode and remains
a toy for isotropic / \(\cos^n\). \(Z\) from `form_green_pec_dipole`, pattern from
`apply_green_pec_pattern`. Controls: \(h,\ell,a\), real \(z_c\), series Self X
(added to \(\mathrm{diag}(Z)\)). Match sets \(z_c=\Re Z_{11}\) and
\(X_\mathrm{self}=-\Im Z_{11}\) of the isolated self kernel. Auto-isolate Green LU
at \(N>1024\); 32×32 stays matched.

**Green slab (implemented):** Element dropdown option **Slab dipole**. Same matching
hide / Match / auto-isolate as PEC. Extra controls \(\varepsilon_r,h_\mathrm{sub},\tan\delta\).
\(Z\) from `form_green_slab_dipole`, pattern from `apply_green_slab_pattern`.
Power rows: isolated Radiated / Surface wave / Substrate loss / Closure at 1 A.
EIRP stays \(D\times P_\mathrm{accepted}\). Surface-wave power is not in the AF.

**Agent rule:** complete **one** work package per session. Do not start the next WP.
End with the handoff block so the following agent can load this file plus the
listed paths only.

Related: [radiated_power_matching_audit.md](radiated_power_matching_audit.md)
(Stage 3, no scan-angle multipliers), [radiated_power_gram.md](radiated_power_gram.md)
(J0 unique-\(\rho\), \(S/T\)), [approximate_matched_basis.md](approximate_matched_basis.md)
(§8–9 mixing).

---

## 1. Goal

For a chosen radiator + infinite laterally invariant environment, at one frequency:

| Output | Source |
| --- | --- |
| \(Z\in\mathbb{C}^{N\times N}\) | reaction \(\langle J_p,G,J_q\rangle\) (full spectrum of \(G\)) |
| \(F^\mathrm{iso}(\hat r)\) | space-wave / saddle-point piece of the same \(G\cdot J\) |
| \(S\), \(T\), \(z_0\) | existing `MatchedS` from that \(Z\) (not from the Gram) |
| Array far field | existing AF of \(w=Ta\), times \(\lvert F^\mathrm{iso}\rvert^2\) (vector \(E_\theta,E_\phi\) later) |

Green mode **replaces** Gram-as-\(R\), phenomenological \(X(\rho)\), and independent
\(\cos^n\). Isolated / Per-port / Common Z0 / Propagation stay as labelled toys.

---

## 2. Non-goals (all WPs)

- No Floquet unit-cell / infinite-array scan solver.
- No rewrite of `wasm/src/kernel.rs` accumulate, worker AF tiling, or Kurokawa \(S,T\).
- Do not set \(\Re(Z)\) from the radiated-power Gram in Green mode.
- Do not multiply \(Z\) or \(S\) by \(1/\cos\theta\).
- Do not treat a complex \(z_0\) as a matching network.
- No finite ground/substrate edge diffraction.
- No full-wave MoM of patch metal (equivalent slots only, later WP).

---

## 3. Architecture (frozen)

Two callbacks, one environment \(G\):

```text
F_iso(θ, φ; radiator, env, f)  →  |F|²   (WP4; Eθ,Eφ optional later)
Z_pair(Δx, Δy; radiator, env, f) →  complex, including self
```

Then:

1. Unique-lag table of \(Z_pair\) (same idea as J0 unique-\(\rho\)).
2. Scatter into dense \(N\times N\) \(Z\).
3. `MatchedS::from_z` (WP2) → \(S,T\).
4. GUI: \(w=Ta\) as today; pattern multiply uses \(F^\mathrm{iso}\), not `PATTERN_COS_N`.

**Do not** keep the Gram as \(\Re(Z)\). PEC dipole: hemispheric Gram of \(F^\mathrm{iso}\)
should **agree** with \(\Re Z\) (test). Slab/patch later: extra \(\Re Z\) is surface-wave
/ dielectric and must not be plotted as space radiation.

Positions stay in wavelengths at \(f_0\); electrical \(k=2\pi\cdot\texttt{frequency_scale}\).

---

## 4. Fast Rust: feasibility and 32×32

Yes. Dipole-over-PEC is closed form (direct spherical wave + image). It belongs in
Rust/WASM next to J0, not in JS.

Cost model (identical co-aligned elements, planar lattice):

| Step | Complexity | 32×32 (\(N=1024\)) |
| --- | --- | --- |
| Unique lags | \(O(U)\) kernel evals | Rectangular: \(U\sim 32\times 32=1024\) distinct \((\Delta i,\Delta j)\), not \(N^2/2\) |
| Closed-form \(Z_pair\) | tens of flops + `exp` | microseconds |
| Scatter into \(Z\) | \(O(N^2)\) | same class as J0 write (**~63 ms** J0 Gram in `radiated_power_gram.md`) |
| \(S,T\) from \(Z\) | \(O(N^3)\) complex LU | **this is the limiter**, already in `form_matched_s` |
| \(F^\mathrm{iso}\) on plot grid | \(O(M)\) once | negligible vs AF |

J0 already showed: at large \(N\) the pair loop that **writes** the matrix dominates,
not the special function. PEC dipole is cheaper per lag than \(J_0\) quadrature, so
**kernel fill at 32×32 should be ≲ J0 (~tens of ms), likely dominated by the \(N^2\)
store.** Memory: one complex f64 \(N\times N\) is 16 MiB; \(Z+S+T\) ~48 MiB. Fine.

**32×32 is therefore “same toughness as today’s matched path,” not harder.**
The Green function is not the issue. Practical constraints that already exist:

- `MATCHED_AUTO_ISOLATE_N = 512` forces Isolated on rebuild for \(N>512\). 32×32
  will not stay in a coupling mode until a later WP raises or special-cases Green.
- Complex LU at \(N=1024\) is \(\sim N^3\) and will likely be **~0.5–few seconds**
  in WASM (not measured here). Acceptable on geometry/frequency change, not per
  steer. Steer only uses \(T\) GEMV + AF, as today.
- 64×64 (\(N=4096\)): \(N^2\) write ~2.5 s (J0 analog); LU \(\sim 64\times\) 32×32
  → tens of seconds and ~0.75 GiB for \(Z,S,T\). Out of interactive scope unless
  a later WP adds a truncated-lag or matrix-free path. **Do not promise 64×64
  Green matching in WP1–5.**

Workers: **do not** worker-split PEC \(Z\). The kernel is too cheap; J0 stayed
single-thread for the same reason. Revisit workers only if a later Sommerfeld
unique-lag table is slow.

Implementation rules for speed:

- f64 kernel + \(Z\) (match path is already f64). f32 pattern/AF as today.
- Cache key \((\Delta x,\Delta y)\) (dipole is **not** isotropic: \(Z(\Delta x,\Delta y)\neq Z(\Delta y,\Delta x)\)).
  Use lattice index keys on rectangular grids; fall back to bit-keyed floats like J0
  for sunflower/irregular.
- Evaluate each unique lag **once**, then write both triangles (reciprocity \(Z=Z^T\)).
- No `HashMap` on the \(N^2\) write hot loop; HashMap only for unique-lag lookup
  (deterministic hasher, same as J0).

---

## 5. WP1 radiator (first physics)

**Finite short dipole**, not a point Hertzian. Uniform current \(I\) on length
\(\ell\ll\lambda\), equivalent cylindrical radius \(a\), \(+\hat x\), center at
height \(h\) over infinite PEC in \(z=0\). All elements co-aligned (WP1–5).

A point Hertzian has no finite self \(Z\). Mutual fields may use the equivalent
moment \(I\ell\); self free-space \(R\) and \(X\) use the short-dipole formulae
below. Do not evaluate \(1/R^3\) at \(R=0\).

### Coordinates and constants

- Positions \((x,y)\) and \((h,\ell,a)\) are in **wavelengths at \(f_0\)**.
- \(k=2\pi\cdot\texttt{freq\_scale}\), \(\eta_0=120\pi\).
- Time convention **\(e^{j\omega t}\)**, outgoing \(e^{-jkR}/R\), \(Z=R+jX\).
- \(\theta=0\) is \(+\hat z\) (boresight); \(\phi=0\) is \(+\hat x\).
- Defaults: \(h=0.25\), \(\ell=0.1\), \(a=0.001\), \(\texttt{freq\_scale}=1\).

Image: \(-\hat x\) current at \(z=-h\) (tangential \(E\) reversed on the PEC).

### Mutual \(Z\) (including the image of a pair)

For two centers separated by \(\mathbf{R}=( \Delta x,\Delta y,\Delta z )\) with
\(R=\lvert\mathbf{R}\rvert>0\), take the \(x\)-directed Hertzian field of moment
\(I\ell\) (Balanis infinitesimal dipole \(E_R,E_\theta\), then Cartesian \(E_x\)).
Induced EMF at unit current:

\[
Z(\mathbf{R})=-\frac{E_x(\mathbf{R})\,\ell}{I}.
\]

A pair at the same height is the direct term \(\mathbf{R}=(\Delta x,\Delta y,0)\)
plus the image term of a **reversed** moment at \(\mathbf{R}=(\Delta x,\Delta y,2h)\):

\[
Z_{pq}=Z(\Delta x,\Delta y,0)-Z(\Delta x,\Delta y,2h).
\]

(The minus is the reversed image current.) Center-sampled EMF is the WP1
approximation; valid when \(\ell\) is well below element spacing (default
\(\ell=0.1\), typical \(d=0.5\)).

### Self \(Z\)

\[
Z_{\mathrm{self}}=Z_{\mathrm{fs}}(\ell,a)+Z_{\mathrm{image}},
\qquad
Z_{\mathrm{image}}=-Z(0,0,2h).
\]

\(Z_{\mathrm{image}}\) is finite (\(R=2h>0\)). Free-space short-dipole (uniform current):

\[
R_{\mathrm{fs}}=\frac{\eta_0(k\ell)^2}{6\pi},\qquad
X_{\mathrm{fs}}=-\frac{\eta_0}{\pi\,k\ell}\left(\ln\frac{\ell}{2a}-1\right).
\]

\(Z_{\mathrm{fs}}=R_{\mathrm{fs}}+j X_{\mathrm{fs}}\). Require \(a>0\), \(\ell>2a\),
\(h>0\). \(\mathrm{Re}Z_{\mathrm{self}}>0\) and \(\mathrm{Im}Z_{\mathrm{self}}\) finite.

### Isolated pattern

Far-field angular factor of the **same** unit-current short dipole plus image
(space-wave only). Return \((E_\theta,E_\phi)\) and
\(\lvert F\rvert^2=\lvert E_\theta\rvert^2+\lvert E_\phi\rvert^2\). Scale: \(I=1\),
the usual \(e^{-jkr}/r\) stripped so \(E_{\theta,\phi}\) are the \(r\)-independent
factors in V (per the same \(\eta_0,I\ell\) constants as \(Z\)).

**Zero the back hemisphere** \(\theta>\pi/2\) and UV-invisible directions, same as
`PATTERN_COS_N`. No \(\cos^n\) input. Do not renormalize to \(\int D\,d\Omega=4\pi\)
in WP1; later Gram \(\leftrightarrow\operatorname{Re}Z\) needs this raw current basis.

Horizontal-\(x\) over PEC: image interference is a \(\sin(kh\cos\theta)\)-type
factor on the parallel field (test at \(h=1/4\)).

---

## 6. Work package index

| WP | Title | Next agent reads |
| --- | --- | --- |
| 1 | Short-dipole PEC kernel (Rust): \(Z_pair\), \(F^\mathrm{iso}\), tests | this file §§5 and 7 |
| 2 | `MatchedS::from_z` + wasm fill of dense \(Z\) | §8 + WP1 API |
| 3 | Unique-lag array builder + benches through 32×32 | §9 |
| 4 | Plug \(F^\mathrm{iso}\) into AF pattern multiply | §10 |
| 5 | GUI Green/PEC mode, size safeguard, CHANGE.md | §11 |
| 6 | Spectral PEC check (optional; gate for slab) | §12 |
| 7 | Grounded-slab \(G\) (Sommerfeld) | §13 |
| 8 | Two-slot patch on slab | §14 |

---

## 7. WP1 — Short-dipole PEC kernel (Rust only)

**Owner intent:** implement §5 as closed-form Rust functions. No GUI, no \(S,T\),
no \(N\times N\) scatter. **Start prompt for a fresh agent:** *WP1 only, §§1–5
and §7 of `docs/green_function_plan.md`. Finite short dipole over PEC, not a
point Hertzian.*

### Read

- This file §§1–5 (the formulae in §5 are the spec; do not substitute another
  self-term).
- `wasm/src/lib.rs` module list (add `mod green;` only).
- `wasm/src/element.rs` for \(\theta,\phi\) / back-hemisphere convention (do not change it).

### Do

1. Add **`wasm/src/green.rs`**. No wasm-bindgen. `pub` functions + `#[cfg(test)]`.
   `lib.rs` gets `mod green;` and nothing else.
2. Constants: `ETA0 = 120.0 * PI`, defaults \(h=0.25\), \(\ell=0.1\), \(a=0.001\).
3. `z_pair_pec_dipole(dx, dy, h, ell, a, freq_scale) -> (re, im)` in ohms, §5
   mutual + self when `dx==0 && dy==0`.
4. `f_iso_pec_dipole(theta, phi, h, ell, freq_scale) -> (e_th_re, e_th_im, e_ph_re, e_ph_im)`
   and a helper that returns \(\lvert F\rvert^2\). Zero \(\theta>\pi/2\).
5. Reciprocity \(Z(\Delta x,\Delta y)=Z(-\Delta x,-\Delta y)\); \(\mathrm{Re}Z_{\mathrm{self}}>0\);
   \(\mathrm{Im}Z_{\mathrm{self}}\) finite; no NaNs.
6. Tests (Rust only):
   - coincident pair uses the self formula, not \(1/R^3\);
   - \(h\to 0\) drives \(\lvert F\rvert^2\to 0\) for \(\theta<\pi/2\) (shorted parallel \(E\));
   - large \(\rho\): mutual phase tracks \(e^{-jk\rho}\);
   - \(h=0.25\), \(\phi=0\) vs \(\phi=\pi/2\): E/H-plane shapes of an \(x\)-dipole over PEC;
   - \(Z_{fs}\) radiation resistance matches \(\eta_0(k\ell)^2/(6\pi)\) when the image is
     omitted in a unit test helper.

### Do not

- Do not touch `match_s.rs`, `prad.rs` Gram, JS, or `index.html`.
- Do not implement Sommerfeld, slab, or a point-Hertzian self.
- Do not fill \(N\times N\) matrices yet.
- Do not power-normalize \(F^\mathrm{iso}\) to 4π.

### Done when

- `cargo test --manifest-path wasm/Cargo.toml` covers the kernel.
- File-top comment restates the §5 signatures (WP2 contract).

### Handoff

```text
WP1 done. Next: WP2.
API: z_pair_pec_dipole(dx, dy, h, ell, a, freq_scale) -> (re, im)
     f_iso_pec_dipole(theta, phi, h, ell, freq_scale) -> (Eθ, Eφ)
File: wasm/src/green.rs  (mod green in lib.rs only)
Do not change GUI. Read match_s.rs from_z_common.
```

### WP1 results (implemented)

- Files: `wasm/src/green.rs` (kernel + tests); `wasm/src/lib.rs` is `mod green;` only (no bindgen, no `pub use`).
- Signatures (f64, ohms / stripped \(E_{\theta,\phi}\)):

```text
ETA0 = 120 * PI
DEFAULT_H = 0.25, DEFAULT_ELL = 0.1, DEFAULT_A = 0.001
z_pair_pec_dipole(dx, dy, h, ell, a, freq_scale) -> (re, im)
f_iso_pec_dipole(theta, phi, h, ell, freq_scale) -> (e_th_re, e_th_im, e_ph_re, e_ph_im)
f_iso_pec_dipole_power(...) -> |F|²
```

- Physics for WP2: coincident pair uses \(Z_\mathrm{fs}(\ell,a)-Z(0,0,2h)\), not \(1/R^3\); mutual is Hertzian \(E_x\) plus reversed image. \(F^\mathrm{iso}\) is real with a \(\sin(kh\cos\theta)\) image factor; back hemisphere \(\theta>\pi/2\) is zero. Do not call the Gram. Do not power-normalize to \(4\pi\).
- Tests: `cargo test --manifest-path wasm/Cargo.toml` (8 `green::tests::*`, 91 total wasm tests).
- Next: WP2. `from_z` is the only match entry for Green mode. Do not change GUI.

```text
WP1 done. Next: WP2.
API: z_pair_pec_dipole(dx, dy, h, ell, a, freq_scale) -> (re, im)
     f_iso_pec_dipole(theta, phi, h, ell, freq_scale) -> (Eθ, Eφ)
File: wasm/src/green.rs  (mod green in lib.rs only)
Do not change GUI. Read match_s.rs from_z_common.
```

---

## 8. WP2 — `from_z` and dense \(Z\) wiring

**Owner intent:** any complex \(Z\) (here PEC dipole, naïve \(O(N^2)\) pair loop is
OK) produces \(S,T\) through existing Kurokawa code.

### Read

- WP1 handoff + `wasm/src/match_s.rs` (`from_gram_coupled`, `from_z_common`).
- `wasm/src/lib.rs` `form_matched_s*` bindgen pattern.
- `docs/approximate_matched_basis.md` §8.

### Do

1. Add `MatchedS::from_z(z_re, z_im, n, z_ref, z_common_re)` that **does not**
   scale a Gram. Reuse `from_z_common`. Default common real \(z_c=Z_\mathrm{ref}\)
   if invalid.
2. Add `RadiatedPowerKernel::form_green_pec_dipole(...)` (name flexible) that
   builds dense \(Z\) with a double loop over elements calling `z_pair_*`, then
   `from_z`. Slow is acceptable; WP3 makes it fast.
3. Tests: \(N=1\) \(\to\) \(S=0\), \(T=1\) when \(Z_{11}=z_c\) (tune \(a,\ell,h\) or
   \(z_c=\Re Z_{11}\) in the test); \(N=2\) reciprocity; `take_s_*` / `take_t_*`
   same buffers as today.
4. **Do not** call `compute_j0` / Gram in this path.

### Do not

- Unique-lag cache (WP3).
- Pattern / GUI (WP4–5).

### Done when

- Node or Rust tests prove \(S,T\) from Green \(Z\) without \(P_H\).
- Rebuild note: `./wasm/build.ps1` if bindgen was added.

### Handoff

```text
WP2 done. Next: WP3.
from_z is the only match entry for Green mode.
Naïve N² z_pair fill exists; replace with unique-lag in WP3.
```

### WP2 results (implemented)

- Files: `wasm/src/match_s.rs` (`MatchedS::from_z`); `wasm/src/prad.rs` (`form_green_pec_dipole` naïve \(N^2\)); `wasm/src/lib.rs` bindgen of the same name. `take_s_*` / `take_t_*` / `take_z_*` unchanged.
- Signatures:

```text
MatchedS::from_z(z_re, z_im, n, z_ref, z_common_re)
  // Z already ohms; no Gram scale. Invalid z_c → z_ref. Always common real z_c.
PradState / RadiatedPowerKernel::form_green_pec_dipole(
  x, y, frequency_scale, h, ell, a, z_ref, z_common_re)
  // Z_pq at p*n+q, dx=x[p]-x[q], dy=y[p]-y[q]; diagonal is WP1 self.
```

- Physics for WP3: Green \(Z\) never uses \(P_H\). \(N=1\) \(S=0,T=1\) only for **real** \(Z_{11}=z_c\); dipole \(X_\mathrm{fs}\) leaves \(|S_{11}|\). Rebuild: `./wasm/build.ps1`.
- Tests: `cargo test --manifest-path wasm/Cargo.toml` (`from_z_n1_real_is_open`, `from_z_invalid_zc_clamps_to_zref`, `green_n1_matches_z_pair_self_no_gram`, `green_n2_reciprocal_matches_z_pair`; 95 wasm tests).
- Next: WP3 unique-lag fill. Do not change GUI.

```text
WP2 done. Next: WP3.
from_z is the only match entry for Green mode.
Naïve N² z_pair fill exists; replace with unique-lag in WP3.
```

---

## 9. WP3 — Unique-lag fill + performance

**Owner intent:** 8×8 and 16×16 interactive; 32×32 kernel fill in the same ballpark
as J0 (~tens of ms), LU extra and measured.

### Read

- `wasm/src/prad.rs` J0 unique-\(\rho\) cache.
- `tests/prad-bench.test.js` timing style.
- This file §4.

### Do

1. Replace naïve pair `z_pair` with unique \((\Delta x,\Delta y)\) (lattice keys on
   rect grids; float-bit keys otherwise).
2. Reciprocal triangle write; diagonal from self kernel once.
3. Bench 8×8, 16×16, 32×32: (a) unique-lag \(Z\) fill, (b) `from_z` LU, (c) total.
   Record in this doc or `CHANGE.md` only if WP5 is the same agent — **prefer a
   table in this file under “WP3 results”** so WP5 can cite it.
4. Assert unique-lag \(Z\) matches naïve fill on 4×4.

### Do not

- Raise `MATCHED_AUTO_ISOLATE_N` (WP5 policy).
- Workers.
- 64×64 as a target.

### Done when

- 32×32 \(Z\) fill is clearly not the LU; numbers written in the handoff.
- Tests: naïve vs unique-lag.

### Handoff

```text
WP3 done. Next: WP4.
Z fill ms @ 8/16/32: ...
from_z LU ms @ 8/16/32: ...
32×32 conclusion: ...
```

### WP3 results (implemented)

- Files: `wasm/src/prad.rs` unique-lag fill (`fill_green_pec_dipole_z`; rect \((\lvert\Delta i\rvert,\lvert\Delta j\rvert)\) table, else bit-keyed HashMap); `form_from_z`; `form_green_pec_dipole` still fill+match. Bindgen of the split APIs in `wasm/src/lib.rs`. Naive loop is `#[cfg(test)]` only.
- 4×4 rect unique-lag \(Z\) matches naïve (max \(|\Delta|<10^{-12}\)); \(U=n_x n_y=16\).
- Bench: `node --test tests/green-bench.test.js` (SIMD, \(0.5\lambda\), default \(h,\ell,a\)).

| array | \(N\) | Z fill ms | `from_z` LU ms | total ms |
| --- | --- | --- | --- | --- |
| 8×8 | 64 | 0.1 | 2.2 | 2.2 |
| 16×16 | 256 | 0.2 | 97 | 98 |
| 32×32 | 1024 | 11 | 21000 | 21000 |

- Next: WP4 \(F^\mathrm{iso}\) pattern. Do not change GUI. Do not raise `MATCHED_AUTO_ISOLATE_N` here.

```text
WP3 done. Next: WP4.
Z fill ms @ 8/16/32: 0.1 / 0.2 / 11
from_z LU ms @ 8/16/32: 2.2 / 97 / 21000
32×32 conclusion: fill is not the LU (~11 ms vs ~21 s).
```

---

## 10. WP4 — \(F^\mathrm{iso}\) as element pattern

**Owner intent:** AF intensity × PEC-dipole \(|F^\mathrm{iso}(\hat r)|^2\) instead of
`PATTERN_COS_N` when Green PEC mode is active. Still one template × geometric phase.

### Read

- `wasm/src/element.rs` `apply_element_pattern`
- `js/phasedarray/farfield.js` `apply_element_pattern`
- `js/phasedarray/phasedarray.js` `create_farfield_vectors`
- WP1 `f_iso_pec_dipole`

### Do

1. Add a pattern kind or a parallel `apply_green_pec_pattern(domain, ax1, ax2, total, h, ell, freq_scale)` that **does not** remove isotropic/cos^n.
2. Same hemisphere / UV-exterior conventions as current cos^n (zero back / invisible) unless WP1 documented otherwise — then match WP1.
3. Tests: boresight vs horizon shape for \(h=\lambda/4\); UV outside unit circle is 0.
4. Wire farfield JS to call it when a kernel flag / element kind says Green PEC.
   If GUI flag does not exist yet, add a **minimal** kind constant and a test hook;
   WP5 owns the dropdown.

### Do not

- Polarized \(E_\theta,E_\phi\) AF rewrite (optional later).
- \(F^\mathrm{emb}\) on the quadrature grid (still \(w=Ta\)).

### Done when

- Isolated Green PEC (T = I) plot uses dipole-over-ground pattern, not \(\cos^n\).
- Matched Green PEC still uses \(w=Ta\) then the **same** \(F^\mathrm{iso}\).

### Handoff

```text
WP4 done. Next: WP5.
Pattern apply: apply_green_pec_pattern (wasm); FarfieldABC.apply_element_pattern
  dispatches on PATTERN_GREEN_PEC.
GUI still needs Matching=Green PEC + h,ℓ,a controls.
```

### WP4 results (implemented)

- Files: `wasm/src/metrics.rs` (`direction_theta_phi`); `wasm/src/element.rs` (`apply_green_pec_pattern`); `js/phasedarray/element.js` (`PATTERN_GREEN_PEC`, `ElementGreenPec` test hook, not in `ElementTypes`); `js/phasedarray/farfield.js` dispatch; `js/wasm/init.js` wrapper. Rebuild: `./wasm/build.ps1`.
- Signature:

```text
apply_green_pec_pattern(domain, ax1, ax2, total, h, ell, freq_scale) -> peak
  // |F^iso|² of WP1 f_iso_pec_dipole_power after look()→(θ,φ)
  // UV exterior / unknown domain: factor 0. Does not replace isotropic/cos^n.
```

- Isolated \(T=I\) and matched \(w=Ta\) both use the same post-AF multiply (`create_farfield_vectors` unchanged). No GUI. No Gram.
- Tests: `cargo test --manifest-path wasm/Cargo.toml` (4 new `element::tests::green_pec_*`); `node --test tests/farfield-directivity.test.js`.
- Next: WP5 Matching=Green PEC + \(h,\ell,a\) controls. Do not change steer illumination-after-\(T\).

```text
WP4 done. Next: WP5.
Pattern apply: apply_green_pec_pattern (wasm); FarfieldABC.apply_element_pattern
  dispatches on PATTERN_GREEN_PEC.
GUI still needs Matching=Green PEC + h,ℓ,a controls.
```

---

## 11. WP5 — GUI, safeguard, docs

**Owner intent:** a Matching option that uses WP2–4 end-to-end.

### Read

- `js/index-scenes.js` coupling dropdown, `applyCouplingSizeSafeguard`
- `js/phasedarray/matched.js` `MATCHED_AUTO_ISOLATE_N`
- `CHANGE.md`, `docs/radiated_power_gram.md` Matching paragraph
- WP3 timing table

### Do

1. Matching value e.g. `green-pec` (keep isolated / per-port / common / propagation).
2. Controls: \(h\), \(\ell\), \(a\); hide \(X_{nn},\alpha,\beta,\varepsilon_x\) in this mode.
3. Common real \(z_c\) as today (do not revive complex \(z_0\)).
4. Size policy: either raise auto-isolate for Green PEC to 1024 with a warning, or
   leave 512 and document that 32×32 is a **manual** override. Prefer: auto-isolate
   at 1024, warn in UI that first \(S,T\) build may take ~seconds.
5. `CHANGE.md` + a short “Green PEC (implemented)” section at the top of this file.
6. Tests: `tests/matched-s.test.js` analogue for green-pec 8×8 \(S_{ii}\).

### Do not

- Slab/patch.
- Changing steer illumination-after-\(T\) (known P2; out of scope unless trivial).

### Handoff

```text
WP5 done. First user-visible Green mode is PEC dipole (Element dropdown).
Matching stays Isolated / Per-port / Common Z0 / Propagation for isotropic/cos^n.
Next optional: WP6 spectral validation, then WP7 slab.
```

### WP5 results (implemented)

- Files: `js/phasedarray/element.js` (`ElementGreenPec` in `ElementTypes`, title PEC dipole, controls \(h,\ell,a,z_c\), series Self X); `js/index-scenes.js` hides Matching when PEC dipole is selected, `form_green_pec_dipole` for \(S,T\), Match button from `z_self_pec_dipole`; `js/phasedarray/matched.js` `GREEN_PEC_AUTO_ISOLATE_N = 1024` (does not change `MATCHED_AUTO_ISOLATE_N = 512`).
- GUI: Matching select is **not** rewritten. Switching back to isotropic/\(\cos^n\) restores the previous toy. Note under Element warns that 32×32 LU can take seconds. Self X **adds** to physical \(\Im Z_{ii}\) (unlike the Gram toys, which write the diagonal). Match is single-element, not scan-matched.
- Tests: `node --test tests/matched-s.test.js` (Green N=1 leftover \(|S_{11}|\); N=1 cancelled \(X_{11}\) opens \(S\); 8×8 reciprocal finite \(S_{ii}\), Gram unused).

```text
WP5 done. First user-visible Green mode is PEC dipole (Element dropdown).
Matching stays Isolated / Per-port / Common Z0 / Propagation for isotropic/cos^n.
Next optional: WP6 spectral validation, then WP7 slab.
```

---

## 12. WP6 — Spectral PEC validation (optional gate)

**Owner intent:** prove the Sommerfeld path against WP1 closed form **before** slab.

### Do

- Visible + evanescent \(k_\rho\) integral of free-space/PEC spectral \(G\) for the
  same dipole; compare \(Z_pair\) and \(F^\mathrm{iso}\) to WP1.
- Same unique-lag driver as WP3.
- Stop. Do not add \(\varepsilon_r\).

If this is slow, fix quadrature here (not in WP7).

### WP6 results (implemented)

- Files: `wasm/src/green_spectral.rs` (Sommerfeld \(Z\) + saddle-point \(F^\mathrm{iso}\) + tests); `mod green_spectral;` in `wasm/src/lib.rs`; `wasm/src/prad.rs` unique-lag fill takes a \(z_\mathrm{pair}\) callback so the WP3 driver can call the spectral kernel in tests. `green.rs` closed form is still the production path. No wasm-bindgen, no GUI, no \(\varepsilon_r\).
- Signatures:

```text
SpectralQuadConfig { n_k_prop, n_k_evan, k_evan_max_over_k, n_lobes }
DEFAULT: 48, 32, 12, 40
z_pair_pec_dipole_spectral(dx, dy, h, ell, a, freq_scale) -> (re, im)
f_iso_pec_dipole_spectral(theta, phi, h, ell, freq_scale) -> (Eθ, Eφ)
```

- Quadrature: propagating \(\alpha\) (\(k_\rho=k\sin\alpha\)) and evanescent \(\beta\) (\(k_\rho=k\cosh\beta\)) remove the \(k_z=0\) branch point. Angular integral is \(J_0,J_2\) (no \(\phi_k\) grid). Coplanar \(z=0\) tail uses Bessel lobes + Wynn \(\varepsilon\). Mutual is Hertzian (moment \(I\ell\)), same as WP1; self is \(Z_\mathrm{fs}-Z(0,0,2h)\). \(F^\mathrm{iso}\) is the space-wave saddle of the same spectral \(\tilde E\).
- Errors vs `green.rs` (release): self \(|\Delta Z|\sim 3\times 10^{-14}\); mutual \(\sim 10^{-8}\)–\(10^{-10}\,\Omega\) on the WP6 lag grid. \(|F|^2\) matches to \(10^{-12}\). 4×4 unique-lag \(Z\) and `form_from_z` \(S,T\) (faer) agree to \(<10^{-3}\).
- Per-lag (native release): closed \(\sim 0.001\,\mu\mathrm{s}\); spectral \(\sim 50\,\mu\mathrm{s}\) (self) / \(\sim 200\)–\(300\,\mu\mathrm{s}\) (mutual). Production stays closed-form.

```text
WP6 done. Next: WP7 grounded slab.
Spectral gate: green_spectral.rs matches green.rs within 1e-8 Ω typical (self 3e-14).
Defaults: n_k_prop=48, n_k_evan=32, k_evan_max/k=12, n_lobes=40
Per-lag: closed ~0.001 µs, spectral ~50–300 µs.
Production: still z_pair_pec_dipole (closed form).
WP7: reuse green_spectral.rs integrator; add slab Γ(kρ) + TM0 residue. Do not add ε_r here.
```

---

## 13. WP7 — Grounded dielectric slab \(G\)

Same `Z_pair` / `F_iso` interface. Sommerfeld + TM0 residue. New env params
\(\varepsilon_r,h_\mathrm{sub},\tan\delta\). Radiator still the WP1 dipole **or** a
magnetic Hertzian slot (choose one and document). Gram still not \(\Re(Z)\).
Surface-wave power is not AF intensity.

**Owner intent:** slab spectral \(G\) for the WP1 \(+\hat x\) short dipole (not
the magnetic slot). Power budget split from the same \(G\). No GUI (WP7b).

### WP7 results (implemented)

- Radiator: WP1 finite short dipole. Environment: air \(z>0\), dipole at \(z=h\),
  dielectric \(-h_\mathrm{sub}<z<0\) with \(\varepsilon=\varepsilon_r(1-j\tan\delta)\),
  PEC at \(z=-h_\mathrm{sub}\).
- Files: `wasm/src/green_slab.rs`; `mod green_slab;` in `wasm/src/lib.rs`;
  `wasm/src/prad.rs` test helpers `fill_green_slab_dipole_z` /
  `form_green_slab_dipole` (unique-lag callback, no bindgen). Production PEC
  path is still `green.rs`. No GUI, no wasm-bindgen.
- Signatures:

```text
SlabEnv { eps_r, h_sub, tan_delta }
DEFAULT: 10, 0.05, 0     PEC_LIMIT: 12, 1e-6, 0
z_pair_slab_dipole(dx, dy, h, ell, a, freq_scale, env) -> (re, im)
f_iso_slab_dipole(theta, phi, h, ell, freq_scale, env) -> (Eθ, Eφ)
f_iso_slab_dipole_power(...) -> |F|²
SlabPowerBudget { re_z_self, p_rad, p_sw, p_diss, closure_residual }
slab_dipole_power_budget(h, ell, a, freq_scale, env)
```

- \(Z\): free-space Hertzian (self \(Z_\mathrm{fs}\), mutual coplanar Sommerfeld)
  plus slab reflection \(\Gamma_\mathrm{TE}(k_\rho),\Gamma_\mathrm{TM}(k_\rho)\)
  with analytic TM0 residue. \(h_\mathrm{sub}\to 0\) recovers WP1 PEC
  (\(|\Delta Z|<10^{-3}\,\Omega\)).
- \(F^\mathrm{iso}\): space-wave saddle of the same spectral \(\tilde E\)
  (\(e^{j k_z h}+\Gamma e^{-j k_z h}\), TE/TM split). Back hemisphere zero.
  No surface-wave term in the pattern.
- Power (unit current; ohm \(\equiv\) watt in this basis):
  \(P_\mathrm{rad}\) from the propagating \(k_\rho<k\) branch,
  \(P_\mathrm{sw}\) from the TM0 indentation jump \(-j\pi\operatorname{Res}\),
  \(P_\mathrm{diss}=\operatorname{Re}Z-P_\mathrm{rad}-P_\mathrm{sw}\)
  (evan continuum / \(\tan\delta\)). Lossless: \(P_\mathrm{diss}\approx 0\),
  \(P_\mathrm{sw}>0\) for the default slab. **EIRP stays**
  \(D\times P_\mathrm{stimulated}\) (WP7b must not scale EIRP by \(P_\mathrm{rad}\)).
- Per-lag (native release, default env): self \(\sim 290\,\mu\mathrm{s}\);
  mutual \(\sim 0.9\,\mathrm{ms}\). 32×32 unique-lag fill \(\sim 0.3\)–\(1\,\mathrm{s}\)
  native; WASM will be slower. Geometry-change only, not per-steer.

```text
WP7 done. Next: WP7b GUI slab element.
API: z_pair_slab_dipole / f_iso_slab_dipole / slab_dipole_power_budget
File: wasm/src/green_slab.rs
Radiator: WP1 +x short dipole (not magnetic slot). Slot is WP8.
EIRP: P_stimulated, not P_rad.
No bindgen, no GUI.
WP7b: Element "Slab dipole", controls (h,ℓ,a,ε_r,h_sub,tanδ,z0,xself),
  form_green_slab_dipole bindgen, apply_green_slab_pattern,
  power rows Radiated / Surface wave / Substrate loss / Closure residual.
  Do not change eirp = dirMax * pAcc.
```

### WP7b results (implemented)

- Files: `js/phasedarray/element.js` (`ElementGreenSlab`, `PATTERN_GREEN_SLAB`); `js/index-scenes.js` (`isGreenElement`, `form_green_slab_dipole`); `js/phasedarray/farfield.js` + `js/wasm/init.js`; `wasm/src/element.rs` (`apply_green_slab_pattern`, `z_self_slab_dipole`, `slab_dipole_power_budget_wasm`); `wasm/src/prad.rs` / `lib.rs` bindgen of `form_green_slab_dipole`. Power rows in `index.html` / `js/index.js`. Rebuild: `./wasm/build.ps1`.
- GUI: Element **Slab dipole** with \(h,\ell,a,\varepsilon_r,h_\mathrm{sub},\tan\delta,z_c,X_\mathrm{self}\). Matching hidden. Match sets \(z_c=\Re Z_{11}\) and \(X_\mathrm{self}=-\Im Z_{11}\) of the isolated slab self kernel. Auto-isolate at \(N>1024\) (same as PEC). Surface-wave power is not in the AF pattern.
- Power: Radiated / Surface wave / Substrate loss / Closure residual from `slab_dipole_power_budget` at \(|I|=1\,\mathrm{A}\). **EIRP stays** `dirMax * pAcc`.
- Tests: `cargo test --manifest-path wasm/Cargo.toml` (`element::tests::green_slab_*`); `node --test tests/matched-s.test.js`; `node --test tests/farfield-directivity.test.js`.

```text
WP7b done. Next: WP8 two-slot patch.
API: form_green_slab_dipole / apply_green_slab_pattern / z_self_slab_dipole
     slab_dipole_power_budget_wasm -> [re_z_self, p_rad, p_sw, p_diss, closure]
Element: Slab dipole (PATTERN_GREEN_SLAB)
EIRP: unchanged dirMax * pAcc
Radiator: WP1 +x short dipole. Slot is WP8.
```

---

## 14. WP8 — Patch as two slots on WP7 \(G\)

Cavity two-slot equivalent \(\mathbf{M}\), one port per element, same dense \(Z\)
and \(F^\mathrm{iso}\). No metal MoM. Same \(S,T\), AF path.

---

## 15. Context pack (every agent)

Minimum files to open:

1. `docs/green_function_plan.md` (this file) — assigned WP only.
2. Previous WP handoff (last section of that session / git diff of listed files).
3. WP-specific list above.

Do not re-read the full audit unless the WP mentions physics policy.
