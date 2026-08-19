# Green-function coupling: sequential work packages

Handoff plan for finite-array Green coupling. **No Floquet / infinite-array path.**
\(F^\mathrm{iso}\) and complex \(Z\) come from the same environment Green function;
existing \(S\), \(T\), and array-factor machinery stay. First physics: **horizontal
short dipole over infinite PEC** (finite \(\ell\), not a point Hertzian). Later:
patch over grounded slab, same interface.

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
Pattern apply: <function names>.
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
WP5 done. First user-visible Green mode is PEC dipole.
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

---

## 13. WP7 — Grounded dielectric slab \(G\)

Same `Z_pair` / `F_iso` interface. Sommerfeld + TM0 residue. New env params
\(\varepsilon_r,h_\mathrm{sub},\tan\delta\). Radiator still the WP1 dipole **or** a
magnetic Hertzian slot (choose one and document). Gram still not \(\Re(Z)\).
Surface-wave power is not AF intensity.

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
