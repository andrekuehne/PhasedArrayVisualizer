# Change Log
## 8/19/2026
- Element has a Slab dipole option: a finite short dipole over a grounded dielectric slab. Mutual Z and the isolated pattern come from the same Sommerfeld Green function; Matching is hidden as with PEC dipole. Extra controls are ε_r, h_sub, and tan δ. Power shows isolated-element Radiated / Surface wave / Substrate loss / Closure residual at 1 A; Peak EIRP stays directivity × accepted power.
- Matched \(S,T\) now uses faer (blocked complex LU and Hermitian Cholesky, single-threaded). The GUI matching path (`compute_matched_basis`: Per-port, Common Z0, Propagation, and PEC dipole) always calls the WASM kernel, which always uses faer. Textbook LU/Cholesky stays in `wasm/src/legacy_linalg.rs` for parity tests only. On a 32×32 Green PEC array (\(N=1024\)), WASM `from_z` dropped from 16.8 s to 1.8 s (~9.5×); native LU is 1.80 s → 0.26 s (~7×). SIMD wasm grows from 152 KB to 488 KB.
- Element has a PEC dipole option: a finite short dipole over infinite ground. Mutual Z and the isolated pattern come from the same Green function; Matching (Isolated / Per-port / Common Z0 / Propagation) is hidden in this mode. Controls are height h, length ℓ, radius a (wavelengths at f₀), a real port Z0 (default 50 Ω), and Self X (series reactance added to diag(Z), default 0). Match sets Z0 = Re(Z11) and Self X = −Im(Z11) of the isolated element. First S,T build is about 2 s at 32×32; arrays larger than 1024 elements keep T = I and still apply the dipole pattern.
- Matching is a single dropdown: Isolated, Per-port, Common Z0, and Propagation. Legacy URLs with Coupling=Matched map Per-port / Common Z0 from the old match-style param.
- Propagation is a new matched model: mutual X from effective \(\varepsilon_x,\varepsilon_y\) phase and wavelength decay \(\alpha_\lambda\), always at a common real Z0. Per-port and Common Z0 keep the existing Mutual X / \(\alpha\) / \(\beta\) / A kernel.
- Matched coupling now keeps a purely real port reference. Self X (Ω) is a common diagonal residual on Z (the old Im(Z0) control); Mutual X is the pairwise overlay. Per-port matches Re(Z_in) only, so leftover reactance stays in Z and in S_ii. Common Z0 is a single real z_c (default 45 Ω).
- Matched coupling's optional mutual reactance is now a pairwise kernel on (Δx, Δy): the old X_nn (d_min/ρ)^α tail plus Oscillation β (sign-changing with distance) and Location A (x vs y at the same distance). Defaults β = 0, A = 0 keep the previous overlay. Works for sunflower and other irregular layouts; β is scaled by frequency.

## 8/18/2026
- Added an Element control: isotropic (default, +3 dBi hemisphere) or a power-conserving cos^n pattern from a user peak gain in dBi. The exponent n = 10^(G/10)/2 − 1 is shown read-only. The pattern multiplies far-field intensity after the array-factor sum, so calculated directivity, EIRP, and beam metrics include the element. Cos^n gain floors at 3.01 dBi (10·log10 2) and zeros the back hemisphere / U-V exterior.
- Frequency scale now also scales geometric k in U/V and Ludwig-3 (same 2π·f/f₀ as spherical). Commanded element phase is still unscaled, so phase-shifter beam squint appears in all domains.

## 8/17/2026
- Added a Power control block: value + unit dropdown (default 20 dBm), per-element vs array radio, and a separate EIRP unit (default dBW).
- Element hover shows per-element power. 2-D pattern header, hover, and metrics show EIRP from directivity times taper-weighted array power.

## 8/8/2026
- Added U/V bounds so that plot can show outside +/- 1. Added visible circle to U/V plots.
- Added illumination feeds and illumination plots.
