# Change Log
## 8/19/2026
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
