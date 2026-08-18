# Change Log
## 8/18/2026
- Added an Element control: isotropic (default, +3 dBi hemisphere) or a power-conserving cos^n pattern from a user peak gain in dBi. The exponent n = 10^(G/10)/2 − 1 is shown read-only. The pattern multiplies far-field intensity after the array-factor sum, so calculated directivity, EIRP, and beam metrics include the element. Cos^n gain floors at 3.01 dBi (10·log10 2) and zeros the back hemisphere / U-V exterior.

## 8/17/2026
- Added a Power control block: value + unit dropdown (default 20 dBm), per-element vs array radio, and a separate EIRP unit (default dBW).
- Element hover shows per-element power. 2-D pattern header, hover, and metrics show EIRP from directivity times taper-weighted array power.

## 8/8/2026
- Added U/V bounds so that plot can show outside +/- 1. Added visible circle to U/V plots.
- Added illumination feeds and illumination plots.
