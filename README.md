# PhasedArrayVisualizer

[`PhasedArrayVisualizer`](https://jasondurbin.github.io/PhasedArrayVisualizer/) is a javascript application for simulating and visualizing phased array radiation patterns. This visualizer was created by [Jason Durbin](https://www.linkedin.com/in/jasondurbin/).

[Try it!](https://jasondurbin.github.io/PhasedArrayVisualizer/)

While this is intended to be a demo, you're welcome to view the source code.

## Rust / WebAssembly far-field kernel

The far-field array-factor sum (the nested element × grid loop in spherical, U-V, and Ludwig3 domains) is implemented in Rust and compiled to WebAssembly. Geometry, illumination, steering, tapers, quantization, and plotting remain in JavaScript. JS still builds per-element positions plus polar weights `(magnitude, phase)` and the observation grid; WASM accumulates the complex field and returns intensity.

The crate lives in [`wasm/`](wasm/). Two artifacts are committed under [`js/wasm/`](js/wasm/): a `simd128` build and a scalar fallback. The page picks SIMD when the browser supports it.

To rebuild after changing the Rust kernel (Rust stable, `wasm32-unknown-unknown`, and [`wasm-pack`](https://rustwasm.github.io/wasm-pack/) required):

```powershell
rustup target add wasm32-unknown-unknown
cargo install wasm-pack
./wasm/build.ps1
```

Equivalence tests compare the original JavaScript loops to both WASM builds:

```powershell
node --test tests/farfield-equivalence.test.js
```

## License

You are free to edit, share, change, iterate for **educational and non-commercial** uses. For commercial uses (such as embedding on your site, sharing with customers to simulate your own array, etc), please contact hello@neonphysics.com.

## Donation

If you enjoy the Visualizer, consider donating [using PayPal](https://www.paypal.com/donate/?business=D7S3JKRAAKUNQ&no_recurring=0&currency_code=USD).
