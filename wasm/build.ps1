$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

Write-Host "Building SIMD wasm..."
$env:RUSTFLAGS = "-C target-feature=+simd128"
wasm-pack build --target web --release --no-pack --out-dir ../js/wasm/simd

Write-Host "Building scalar wasm..."
Remove-Item Env:RUSTFLAGS -ErrorAction SilentlyContinue
cargo clean -p farfield_kernel
wasm-pack build --target web --release --no-pack --out-dir ../js/wasm/scalar

Remove-Item ../js/wasm/simd/.gitignore -ErrorAction SilentlyContinue
Remove-Item ../js/wasm/scalar/.gitignore -ErrorAction SilentlyContinue
Write-Host "Done. Outputs in js/wasm/simd and js/wasm/scalar."
