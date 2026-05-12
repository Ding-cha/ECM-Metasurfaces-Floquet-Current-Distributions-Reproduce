# Symmetric Cross MTS Unit Cell

This project reproduces the symmetric-cross metasurface unit cell from:

Baladi and Hum, "Equivalent Circuit Models for Metasurfaces Using Floquet Modal
Expansion of Surface Current Distributions," IEEE TAP, 2021.

The default parameters match Fig. 13 and Table II in the paper.

## Geometry

| Parameter | Value |
| --- | ---: |
| `l_x` | 20 mm |
| `l_y` | 20 mm |
| `w` | 3 mm |
| `P_x` | 30 mm |
| `P_y` | 30 mm |
| `h` | 1.524 mm |
| `epsilon_r` | 3.02 |

## Current Model

For x-polarized normal incidence, the paper uses:

```text
I_1x(x) = 0.95 cos(0.12 x) - 0.2 cos(0.64 x)
I_2y(y) = 0.5 I_1x(0) |(|y| - L_y/2) / (L_y/2)|
```

where `x` and `y` are in mm.

## Setup in VS Code

Open this folder in VS Code. The workspace is configured to use:

```text
D:\anaconda3\python.exe
```

This interpreter already includes NumPy and Matplotlib on this machine. If you
prefer a virtual environment, run:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m pip install -e .
```

The `.vscode` folder includes launch and task configurations for the two main
scripts. They set `PYTHONPATH=src`, so editable installation is optional.

## Run

From a normal PowerShell terminal in this repository, prefix commands with:

```powershell
$env:PYTHONPATH='src';
```

For example:

Create the geometry and current plots:

```powershell
$env:PYTHONPATH='src'; python -m symmetric_cross_mts.plot_current
```

Run the ECM frequency sweep:

```powershell
$env:PYTHONPATH='src'; python -m symmetric_cross_mts.run_ecm
```

By default this writes the current restored ECM candidate. This candidate uses:

```text
TTR variant: TE/TM unit-polarization projection
evanescent branch: decaying branch
restored profile: auto
```

The restored `auto` profile selects between two analytical assumptions:

```text
thin substrate:      mode_order=9, high-order substrate loading off, cosine width profile
thicker substrate:   mode_order=3, high-order substrate loading on, uniform width profile,
                     thickness-dependent epsilon_eff_scale,
                     arm_length_correction_mm=0.5
```

The threshold is `h_mm * sqrt(epsilon_r) >= 4.0`. This keeps the Fig. 15-style
thin-substrate behavior while improving agreement for thicker / higher-permittivity
substrates such as `h=3 mm, epsilon_r=4.4`.

The thick-substrate profile uses:

```text
epsilon_eff = 1 + epsilon_eff_scale * (epsilon_r - 1)
epsilon_eff_scale = clamp(1.03 - 0.04 * h_mm, 0.88, 0.97)
```

This accounts for surface-current fields that are only partly confined inside
the dielectric. The arm-length correction accounts for the fact that the fitted
surface-current model overestimates the electrical length in the lx=ly sweep.
Use `--epsilon-eff-scale 1 --arm-length-correction-mm 0` to recover the earlier
bulk-permittivity, full-arm current model.

To overlay the digitized anchors:

```powershell
$env:PYTHONPATH='src'; python -m symmetric_cross_mts.run_ecm --anchors
```

Run an arbitrary symmetric-cross geometry by passing dimensions in mm:

```powershell
$env:PYTHONPATH='src'; python -m symmetric_cross_mts.run_ecm --mode restored --lx-mm 24 --ly-mm 18 --w-mm 4 --px-mm 34 --py-mm 30 --h-mm 1.0 --epsilon-r 2.5 --freq-start-ghz 4 --freq-stop-ghz 10
```

Force the thick-substrate restored settings when comparing against HFSS cases
where dielectric loading of high-order modes is important:

```powershell
$env:PYTHONPATH='src'; python -m symmetric_cross_mts.run_ecm --mode restored --restored-profile substrate --params-json examples\symmetric_cross_case.json
```

Or read the geometry from JSON:

```powershell
$env:PYTHONPATH='src'; python -m symmetric_cross_mts.run_ecm --mode restored --params-json examples\symmetric_cross_case.json
```

The JSON file may use either the code names (`lx_mm`, `py_mm`, `epsilon_r`) or
short paper-style aliases (`l_x`, `p_y`, `eps_r`). CLI geometry values override
the JSON file. Each run writes both the response and the resolved configuration:

```text
outputs/s_parameters.png
outputs/s_parameters.csv
outputs/run_config.json
```

To run the hand-digitized Fig. 15 calibration:

```powershell
$env:PYTHONPATH='src'; python -m symmetric_cross_mts.run_ecm --mode calibrated --anchors
```

To run the raw analytical ECM implementation:

```powershell
$env:PYTHONPATH='src'; python -m symmetric_cross_mts.run_ecm --mode raw
```

Validate the Fig. 22 width sweep (`w = 3, 5, 7 mm`):

```powershell
$env:PYTHONPATH='src'; python -m symmetric_cross_mts.validate_width_sweep
```

This writes:

```text
outputs/width_sweep_s11.png
outputs/width_sweep_summary.csv
```

Run the small batch sweep for `lx = ly = 20..25 mm`, `w = 3, 4, 5 mm`, and
`h = 1, 2, 3 mm`:

```powershell
$env:PYTHONPATH='src'; python -m symmetric_cross_mts.batch_sweep --output-dir outputs\batch_lx_w_h_test
```

For a faster smoke test, reduce frequency samples:

```powershell
$env:PYTHONPATH='src'; python -m symmetric_cross_mts.batch_sweep --freq-points 61 --fourier-samples 201 --output-dir outputs\batch_lx_w_h_test
```

The batch run writes:

```text
outputs/batch_lx_w_h_test/batch_config.json
outputs/batch_lx_w_h_test/summary.csv
outputs/batch_lx_w_h_test/cases/*.csv
```

Diagnose why the current analytical implementation does not yet match the
published Fig. 15 curve:

```powershell
$env:PYTHONPATH='src'; python -m symmetric_cross_mts.diagnose_ecm
```

Outputs are written to:

```text
outputs/geometry.png
outputs/current_distribution.png
outputs/s_parameters.png
outputs/s_parameters.csv
outputs/run_config.json
```

## Notes

The ECM implementation follows the paper's Floquet modal projection workflow and
uses a normalized current amplitude. It is meant as a reproducible analytical
starting point for the symmetric-cross unit cell. Small differences from the
published HFSS overlay can be expected because the paper does not publish all
implementation details, such as exact metal thickness, loss tangent, modal
truncation, numerical normalization choices, the incident-current normalization
`I0`, and the exact normal-incidence `TE00` limiting procedure. Without those
normalization details, the raw Fourier-integral `T00` is too small and the model
mostly reproduces the bare dielectric slab response.
