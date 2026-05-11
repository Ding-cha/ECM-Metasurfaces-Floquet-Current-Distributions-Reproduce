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

Create the geometry and current plots:

```powershell
python -m symmetric_cross_mts.plot_current
```

Run the ECM frequency sweep:

```powershell
python -m symmetric_cross_mts.run_ecm
```

By default this writes the current restored ECM candidate. This candidate uses:

```text
TTR variant: TE/TM unit-polarization projection
evanescent branch: decaying branch
high-order substrate loading: disabled for evanescent modal sum
mode order: 9
```

This is the closest fully analytical reconstruction found so far. It naturally
reproduces the main Fig. 15 features: the reflection peak / transmission null
near 5.8 GHz and the reflection minimum near 8 GHz.

To overlay the digitized anchors:

```powershell
python -m symmetric_cross_mts.run_ecm --anchors
```

To run the hand-digitized Fig. 15 calibration:

```powershell
python -m symmetric_cross_mts.run_ecm --mode calibrated --anchors
```

To run the raw analytical ECM implementation:

```powershell
python -m symmetric_cross_mts.run_ecm --mode raw
```

Validate the Fig. 22 width sweep (`w = 3, 5, 7 mm`):

```powershell
python -m symmetric_cross_mts.validate_width_sweep
```

This writes:

```text
outputs/width_sweep_s11.png
outputs/width_sweep_summary.csv
```

Diagnose why the current analytical implementation does not yet match the
published Fig. 15 curve:

```powershell
python -m symmetric_cross_mts.diagnose_ecm
```

Outputs are written to:

```text
outputs/geometry.png
outputs/current_distribution.png
outputs/s_parameters.png
outputs/s_parameters.csv
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
