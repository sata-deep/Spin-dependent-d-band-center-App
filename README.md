# Spin Dependent d-band Center 

A simple Python/Tkinter GUI to compute spin-resolved and effective d-band centers directly from a VASP `DOSCAR` (ISPIN=2, LORBIT=11). 

## Features
- Load a `DOSCAR` (ISPIN=2, LORBIT=11).
- Select atom range (Group of atoms using comma, such as 2-5,7-8,.. etc)
- Compute:
  - Spin-up d-band center ($\varepsilon_{d\uparrow}$)
  - Spin-down d-band center ($\varepsilon_{d\downarrow}$)
  - Effective d-band center ($\varepsilon_{\mathrm{eff}}$)
  - Hammer–Nørskov d-band center, $\varepsilon_{\mathrm{HN}}$
  - Fractional occupations $f_\uparrow$, $f_\downarrow$

## Equations
Let energies be aligned such that \(E_F = 0\), i.e., we use \(E - E_F\) internally.

- Spin-resolved d-band centers (integrated over the full energy grid):
  
$$
  \varepsilon_{d\sigma} 
  = \frac{\int_{-\infty}^{+\infty} E D_{d\sigma}(E - E_F) dE}
         {\int_{-\infty}^{+\infty} D_{d\sigma}(E - E_F) dE}
$$

- Fractional occupations (integrated up to the Fermi level):

 
 $$
  f_\sigma = \frac{1}{N_{\text{atoms}}\times 5} \int_{-\infty}^{E_F} D_{d\sigma}(E) dE
  \quad (E_F = 0)
  $$

- Effective d-band center (**Bhattacharjee, S.**, Waghmare, U. & Lee, SC. _An improved d-band model of the catalytic activity of magnetic transition metal surfaces_. Sci Rep 6, 35916 (2016). https://doi.org/10.1038/srep35916):

$$
  \varepsilon_{\mathrm{eff}} =
  \frac{f_\uparrow\varepsilon_{d\uparrow} + f_\downarrow \varepsilon_{d\downarrow}}
       {f_\uparrow + f_\downarrow}-(\varepsilon_{d\downarrow} - \varepsilon_{d\uparrow})
    \frac{f_\uparrow - f_\downarrow}{f_\uparrow + f_\downarrow}
  $$

- Hammer–Nørskov d-band center (spin-summed, full energy range):

$$
  \varepsilon_{\mathrm{HN}} =
  \frac{\int E \big(D_{d\uparrow}(E)+D_{d\downarrow}(E)\big) dE}
       {\int \big(D_{d\uparrow}(E)+D_{d\downarrow}(E)\big) dE}
  $$



Optionally, the reduced fractional occupation is

$$
\mu = \frac{f_\uparrow - f_\downarrow}{f_\uparrow + f_\downarrow}.
$$

## Installation
- Python 3.8+ (standard library only).
- Tkinter must be available (usually bundled with Python). On Debian/Ubuntu:
  - `sudo apt-get install python3-tk`

## For ubuntu
- Just type in your terminal:
  
  **sudo apt update**
  
  **pip install --upgrade Spin-dependent-d-band-center**
- Then in bash type:
  
  **d-band-satadeep**
  
 --> This  will launch the application (see the pic below)

## Usage: If you don't want to use pip install
```
python3 d-band.py
```
- Click “Browse…” and select your `DOSCAR` (generated with `ISPIN=2`, `LORBIT=11`).
- Enter the atom ranges (e.g., 3-5,8-7,22-25 etc).
- Click “Compute”.
- Once executed the window will look like the following:
<img width="907" height="795" alt="image" src="https://github.com/user-attachments/assets/5ad4dd97-c52c-455b-ab5c-ef17da6a979d" />



## Notes
- The app shifts energies so that \(E_F=0\). Integrals for \(f_\sigma\) run from the grid minimum to 0; centers use the full grid range.
- The code aggregates the five d-orbitals (dxy, dyz, dz2, dxz, dx2) for each spin channel.
- Results shown in eV.

## Credits
Created by Satadeep Bhattacharjee, IKST, Bangalore

