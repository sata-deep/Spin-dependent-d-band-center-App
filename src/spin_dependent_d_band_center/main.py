import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from dataclasses import dataclass
from typing import Optional, Tuple

import math


@dataclass
class DOSData:
    # Energies shifted so that E_F = 0
    energies: list
    # Per-atom aggregated d-DOS (sum over five d orbitals) for spin up and spin down
    d_up: list  # shape: [natom][NE]
    d_dn: list  # shape: [natom][NE]

    @property
    def natom(self) -> int:
        return len(self.d_up)

    @property
    def ne(self) -> int:
        return len(self.energies)


class DOSCARParser:
    """Parse VASP DOSCAR (LORBIT=11, ISPIN=2) and build minimal d-DOS arrays per atom.

          d order: dxy, dyz, dz2, dxz, dx2
    """

    def parse(self, path: str) -> DOSData:
        with open(path, "r") as f:
            lines = f.read().strip().splitlines()

        if len(lines) < 7:
            raise ValueError("DOSCAR looks too short.")

        # First line contains natom twice in Fortran reader; be tolerant
        first_tokens = lines[0].split()
        if len(first_tokens) < 1:
            raise ValueError("Unable to parse natom from DOSCAR first line.")
        try:
            natom = int(float(first_tokens[0]))
        except Exception as e:
            raise ValueError("First token of DOSCAR is not numeric; unsupported DOSCAR format.") from e

        # Skip next 4 header lines
        idx = 1 + 4

        # Read global grid line
        try:
            emax, emin, NE_f, efermi, _kn = map(float, lines[idx].split()[:5])
        except Exception as e:
            raise ValueError("Failed parsing energy grid line.") from e
        NE = int(NE_f)
        idx += 1

        # Read NE lines of total DOS (spin polarized)
        total_block = lines[idx: idx + NE]
        if len(total_block) != NE:
            raise ValueError("DOSCAR truncated in total DOS block.")

        energies = []
        for ln in total_block:
            toks = ln.split()
            if len(toks) < 5:
                raise ValueError("Expected 5 columns in total DOS (E, DOSu, DOSd, IDOSu, IDOSd).")
            energies.append(float(toks[0]))
        idx += NE

        # Prepare per-atom arrays: we only need aggregated d-up and d-down
        d_up = [[0.0 for _ in range(NE)] for _ in range(natom)]
        d_dn = [[0.0 for _ in range(NE)] for _ in range(natom)]

        # For each atom: 1 line header + NE lines PDOS
        for atom in range(natom):
            if idx >= len(lines):
                raise ValueError(f"DOSCAR truncated before atom {atom+1} block.")
            # atom header
            _ = lines[idx]
            idx += 1
            pdos_block = lines[idx: idx + NE]
            if len(pdos_block) != NE:
                raise ValueError(f"DOSCAR truncated in atom {atom+1} PDOS block.")
            for ie, ln in enumerate(pdos_block):
                toks = ln.split()
                # Expect: E + 18 spin-projected columns = 19 numbers total
                if len(toks) < 19:
                    raise ValueError(
                        f"Atom {atom+1}, energy index {ie+1}: expected >=19 columns, got {len(toks)}. Ensure LORBIT=11, ISPIN=2."
                    )
                # Column map (0-based):
                # 0: E
                # 1: s_up, 2:s_dn, 3:py_up,4:py_dn, 5:pz_up,6:pz_dn, 7:px_up,8:px_dn,
                # 9:dxy_up,10:dxy_dn, 11:dyz_up,12:dyz_dn, 13:dz2_up,14:dz2_dn,
                # 15:dxz_up,16:dxz_dn, 17:dx2_up,18:dx2_dn
                d_up_val = (
                    float(toks[9]) + float(toks[11]) + float(toks[13]) + float(toks[15]) + float(toks[17])
                )
                d_dn_val = (
                    float(toks[10]) + float(toks[12]) + float(toks[14]) + float(toks[16]) + float(toks[18])
                )
                d_up[atom][ie] = d_up_val
                d_dn[atom][ie] = d_dn_val
            idx += NE

        # Shift energies to E - E_F
        # efermi from the header line above
        energies_rel = [e - efermi for e in energies]
        return DOSData(energies=energies_rel, d_up=d_up, d_dn=d_dn)


def trapz_integral(x: list, y: list, x_min: float, x_max: float) -> float:
    """Trapezoidal integration of y(x) over [x_min, x_max], assuming x is monotonic."""
    # Clip to range, creating an auxiliary list with linear interpolation at boundaries if needed
    # Build segments where x[i]..x[i+1] overlaps the interval
    total = 0.0
    for i in range(len(x) - 1):
        x0, x1 = x[i], x[i + 1]
        y0, y1 = y[i], y[i + 1]
        # Skip segments completely outside
        if x1 <= x_min or x0 >= x_max:
            continue
        # Determine clipped segment
        xa = max(x0, x_min)
        xb = min(x1, x_max)
        if xb <= xa:
            continue
        # Linear interpolate y at xa and xb
        # Avoid division by zero if x1==x0 (shouldn't happen normally)
        if x1 != x0:
            t_a = (xa - x0) / (x1 - x0)
            t_b = (xb - x0) / (x1 - x0)
            ya = y0 + t_a * (y1 - y0)
            yb = y0 + t_b * (y1 - y0)
        else:
            ya = y0
            yb = y1
        total += 0.5 * (ya + yb) * (xb - xa)
    return total


def compute_d_band_centers(dos: DOSData, atom_range: Tuple[int, int]) -> Tuple[float, float, float, float, float, float]:
    """Compute spin-up, spin-down, effective, and H–N d-band centers.

    Returns (momu, momd, effective, hn_center, fup, fdn)
    """
    n1, n2 = atom_range
    if n1 < 1 or n2 > dos.natom or n1 > n2:
        raise ValueError(f"Invalid atom range: {n1}..{n2} for natom={dos.natom}")

    # Aggregate DOS over selected atoms
    NE = dos.ne
    # Ensure energies are ascending with corresponding indices
    order = sorted(range(NE), key=lambda i: dos.energies[i])
    e = [dos.energies[i] for i in order]
    dosu = [0.0] * NE
    dosd = [0.0] * NE
    for a in range(n1 - 1, n2):
        up = dos.d_up[a]
        dn = dos.d_dn[a]
        for j, i in enumerate(order):
            dosu[j] += up[i]
            dosd[j] += dn[i]

    # Integrate occupations from E_min to E_F=0 (occupied states)
    e_min = e[0]
    e_max_for_occ = 0.0
    # For d-band centers, integrate over the full grid [E_min, E_max]
    e_max_full = e[-1]

    # Occupations (areas under DOS up to E_F), divide by number of states per spin: nat*5
    nat = (n2 - n1 + 1)
    denom_states = nat * 5.0

    int_dosu = trapz_integral(e, dosu, e_min, e_max_for_occ)
    int_dosd = trapz_integral(e, dosd, e_min, e_max_for_occ)
    fup = abs(int_dosu / denom_states)
    fdn = abs(int_dosd / denom_states)

    # Spin-resolved centers over full range: <E> = int(E*DOS) / int(DOS)
    e_dosu = [e[i] * dosu[i] for i in range(NE)]
    e_dosd = [e[i] * dosd[i] for i in range(NE)]
    num_up = trapz_integral(e, e_dosu, e_min, e_max_full)
    den_up = trapz_integral(e, dosu, e_min, e_max_full)
    num_dn = trapz_integral(e, e_dosd, e_min, e_max_full)
    den_dn = trapz_integral(e, dosd, e_min, e_max_full)
    momu = num_up / den_up if den_up and not math.isnan(den_up) else float('nan')
    momd = num_dn / den_dn if den_dn and not math.isnan(den_dn) else float('nan')

    # Effective center per README/doscalculator_sp.f90
    # momav = (fup*momu + fdn*momd)/(fup+fdn) - (momd - momu) * ((fup - fdn)/(fup + fdn))
    denom_f = (fup + fdn)
    if denom_f == 0:
        effective = float('nan')
    else:
        effective = ((fup * momu + fdn * momd) / denom_f) - ((momd - momu) * ((fup - fdn) / denom_f))

    # Hammer–Nørskov: spin-summed center over full range
    dos_sum = [dosu[i] + dosd[i] for i in range(NE)]
    e_dos_sum = [e[i] * dos_sum[i] for i in range(NE)]
    hn_num = trapz_integral(e, e_dos_sum, e_min, e_max_full)
    hn_den = trapz_integral(e, dos_sum, e_min, e_max_full)
    hn_center = hn_num / hn_den if hn_den != 0 else float('nan')

    return momu, momd, effective, hn_center, fup, fdn


class DBandApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Spin dependent d-band center calculation")
        self.geometry("680x420")
        try:
            self.iconphoto(False, tk.PhotoImage(file="logo.png"))
        except Exception:
            pass

        self._dos: Optional[DOSData] = None
        self._doscar_path: Optional[str] = None

        self._build_ui()

    def _build_ui(self):
        pad = {"padx": 8, "pady": 6}

        top = ttk.Frame(self)
        top.pack(fill=tk.X, **pad)

        title_lbl = ttk.Label(top, text="Spin dependent d-band center calculation", font=("Segoe UI", 14, "bold"))
        title_lbl.pack(side=tk.LEFT)

        # File chooser
        file_fr = ttk.LabelFrame(self, text="DOSCAR")
        file_fr.pack(fill=tk.X, **pad)
        self.path_var = tk.StringVar()
        path_entry = ttk.Entry(file_fr, textvariable=self.path_var)
        path_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(8, 4), pady=8)
        ttk.Button(file_fr, text="Browse...", command=self._choose_file).pack(side=tk.LEFT, padx=4, pady=8)
        ttk.Button(file_fr, text="Load", command=self._load_file).pack(side=tk.LEFT, padx=(4, 8), pady=8)

        # Atom range
        atoms_fr = ttk.LabelFrame(self, text="Atoms (1-based, inclusive)")
        atoms_fr.pack(fill=tk.X, **pad)
        ttk.Label(atoms_fr, text="from").pack(side=tk.LEFT, padx=(8, 4))
        self.from_var = tk.StringVar(value="1")
        ttk.Entry(atoms_fr, width=8, textvariable=self.from_var).pack(side=tk.LEFT)
        ttk.Label(atoms_fr, text="to").pack(side=tk.LEFT, padx=(12, 4))
        self.to_var = tk.StringVar(value="1")
        ttk.Entry(atoms_fr, width=8, textvariable=self.to_var).pack(side=tk.LEFT)
        self.natom_lbl = ttk.Label(atoms_fr, text="")
        self.natom_lbl.pack(side=tk.LEFT, padx=12)

        # Actions
        act_fr = ttk.Frame(self)
        act_fr.pack(fill=tk.X, **pad)
        self.compute_btn = ttk.Button(act_fr, text="Compute", command=self._compute, state=tk.DISABLED)
        self.compute_btn.pack(side=tk.LEFT, padx=8)

        # Results
        res_fr = ttk.LabelFrame(self, text="Results (eV)")
        res_fr.pack(fill=tk.BOTH, expand=True, **pad)

        self.momu_var = tk.StringVar(value="-")
        self.momd_var = tk.StringVar(value="-")
        self.eff_var = tk.StringVar(value="-")
        self.hn_var = tk.StringVar(value="-")
        self.fup_var = tk.StringVar(value="-")
        self.fdn_var = tk.StringVar(value="-")

        grid = ttk.Frame(res_fr)
        grid.pack(fill=tk.X, padx=8, pady=8)

        def row(r, label, var):
            ttk.Label(grid, text=label, width=34, anchor=tk.W).grid(row=r, column=0, sticky=tk.W, padx=(0, 8), pady=4)
            ttk.Label(grid, textvariable=var, width=18).grid(row=r, column=1, sticky=tk.W)

        row(0, "Spin-up d-band center (ε_d↑):", self.momu_var)
        row(1, "Spin-down d-band center (ε_d↓):", self.momd_var)
        row(2, "Effective d-band center (ε_eff):", self.eff_var)
        row(3, "Hammer–Nørskov d-band center:", self.hn_var)
        row(4, "Fractional occupation f↑:", self.fup_var)
        row(5, "Fractional occupation f↓:", self.fdn_var)

        # Footer with logo and credit (bottom-left)
        footer = ttk.Frame(self)
        footer.pack(fill=tk.X, side=tk.BOTTOM)
        try:
            self.footer_img = tk.PhotoImage(file="logo.png")
            ttk.Label(footer, image=self.footer_img).pack(side=tk.LEFT, padx=8, pady=6)
        except Exception:
            self.footer_img = None
        ttk.Label(footer, text="Created by Satadeep Bhattacharjee, IKST, Bangalore").pack(
            side=tk.LEFT, padx=8, pady=6
        )

        # Status bar (sits just above the footer)
        self.status = ttk.Label(self, relief=tk.SUNKEN, anchor=tk.W)
        self.status.pack(fill=tk.X, side=tk.BOTTOM)
        self._set_status("Ready. Load a DOSCAR (ISPIN=2, LORBIT=11).")

    def _set_status(self, msg: str):
        self.status.config(text=msg)

    def _choose_file(self):
        path = filedialog.askopenfilename(title="Select DOSCAR", filetypes=[("DOSCAR", "DOSCAR"), ("All files", "*.*")])
        if path:
            self.path_var.set(path)

    def _load_file(self):
        path = self.path_var.get().strip()
        if not path:
            messagebox.showwarning("No file", "Please select a DOSCAR file.")
            return
        try:
            self._set_status("Parsing DOSCAR...")
            dos = DOSCARParser().parse(path)
            self._dos = dos
            self._doscar_path = path
            self.natom_lbl.config(text=f"natom = {dos.natom}, NE = {dos.ne}")
            self.compute_btn.config(state=tk.NORMAL)
            self._set_status("Loaded. Set atom range and Compute.")
        except Exception as e:
            self._dos = None
            self.compute_btn.config(state=tk.DISABLED)
            messagebox.showerror("Parse error", str(e))
            self._set_status("Failed to load DOSCAR.")

    def _parse_range(self) -> Optional[Tuple[int, int]]:
        try:
            n1 = int(self.from_var.get())
            n2 = int(self.to_var.get())
            return n1, n2
        except Exception:
            messagebox.showwarning("Invalid atoms", "Please enter valid integers for 'from' and 'to'.")
            return None

    @staticmethod
    def _fmt(x: float) -> str:
        if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
            return "-"
        return f"{x: .4f}"

    def _compute(self):
        if self._dos is None:
            messagebox.showwarning("No data", "Please load a DOSCAR first.")
            return
        rng = self._parse_range()
        if not rng:
            return
        try:
            self._set_status("Computing d-band centers...")
            momu, momd, effective, hn_center, fup, fdn = compute_d_band_centers(self._dos, rng)
            self.momu_var.set(self._fmt(momu))
            self.momd_var.set(self._fmt(momd))
            self.eff_var.set(self._fmt(effective))
            self.hn_var.set(self._fmt(hn_center))
            self.fup_var.set(self._fmt(fup))
            self.fdn_var.set(self._fmt(fdn))
            self._set_status("Done.")
        except Exception as e:
            messagebox.showerror("Computation error", str(e))
            self._set_status("Failed.")


if __name__ == "__main__":
    app = DBandApp()
    app.mainloop()
