"""Linha de comando: python -m bruma_quantum <density|shells|hsites|orbital|fermi> ..."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .entanglement import entanglement_report, localization_entropy, binary_entropy, thermal_entropy
from .hsites import boltzmann_site_entropy, find_shells, resolve_center
from .info_density import analyze_density
from .io_castep import read_bands, read_cell, read_den_fmt
from .weights import read_weights


def _fmt(v) -> str:
    return f"{v:.8g}" if isinstance(v, (float, np.floating)) else str(v)


def write_dat(path, comments, cols, rows):
    with open(path, "w") as f:
        for c in comments:
            f.write(f"# {c}\n")
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(_fmt(v) for v in r) + "\n")
    print(f"  salvo: {path}")


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False))
    print(f"  salvo: {path}")


def _out(a, default_seed: str):
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir, (a.seed or default_seed)


# ---------------------------------------------------------------- density
def cmd_density(a):
    den = read_den_fmt(a.den)
    outdir, seed = _out(a, Path(a.den).stem)
    ref = read_den_fmt(a.ref) if a.ref else None
    center = None
    if a.cell and a.center:
        _, sym, frac = read_cell(a.cell)
        center = frac[resolve_center(sym, a.center)]
    rep = analyze_density(den, ref, center, a.radius, fisher=not a.no_fisher)
    if rep["aviso_densidade_de_diferenca"]:
        print("AVISO: há valores negativos: parece densidade de DIFERENÇA; Shannon exige densidade positiva.")
    for ch in ("total", "up", "down"):
        if ch in rep:
            r = rep[ch]
            print(f"[{ch}] N={r['n_electrons']:.2f}  S_dif={r['S_diferencial_nats']:.4f} nats  "
                  f"V_ef={r['V_efetivo_A3']:.1f} Å3 ({100 * r['V_efetivo_sobre_V_celula']:.1f}% da célula)  "
                  f"Onicescu={r['Onicescu_A-3']:.4g} Å-3" + (f"  Fisher={r['Fisher_A-2']:.4g} Å-2" if 'Fisher_A-2' in r else ""))
    if "JS_up_down_nats" in rep:
        print(f"JS(up,down)={rep['JS_up_down_nats']:.5f} nats  m={rep['momento_magnetico_muB']:.3f} muB")
    if "regiao" in rep:
        print("região:", {k: round(v, 5) if isinstance(v, float) else v for k, v in rep["regiao"]["total"].items()})
    if "vs_ref" in rep:
        print("vs referência:", rep["vs_ref"])
    write_json(outdir / f"{seed}_density_entropy.json", rep)


# ----------------------------------------------------------------- shells
def cmd_shells(a):
    lat, sym, frac = read_cell(a.cell)
    outdir, seed = _out(a, Path(a.cell).stem)
    shells = find_shells(lat, sym, frac, a.center, a.target, a.tol)
    print(f"Cascas de {a.target} em torno de {a.center} (tol={a.tol} Å):")
    for s in shells:
        print(f"  casca {s['shell']:2d}  d={s['distance_A']:.3f} Å  g={s['g']:3d}  (dispersão {s['spread_A']:.3f} Å)")
    write_dat(outdir / f"{seed}_shells.dat", [f"{a.target} em torno de {a.center}", f"tol = {a.tol} Å"],
              ["shell", "distance_A", "spread_A", "g"], [(s["shell"], s["distance_A"], s["spread_A"], s["g"]) for s in shells])
    tpl = outdir / f"{seed}_shells_template.csv"
    with open(tpl, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["shell", "distance_A", "g", "E_eV"])
        for s in shells:
            w.writerow([s["shell"], f"{s['distance_A']:.4f}", s["g"], ""])
    print(f"  template (preencha E_eV com a energia de remoção de H de cada casca): {tpl}")


# ----------------------------------------------------------------- hsites
def cmd_hsites(a):
    rows = list(csv.DictReader(open(a.shells)))
    if any(not r["E_eV"].strip() for r in rows):
        raise SystemExit("Preencha a coluna E_eV de todas as cascas no CSV.")
    g = [float(r["g"]) for r in rows]
    E = [float(r["E_eV"]) for r in rows]
    outdir, seed = _out(a, Path(a.shells).stem)
    temps = [float(t) for t in a.T.split(",")]
    res = [boltzmann_site_entropy(E, g, T) for T in temps]
    for r in res:
        print(f"T={r['T_K']:.0f} K  S_sítios={r['S_sitios_kB']:.4f} kB (máx {r['S_max_kB']:.4f})  "
              f"N_ef={r['N_sitios_efetivo']:.2f}  levantamento={r['levantamento_degenerescencia_kB']:.4f} kB")
    write_dat(outdir / f"{seed}_hsites_entropy.dat", ["entropia configuracional da vacância de H (Boltzmann)"],
              ["T_K", "S_familias_kB", "S_sitios_kB", "S_max_kB", "levantamento_kB", "N_sitios_efetivo", "S_sitios_J_mol_K"],
              [(r["T_K"], r["S_familias_kB"], r["S_sitios_kB"], r["S_max_kB"], r["levantamento_degenerescencia_kB"],
                r["N_sitios_efetivo"], r["S_sitios_J_mol_K"]) for r in res])
    write_dat(outdir / f"{seed}_hsites_prob.dat", ["probabilidade por família vs T"],
              ["T_K"] + [f"P_shell{r['shell']}" for r in rows], [[r["T_K"]] + r["P_familias"] for r in res])


# ---------------------------------------------------------------- orbital
def cmd_orbital(a):
    bands = read_bands(a.bands)
    outdir, seed = _out(a, Path(a.bands).stem)
    wt = read_weights(a.weights, bands.n_spins, bands.n_kpoints, bands.n_bands, a.base)
    if wt.has_atoms:
        if not a.region:
            raise SystemExit("O arquivo tem coluna de átomo: informe --region com os átomos (ex.: 1,5,25).")
        idx = [int(x) - a.base for x in a.region.split(",")]
        W = wt.A[..., idx].sum(axis=-1)
    else:
        W = wt.A[..., 0]
    rep = entanglement_report(bands, W, wt.covered, a.window, a.wmin)
    for msg in rep["aviso"]:
        print("AVISO:", msg)
    for sp in rep["por_spin"]:
        print(f"[{sp['spin']}] S_diag={sp['S_diag_nats']:.4f} nats  N_região={sp['N_eletrons_na_regiao']:.3f} e  "
              f"bandas no espaço ativo: {len(sp['espaco_ativo'])}")
    ea = rep["espaco_ativo_total"]
    print(f"S_diag total={rep['S_diag_total_nats']:.4f} nats (cota superior; compare com o puro) | "
          f"espaço ativo: {ea['spin_orbitais_qubits_JW']} spin-orbitais (≈ qubits JW), {ea['eletrons_ativos']} elétrons")
    names = ("up", "down") if bands.n_spins == 2 else ("unico",)
    for s in range(bands.n_spins):
        kw = bands.kweights[:, None]
        e_rel = (kw * bands.energies[s]).sum(0) - bands.fermi_eV[s]
        n_occ = bands.n_occupied(s)
        h = (kw * binary_entropy(W[s])).sum(0)
        w = (kw * W[s]).sum(0)
        cols = ["band", "E_minus_EF_eV", "occupied", "region_weight", "h_binary_nats"]
        data = [np.arange(1, bands.n_bands + 1), e_rel, (np.arange(bands.n_bands) < n_occ).astype(int), w, h]
        if wt.has_atoms:
            s_loc, n_eff, ipr, tot = localization_entropy(wt.A[s])
            ok = tot > 0
            for name, arr in (("S_loc_atoms", s_loc), ("N_eff_atoms", n_eff), ("IPR_atoms", ipr)):
                v = np.where(ok, arr, np.nan)
                data.append(np.nansum(kw * np.nan_to_num(v), axis=0))
            cols += ["S_loc_atoms", "N_eff_atoms", "IPR_atoms"]
        write_dat(outdir / f"{seed}_orbital_entropy_{names[s]}.dat",
                  [f"spin {names[s]}: médias ponderadas em k", f"E_F = {bands.fermi_eV[s]:.6f} eV"],
                  cols, list(zip(*data)))
    write_json(outdir / f"{seed}_entanglement.json", rep)


# ------------------------------------------------------------------ fermi
def cmd_fermi(a):
    bands = read_bands(a.bands)
    outdir, seed = _out(a, Path(a.bands).stem)
    rows = []
    for T in (float(t) for t in a.T.split(",")):
        S, N = thermal_entropy(bands, T)
        rows.append((T, S, S * 8.617333262e-5 * T, N))
        print(f"T={T:.0f} K  S={S:.4f} kB/célula  TS={S * 8.617333262e-5 * T:.5f} eV  N(T)={N:.3f}")
    write_dat(outdir / f"{seed}_thermal_entropy.dat", ["entropia de ocupação de Fermi-Dirac (autovalores KS fixos)"],
              ["T_K", "S_kB", "TS_eV", "N_eletrons_T"], rows)


def build_parser():
    p = argparse.ArgumentParser(prog="bruma_quantum", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--seed", default=None)
        sp.add_argument("--outdir", default=".")

    d = sub.add_parser("density", help="Shannon/Rényi/Onicescu/Fisher da densidade (.den_fmt)")
    d.add_argument("--den", required=True)
    d.add_argument("--ref", help=".den_fmt de referência (mesma célula e grade)")
    d.add_argument("--cell"); d.add_argument("--center"); d.add_argument("--radius", type=float)
    d.add_argument("--no-fisher", action="store_true")
    common(d); d.set_defaults(fn=cmd_density)

    s = sub.add_parser("shells", help="cascas/multiplicidade de H em torno do dopante (.cell)")
    s.add_argument("--cell", required=True); s.add_argument("--center", required=True)
    s.add_argument("--target", default="H"); s.add_argument("--tol", type=float, default=0.05)
    common(s); s.set_defaults(fn=cmd_shells)

    h = sub.add_parser("hsites", help="entropia configuracional (Boltzmann) a partir do CSV de cascas")
    h.add_argument("--shells", required=True); h.add_argument("--T", default="300,500,700")
    common(h); h.set_defaults(fn=cmd_hsites)

    o = sub.add_parser("orbital", help="emaranhamento (cota diagonal), localização e espaço ativo")
    o.add_argument("--bands", required=True); o.add_argument("--weights", required=True)
    o.add_argument("--region", help="átomos da região A (numeração do arquivo de pesos)")
    o.add_argument("--base", type=int, default=1, choices=(0, 1))
    o.add_argument("--window", type=float, default=2.0); o.add_argument("--wmin", type=float, default=0.2)
    common(o); o.set_defaults(fn=cmd_orbital)

    f = sub.add_parser("fermi", help="entropia térmica de ocupação (Fermi-Dirac)")
    f.add_argument("--bands", required=True); f.add_argument("--T", default="100,300,500,800")
    common(f); f.set_defaults(fn=cmd_fermi)
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    a.fn(a)
