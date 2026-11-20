"""Leitura dos arquivos do CASTEP usados pelo bruma_quantum (somente numpy).

Formatos cobertos: .bands (autovalores, k-pontos, E_F por spin), .cell (rede e
posições) e .den_fmt (densidade de carga/spin formatada).
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HARTREE_TO_EV = 27.211386245988


# --------------------------------------------------------------------- .bands
@dataclass
class BandsData:
    n_kpoints: int
    n_spins: int
    n_electrons: list          # por spin (1 spin: total de elétrons)
    fermi_eV: list             # por spin
    kpoints: np.ndarray        # (nk, 3)
    kweights: np.ndarray       # (nk,), soma 1
    energies: list             # por spin: array (nk, nb) em eV

    def n_occupied(self, spin: int) -> int:
        """Número de bandas ocupadas no canal de spin (arredondado)."""
        if self.n_spins == 1:
            return int(round(self.n_electrons[0] / 2.0))
        return int(round(self.n_electrons[spin]))

    @property
    def n_bands(self) -> int:
        return int(self.energies[0].shape[1])


def read_bands(path) -> BandsData:
    lines = Path(path).read_text(errors="ignore").splitlines()
    n_k = n_s = None
    n_el, n_eig, fermi = [], [], []
    i = 0
    while i < len(lines) and not lines[i].strip().startswith("K-point"):
        s = lines[i].strip()
        ns = n_s or 1
        if s.startswith("Number of k-points"):
            n_k = int(s.split()[-1])
        elif s.startswith("Number of spin components"):
            n_s = int(s.split()[-1])
        elif s.startswith("Number of electrons"):
            n_el = [float(x) for x in s.split()[-ns:]]
        elif s.startswith("Number of eigenvalues"):
            n_eig = [int(x) for x in s.split()[-ns:]]
        elif "Fermi energ" in s:
            fermi = [float(x) * HARTREE_TO_EV for x in s.split()[-ns:]]
        i += 1
    if n_k is None or n_s is None or not n_eig:
        raise ValueError("Cabeçalho do .bands não reconhecido.")

    kpts = np.zeros((n_k, 3))
    kw = np.zeros(n_k)
    en = [np.zeros((n_k, n_eig[s])) for s in range(n_s)]
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith("K-point"):
            p = s.split()
            ik = int(p[1]) - 1
            kpts[ik] = [float(x) for x in p[2:5]]
            kw[ik] = float(p[5])
            i += 1
            for sp in range(n_s):
                i += 1  # linha "Spin component"
                en[sp][ik] = np.array(lines[i:i + n_eig[sp]], dtype=float) * HARTREE_TO_EV
                i += n_eig[sp]
        else:
            i += 1
    kw = kw / kw.sum()
    return BandsData(n_k, n_s, n_el, fermi, kpts, kw, en)


# ---------------------------------------------------------------------- .cell
def _blocks(text: str) -> dict:
    out, cur, name = {}, None, None
    for raw in text.splitlines():
        line = re.split(r"[!#]", raw, maxsplit=1)[0].strip()
        if not line:
            continue
        low = line.lower()
        if low.startswith("%block"):
            name, cur = low.split()[1], []
            continue
        if low.startswith("%endblock"):
            out[name] = cur
            name, cur = None, None
            continue
        if cur is not None:
            cur.append(line)
    return out


def _abc_to_cart(a, b, c, al, be, ga) -> np.ndarray:
    al, be, ga = np.radians([al, be, ga])
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz = np.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return np.array([[a, 0, 0], [b * np.cos(ga), b * np.sin(ga), 0], [cx, cy, cz]])


def _floats(line: str):
    try:
        return [float(x) for x in line.split()]
    except ValueError:
        return None


def read_cell(path):
    """Retorna (rede 3x3 em Å, símbolos, coordenadas fracionárias (n,3))."""
    blocks = _blocks(Path(path).read_text(errors="ignore"))
    if "lattice_cart" in blocks:
        rows = [v for v in (_floats(l) for l in blocks["lattice_cart"]) if v and len(v) == 3]
        lattice = np.array(rows[:3])
    elif "lattice_abc" in blocks:
        rows = [v for v in (_floats(l) for l in blocks["lattice_abc"]) if v and len(v) == 3]
        lattice = _abc_to_cart(*rows[0], *rows[1])
    else:
        raise ValueError("Nem lattice_cart nem lattice_abc encontrados no .cell.")

    sym, pos = [], []
    key = "positions_frac" if "positions_frac" in blocks else "positions_abs"
    if key not in blocks:
        raise ValueError("Nem positions_frac nem positions_abs encontrados no .cell.")
    for line in blocks[key]:
        p = line.split()
        if len(p) < 4 or _floats(" ".join(p[1:4])) is None:
            continue  # linha de unidade (ex.: ANG)
        sym.append(p[0].split(":")[0])
        pos.append([float(x) for x in p[1:4]])
    pos = np.array(pos)
    if key == "positions_abs":
        pos = pos @ np.linalg.inv(lattice)
    return lattice, sym, pos


# ------------------------------------------------------------------- .den_fmt
@dataclass
class DensityData:
    lattice: np.ndarray            # (3,3) em Å
    dims: tuple
    volume: float                  # Å³
    rho: np.ndarray                # e/Å³ (total)
    spin: np.ndarray | None        # e/Å³ (rho_up - rho_down) ou None
    looks_like_difference: bool    # True se há valores negativos relevantes


def read_den_fmt(path) -> DensityData:
    text = Path(path).read_text(errors="ignore")
    head, sep, body = text.partition("END header")
    if not sep:
        raise ValueError("'END header' não encontrado: formato inesperado.")
    body = body.split("\n", 1)[1]          # descarta o resto da linha do 'END header'
    lattice, dims = None, None
    lines = head.splitlines()
    for i, line in enumerate(lines):
        if "Real Lattice" in line:
            lattice = np.array([[float(x) for x in lines[i + j].split()[:3]] for j in (1, 2, 3)])
        m = re.match(r"^\s*(\d+)\s+(\d+)\s+(\d+)\b", line)
        if m:
            dims = tuple(int(x) for x in m.groups())
    if lattice is None or dims is None:
        raise ValueError("Rede ou dimensões da grade não encontradas no cabeçalho.")
    ncols = len(body.strip().split("\n", 1)[0].split())
    arr = np.fromstring(body, sep=" ", dtype=float).reshape(-1, ncols)
    if arr.shape[0] != int(np.prod(dims)):
        raise ValueError(f"{arr.shape[0]} pontos lidos, esperado {int(np.prod(dims))}.")
    volume = abs(float(np.linalg.det(lattice)))
    idx = arr[:, :3].astype(int) - 1
    rho = np.zeros(dims)
    rho[idx[:, 0], idx[:, 1], idx[:, 2]] = arr[:, 3] / volume   # valores brutos = rho * V
    spin = None
    if ncols >= 5:
        spin = np.zeros(dims)
        spin[idx[:, 0], idx[:, 1], idx[:, 2]] = arr[:, 4] / volume
    # Ruído numérico dá poucos pontos levemente negativos (fração ~1e-6). Densidade de
    # DIFERENÇA tem parte negativa comparável à positiva: usa-se a fração integrada.
    neg_frac = float(-rho[rho < 0].sum() / max(np.abs(rho).sum(), 1e-30))
    return DensityData(lattice, dims, volume, rho, spin, neg_frac > 0.05)
