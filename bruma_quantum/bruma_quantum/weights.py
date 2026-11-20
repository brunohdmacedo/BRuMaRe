"""Leitor de pesos de projeção em formato 'long' (uma linha por spin, k, banda[, átomo]).

Formato nativo esperado (cabeçalho com nomes; separador: espaço ou TAB):

    # spin  k  band  atom  weight
    up      1  1     1     0.0123

Sinônimos aceitos nos nomes das colunas estão em SYN. 'atom' é opcional: sem ela,
'weight' é tratado como peso já somado na região. Linhas repetidas com a mesma
(spin, k, banda, átomo) são SOMADAS (útil quando há uma linha por orbital).
Se o .dat exportado pelo Bruma tiver outro layout, ajuste só este módulo.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

SYN = {
    "spin": {"spin", "s", "sp"},
    "k": {"k", "kpt", "kpoint", "k_index", "kindex", "ik"},
    "band": {"band", "banda", "n", "ib", "nband"},
    "atom": {"atom", "atomo", "átomo", "at", "ion", "ia"},
    "weight": {"weight", "w", "peso", "pesos", "proj", "projection"},
}


@dataclass
class WeightTensors:
    A: np.ndarray          # (ns, nk, nb, na) pesos
    covered: np.ndarray    # (ns, nk, nb) bool: há peso no arquivo
    n_atoms: int
    has_atoms: bool


def _canon(tok: str):
    t = tok.strip().lower()
    for key, names in SYN.items():
        if t in names:
            return key
    return None


def _split(line: str):
    return line.strip().lstrip("#").replace(",", " ").replace(";", " ").split()


def _spin_index(vals) -> np.ndarray:
    out = []
    for v in vals:
        s = str(v).strip().lower()
        if s in ("up", "u", "+"):
            out.append(0.0)
        elif s in ("down", "dn", "d", "-"):
            out.append(1.0)
        else:
            out.append(float(s))
    arr = np.array(out)
    if arr.min() >= 1:                 # 1/2 -> 0/1
        arr = arr - 1
    return arr.astype(int)


def read_weights(path, n_spins: int, n_kpoints: int, n_bands: int, base: int = 1) -> WeightTensors:
    lines = Path(path).read_text(errors="ignore").splitlines()
    cols, start = None, 0
    for i, raw in enumerate(lines):
        s = raw.strip()
        if not s:
            continue
        canon = [_canon(t) for t in _split(s)]
        if cols is None and sum(c is not None for c in canon) >= 3:
            cols = canon
            continue
        if s.startswith("#"):
            continue
        start = i
        break
    first = next(l for l in lines[start:] if l.strip() and not l.strip().startswith("#"))
    ncol = len(_split(first))
    if cols is None:
        cols = ["spin", "k", "band", "atom", "weight"] if ncol == 5 else ["spin", "k", "band", "weight"]
    if len(cols) != ncol:
        raise ValueError(f"Cabeçalho tem {len(cols)} colunas reconhecidas, mas os dados têm {ncol}.")
    if not {"spin", "k", "band", "weight"} <= set(cols):
        raise ValueError("Colunas obrigatórias: spin, k, band, weight (atom é opcional).")

    rows = [_split(l) for l in lines[start:] if l.strip() and not l.strip().startswith("#")]
    col = {name: [r[cols.index(name)] for r in rows] for name in set(c for c in cols if c)}
    spin = _spin_index(col["spin"])
    k = np.array(col["k"], dtype=float).astype(int) - base
    b = np.array(col["band"], dtype=float).astype(int) - base
    w = np.array(col["weight"], dtype=float)
    has_atoms = "atom" in col
    a = (np.array(col["atom"], dtype=float).astype(int) - base) if has_atoms else np.zeros(len(w), int)

    for name, arr, lim in (("spin", spin, n_spins), ("k", k, n_kpoints), ("band", b, n_bands)):
        if arr.min() < 0:
            raise ValueError(f"Índice de '{name}' negativo após base={base}: tente --base 0.")
        if arr.max() >= lim:
            raise ValueError(f"Índice de '{name}' = {arr.max() + 1} excede o .bands ({lim}).")
    if has_atoms and a.min() < 0:
        raise ValueError(f"Índice de 'atom' negativo após base={base}: tente --base 0.")
    na = int(a.max()) + 1
    A = np.zeros((n_spins, n_kpoints, n_bands, na))
    np.add.at(A, (spin, k, b, a), w)
    covered = np.zeros((n_spins, n_kpoints, n_bands), dtype=bool)
    covered[spin, k, b] = True
    return WeightTensors(A, covered, na, has_atoms)
