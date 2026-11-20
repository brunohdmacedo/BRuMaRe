"""Degenerescência de sítios (ex.: H em torno de um dopante) e entropia configuracional.

1) find_shells: agrupa átomos-alvo em cascas por distância ao dopante (imagem
   mínima) e devolve a multiplicidade g de cada casca.
2) boltzmann_site_entropy: com a energia E_i de cada família (ex.: remoção de H),
   P_i = g_i exp(-E_i/kT)/Z e entropias de Shannon sobre famílias e sobre sítios.

Agrupar só por distância é uma triagem: sítios a mesma distância podem ser
inequivalentes por simetria; confira com a simetria da supercélula (ex.: spglib).
"""
from __future__ import annotations

import numpy as np

KB_EV = 8.617333262e-5
R_J_MOL_K = 8.314462618


def min_image_distances(lattice: np.ndarray, frac_center: np.ndarray, frac_pts: np.ndarray) -> np.ndarray:
    shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    d = frac_pts - frac_center
    d = d - np.round(d)
    cart = (d[:, None, :] + shifts[None, :, :]) @ lattice
    return np.linalg.norm(cart, axis=-1).min(axis=1)


def resolve_center(symbols: list, center) -> int:
    """center: símbolo único (ex.: 'Ni') ou índice 1-based. Retorna índice 0-based."""
    try:
        idx = int(center) - 1
        if not 0 <= idx < len(symbols):
            raise ValueError(f"Índice {center} fora do intervalo 1..{len(symbols)}.")
        return idx
    except (TypeError, ValueError) as exc:
        if str(center).isdigit():
            raise exc
    hits = [i for i, s in enumerate(symbols) if s == center]
    if len(hits) != 1:
        raise ValueError(f"'{center}' aparece {len(hits)} vezes: informe o índice (1-based).")
    return hits[0]


def find_shells(lattice, symbols, frac, center, target: str = "H", tol: float = 0.05) -> list:
    c = resolve_center(symbols, center)
    tgt = [i for i, s in enumerate(symbols) if s == target and i != c]
    if not tgt:
        raise ValueError(f"Nenhum átomo '{target}' encontrado.")
    dist = min_image_distances(np.asarray(lattice), np.asarray(frac)[c], np.asarray(frac)[tgt])
    order = np.argsort(dist)
    shells, cur = [], [order[0]]
    for a, b in zip(order[:-1], order[1:]):
        if dist[b] - dist[a] > tol:
            shells.append(cur)
            cur = []
        cur.append(b)
    shells.append(cur)
    return [{"shell": n + 1, "distance_A": float(dist[s].mean()), "spread_A": float(np.ptp(dist[s])),
             "g": len(s), "atoms_1based": sorted(int(tgt[j]) + 1 for j in s)}
            for n, s in enumerate(shells)]


def boltzmann_site_entropy(E_eV, g, T: float) -> dict:
    """E_eV: energia de cada família; g: multiplicidades; T em K (>0)."""
    E, g = np.asarray(E_eV, dtype=float), np.asarray(g, dtype=float)
    w = g * np.exp(-(E - E.min()) / (KB_EV * T))
    P = w / w.sum()
    m = P > 0
    s_fam = float(-(P[m] * np.log(P[m])).sum())
    s_sites = float(-(P[m] * np.log(P[m] / g[m])).sum())     # distribuição sobre sítios individuais
    s_max = float(np.log(g.sum()))
    return {"T_K": float(T), "P_familias": P.tolist(), "S_familias_kB": s_fam, "S_sitios_kB": s_sites,
            "S_max_kB": s_max, "levantamento_degenerescencia_kB": s_max - s_sites,
            "N_sitios_efetivo": float(np.exp(s_sites)), "S_sitios_J_mol_K": s_sites * R_J_MOL_K}
