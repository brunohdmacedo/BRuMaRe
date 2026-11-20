"""Entropias de von Neumann (emaranhamento) e descritores para espaço ativo.

IMPORTANTE: o determinante de Slater de Kohn-Sham é um estado puro de um sistema
auxiliar de férmions livres. Aqui S mede o emaranhamento ENTRE UMA REGIÃO E O
RESTO NESSE SISTEMA AUXILIAR (fórmula de Peschel), não o emaranhamento do
sistema real interagente. Use para tendências relativas (defeito x puro).
"""
from __future__ import annotations

import numpy as np

from .io_castep import BandsData

KB_EV = 8.617333262e-5


def binary_entropy(x) -> np.ndarray:
    """h(x) = -x ln x - (1-x) ln(1-x), com h(0)=h(1)=0."""
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    a = np.clip(x, 1e-300, 1.0)
    b = np.clip(1.0 - x, 1e-300, 1.0)
    return -(a * np.log(a)) - (b * np.log(b))


def peschel_entropy(C: np.ndarray):
    """S exata de uma região para um determinante de Slater.

    C_ij = <psi_i|P_A|psi_j> sobre os orbitais OCUPADOS (matriz de overlap restrita
    à região A). S = sum h(nu), nu = autovalores de C. Retorna (S, nu).
    """
    C = np.asarray(C)
    nu = np.linalg.eigvalsh(0.5 * (C + C.conj().T))
    return float(binary_entropy(nu).sum()), nu


def diagonal_entropy(W: np.ndarray, kweights: np.ndarray, n_occ: int):
    """Aproximação diagonal a partir de pesos de projeção w_nk (nk, nb).

    S ~ sum_k omega_k sum_n h(w_nk), bandas ocupadas. Como a diagonal de C é
    majorizada pelos autovalores, isto é uma COTA SUPERIOR da entropia exata
    (por k, com projeção completa). Para estados estendidos (w~fração da região)
    a cota é frouxa: compare SEMPRE com o sistema puro e olhe s_n por banda.
    """
    h = binary_entropy(W[:, :n_occ])
    s_band = (kweights[:, None] * h).sum(axis=0)
    return float(s_band.sum()), s_band


def localization_entropy(A: np.ndarray):
    """Entropia de Shannon da distribuição de cada estado entre átomos.

    A: (nk, nb, na) pesos por átomo. Retorna S_loc, N_eff=exp(S_loc), IPR e o
    peso total (nk, nb). N_eff ~ número efetivo de átomos que o estado ocupa.
    """
    tot = A.sum(axis=-1)
    p = np.divide(A, tot[..., None], out=np.zeros_like(A), where=tot[..., None] > 0)
    plogp = np.zeros_like(p)
    m = p > 0
    plogp[m] = p[m] * np.log(p[m])
    s = -plogp.sum(axis=-1)
    return s, np.exp(s), (p ** 2).sum(axis=-1), tot


def suggest_active_space(en_eV: np.ndarray, fermi_eV: float, W: np.ndarray, kweights: np.ndarray,
                         n_occ: int, window_eV: float = 2.0, wmin: float = 0.2) -> list:
    """Triagem heurística de bandas candidatas a espaço ativo (cluster/embedding).

    Seleciona bandas com |<E>_k - E_F| <= janela e peso médio na região >= wmin.
    Não substitui AVAS/DMET/entropia de orbital com matriz completa.
    """
    e_rel = (kweights[:, None] * en_eV).sum(axis=0) - fermi_eV
    w = (kweights[:, None] * W).sum(axis=0)
    sel = np.where((np.abs(e_rel) <= window_eV) & (w >= wmin))[0]
    return [{"band": int(b) + 1, "E_minus_EF_eV": float(e_rel[b]), "weight": float(w[b]),
             "occupied": bool(b < n_occ)} for b in sel]


def entanglement_report(bands: BandsData, W: np.ndarray, covered: np.ndarray | None = None,
                        window_eV: float = 2.0, wmin: float = 0.2) -> dict:
    """W: (ns, nk, nb) pesos da região. covered: máscara (ns, nk, nb) de pesos presentes."""
    out = {"por_spin": [], "aviso": []}
    s_total = n_reg_total = 0.0
    spin_orbitals = n_elec = 0
    for s in range(bands.n_spins):
        n_occ = bands.n_occupied(s)
        S, s_band = diagonal_entropy(W[s], bands.kweights, n_occ)
        n_reg = float((bands.kweights[:, None] * np.clip(W[s][:, :n_occ], 0, 1)).sum())
        miss = 0 if covered is None else int((~covered[s][:, :n_occ]).any(axis=0).sum())
        if miss:
            out["aviso"].append(f"spin {s + 1}: {miss} banda(s) ocupada(s) sem peso no arquivo; S subestimada.")
        act = suggest_active_space(bands.energies[s], bands.fermi_eV[s], W[s], bands.kweights,
                                   n_occ, window_eV, wmin)
        mult = 2 if bands.n_spins == 1 else 1     # sem polarização: os dois spins são idênticos
        out["por_spin"].append({"spin": "unico" if bands.n_spins == 1 else ("up", "down")[s],
                                "S_diag_nats": S, "N_eletrons_na_regiao": n_reg,
                                "n_bandas_ocupadas": n_occ, "bandas_sem_peso": miss,
                                "espaco_ativo": act})
        s_total += mult * S
        n_reg_total += mult * n_reg
        spin_orbitals += mult * len(act)
        n_elec += mult * sum(a["occupied"] for a in act)
    out.update({"S_diag_total_nats": s_total, "N_eletrons_regiao_total": n_reg_total,
                "espaco_ativo_total": {"spin_orbitais_qubits_JW": int(spin_orbitals),
                                       "eletrons_ativos": int(n_elec),
                                       "janela_eV": window_eV, "peso_minimo": wmin}})
    return out


# ---------------------------------------------------------------- térmica
def fermi_dirac(E: np.ndarray, mu: float, T: float) -> np.ndarray:
    return 0.5 * (1.0 - np.tanh((E - mu) / (2.0 * KB_EV * T)))


def thermal_entropy(bands: BandsData, T: float):
    """Entropia de ocupação de Fermi-Dirac dos autovalores KS congelados (em k_B por célula).

    Retorna (S/k_B, N(T)). Mantém mu = E_F do .bands; confira N(T) ~ N_elétrons.
    """
    S = N = 0.0
    for s in range(bands.n_spins):
        f = fermi_dirac(bands.energies[s], bands.fermi_eV[s], T)
        S += float((bands.kweights[:, None] * binary_entropy(f)).sum())
        N += float((bands.kweights[:, None] * f).sum())
    if bands.n_spins == 1:
        S, N = 2 * S, 2 * N
    return S, N
