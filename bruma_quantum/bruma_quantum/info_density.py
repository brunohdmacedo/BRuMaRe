"""Descritores de teoria da informação (ITA) da densidade eletrônica.

Todas as grandezas usam a densidade normalizada p_i = rho_i dV / N_e sobre a grade.
Unidades: nats; volumes em Å³; Fisher em Å⁻².
"""
from __future__ import annotations

import numpy as np

from .io_castep import DensityData


def _prob(rho: np.ndarray, dV: float):
    r = np.clip(rho, 0.0, None)
    n = float(r.sum() * dV)
    if n <= 0:
        raise ValueError("Densidade sem carga positiva.")
    return r * dV / n, n


def kl_divergence(p: np.ndarray, q: np.ndarray, eps: float = 1e-300) -> float:
    m = p > 0
    return float(np.sum(p[m] * np.log(p[m] / np.maximum(q[m], eps))))


def jensen_shannon(p: np.ndarray, q: np.ndarray) -> float:
    m = 0.5 * (p + q)
    return 0.5 * kl_divergence(p, m) + 0.5 * kl_divergence(q, m)


def shannon_report(rho: np.ndarray, volume: float) -> dict:
    """Shannon (discreta e diferencial), Rényi-2 e energia de Onicescu."""
    dV = volume / rho.size
    p, n_el = _prob(rho, dV)
    m = p > 0
    s_disc = float(-(p[m] * np.log(p[m])).sum())
    s_dif = s_disc + np.log(dV)                  # S = -sum p ln(p/dV)
    sum_p2 = float((p ** 2).sum())
    return {
        "n_electrons": n_el,
        "S_discreta_nats": s_disc,
        "S_discreta_normalizada": s_disc / float(np.log(rho.size)),
        "S_diferencial_nats": float(s_dif),
        "V_efetivo_A3": float(np.exp(s_dif)),
        "V_efetivo_sobre_V_celula": float(np.exp(s_dif) / volume),
        "Renyi2_nats": float(-np.log(sum_p2) + np.log(dV)),
        "Onicescu_A-3": sum_p2 / dV,
    }


def fisher_information(rho: np.ndarray, lattice: np.ndarray) -> float:
    """I = 4 ∫|∇ sqrt(rho~)|² dr via FFT (periódico, rede geral). Unidade: Å⁻²."""
    nx, ny, nz = rho.shape
    volume = abs(float(np.linalg.det(lattice)))
    dV = volume / rho.size
    p, _ = _prob(rho, dV)
    psi = np.sqrt(p / dV)                         # sqrt da densidade normalizada
    F = np.fft.fftn(psi)
    B = np.linalg.inv(lattice).T                  # linhas b_i com b_i·a_j = delta_ij
    grads = [np.zeros(rho.shape) for _ in range(3)]
    for i, n in enumerate((nx, ny, nz)):
        shape = [1, 1, 1]
        shape[i] = -1
        f = np.fft.fftfreq(n, d=1.0 / n).reshape(shape)       # ciclos por célula
        d_si = np.fft.ifftn(2j * np.pi * f * F).real           # d psi / d s_i
        for c in range(3):
            grads[c] += B[i, c] * d_si
    return float(4.0 * sum((g ** 2).sum() for g in grads) * dV)


def sphere_mask(lattice: np.ndarray, dims, center_frac, radius: float) -> np.ndarray:
    """Máscara booleana dos pontos da grade a até 'radius' Å do centro (imagem mínima)."""
    nx, ny, nz = dims
    d = np.empty((nx, ny, nz, 3))
    d[..., 0] = (np.arange(nx) / nx - center_frac[0])[:, None, None]
    d[..., 1] = (np.arange(ny) / ny - center_frac[1])[None, :, None]
    d[..., 2] = (np.arange(nz) / nz - center_frac[2])[None, None, :]
    d -= np.round(d)
    return np.linalg.norm(d @ lattice, axis=-1) <= radius


def region_report(rho: np.ndarray, mask: np.ndarray, dV: float, ref: np.ndarray | None = None) -> dict:
    r = np.clip(rho, 0, None)
    n_r = float(r[mask].sum() * dV)
    p = r[mask] * dV / n_r
    m = p > 0
    s = float(-(p[m] * np.log(p[m] / dV)).sum())
    out = {"n_pontos": int(mask.sum()), "N_eletrons": n_r,
           "S_diferencial_nats": s, "V_efetivo_A3": float(np.exp(s))}
    if ref is not None:
        rr = np.clip(ref, 0, None)
        n_ref = float(rr[mask].sum() * dV)
        q = rr[mask] * dV / n_ref
        mq = q > 0
        s_ref = float(-(q[mq] * np.log(q[mq] / dV)).sum())
        out.update({"N_eletrons_ref": n_ref, "dN_regiao": n_r - n_ref,
                    "S_ref_nats": s_ref, "dS_nats": s - s_ref,
                    "KL_defeito_ref": kl_divergence(p, q),
                    "JS_defeito_ref": jensen_shannon(p, q)})
    return out


def analyze_density(den: DensityData, ref: DensityData | None = None,
                    center_frac=None, radius: float | None = None, fisher: bool = True) -> dict:
    """Relatório completo: total, spin up/down, região esférica e comparação com referência."""
    out = {"grade": list(den.dims), "volume_A3": den.volume,
           "aviso_densidade_de_diferenca": den.looks_like_difference}
    channels = {"total": den.rho}
    if den.spin is not None:
        channels["up"] = 0.5 * (den.rho + den.spin)
        channels["down"] = 0.5 * (den.rho - den.spin)
    for name, ch in channels.items():
        rep = shannon_report(ch, den.volume)
        if fisher:
            rep["Fisher_A-2"] = fisher_information(ch, den.lattice)
        out[name] = rep
    dV = den.volume / den.rho.size
    if den.spin is not None:
        pu, _ = _prob(channels["up"], dV)
        pd, _ = _prob(channels["down"], dV)
        out["JS_up_down_nats"] = jensen_shannon(pu, pd)
        out["momento_magnetico_muB"] = float(den.spin.sum() * dV)

    mask = None
    if center_frac is not None and radius:
        mask = sphere_mask(den.lattice, den.dims, center_frac, radius)
        out["regiao"] = {"raio_A": radius, "total": region_report(
            den.rho, mask, dV, ref.rho if ref is not None else None)}

    if ref is not None:
        if ref.dims != den.dims or not np.allclose(ref.lattice, den.lattice, atol=1e-3):
            out["vs_ref"] = {"erro": "grade/célula diferentes: KL/JS exigem a mesma grade."}
        else:
            p, _ = _prob(den.rho, dV)
            q, _ = _prob(ref.rho, dV)
            s_ref = shannon_report(ref.rho, ref.volume)
            out["vs_ref"] = {"dS_diferencial_nats": out["total"]["S_diferencial_nats"] - s_ref["S_diferencial_nats"],
                             "KL_defeito_ref": kl_divergence(p, q),
                             "JS_defeito_ref": jensen_shannon(p, q)}
    return out
