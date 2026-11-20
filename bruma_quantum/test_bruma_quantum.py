"""Validações com respostas conhecidas. Rode: python test_bruma_quantum.py"""
import tempfile
from pathlib import Path

import numpy as np

from bruma_quantum import (binary_entropy, boltzmann_site_entropy, diagonal_entropy, find_shells,
                           fisher_information, jensen_shannon, localization_entropy, peschel_entropy,
                           read_weights, shannon_report)


def test_gaussian_shannon_fisher():
    L, n, sigma = 16.0, 80, 1.0
    x = (np.arange(n) / n - 0.5) * L
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    rho = 10.0 * np.exp(-(X**2 + Y**2 + Z**2) / (2 * sigma**2)) / ((2 * np.pi) ** 1.5 * sigma**3)
    lattice = np.eye(3) * L
    rep = shannon_report(rho, L**3)
    s_ref = 1.5 * np.log(2 * np.pi * np.e * sigma**2)        # Shannon de uma Gaussiana 3D
    assert abs(rep["S_diferencial_nats"] - s_ref) < 0.02, rep["S_diferencial_nats"]
    assert abs(rep["n_electrons"] - 10.0) < 1e-3
    assert abs(fisher_information(rho, lattice) - 3.0 / sigma**2) < 0.02       # I = 3/sigma^2
    # Onicescu analítico: ∫rho~² = 1/((4 pi sigma²)^(3/2))
    assert abs(rep["Onicescu_A-3"] - (4 * np.pi * sigma**2) ** -1.5) < 1e-3


def test_fisher_nonorthogonal_cell():
    """Mesma Gaussiana em célula oblíqua deve dar o mesmo Fisher (testa a métrica)."""
    n, sigma = 72, 1.0
    A = np.array([[14.0, 0, 0], [3.0, 14.0, 0], [2.0, 2.0, 14.0]])
    s = (np.arange(n) / n)[:, None] * np.ones((1, n))
    f = np.stack(np.meshgrid(*(3 * [np.arange(n) / n]), indexing="ij"), axis=-1)   # (n,n,n,3)
    f = f - np.round(f - 0.5)                                                      # centro no meio
    r = (f - 0.5) @ A
    rho = np.exp(-(r**2).sum(-1) / (2 * sigma**2))
    assert abs(fisher_information(rho, A) - 3.0) < 0.05


def test_peschel_two_sites_ln2():
    psi = np.array([1.0, 1.0]) / np.sqrt(2)                  # orbital ligante
    P = np.diag([1.0, 0.0])
    C = np.array([[psi @ P @ psi]])
    S, nu = peschel_entropy(C)
    assert abs(S - np.log(2)) < 1e-12 and abs(nu[0] - 0.5) < 1e-12


def test_peschel_vs_diagonal_bound_chain():
    L, ell = 20, 10
    H = -(np.eye(L, k=1) + np.eye(L, k=-1))
    _, V = np.linalg.eigh(H)
    occ = V[:, : L // 2]                                      # meio preenchimento
    C = occ[:ell].T @ occ[:ell]
    S_exact, _ = peschel_entropy(C)
    S_diag, _ = diagonal_entropy(np.diag(C)[None, :], np.array([1.0]), L // 2)
    assert 0 < S_exact < S_diag                               # cota superior (e frouxa p/ estados estendidos)
    print(f"   cadeia: S_exata={S_exact:.3f}  S_diag={S_diag:.3f}")


def test_boltzmann_sites():
    r = boltzmann_site_entropy([0.0, 0.0], [1, 3], 300.0)    # energias iguais -> uniforme em 4 sítios
    assert abs(r["S_sitios_kB"] - np.log(4)) < 1e-12 and abs(r["levantamento_degenerescencia_kB"]) < 1e-12
    r2 = boltzmann_site_entropy([0.0, 1.0], [1, 3], 300.0)   # 1 eV de diferença -> só a família 1
    assert r2["S_sitios_kB"] < 1e-6 and abs(r2["P_familias"][0] - 1) < 1e-6


def test_shells_multiplicity():
    L = 10.0
    lat = np.eye(3) * L
    sym, frac = ["Ni"], [[0.0, 0.0, 0.0]]
    pts = []
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            for sz in (-1, 0, 1):
                if (sx, sy, sz) != (0, 0, 0):
                    pts.append(np.array([sx, sy, sz]) * 0.2)
    for p in pts:
        sym.append("H"); frac.append(list(p % 1.0))
    sh = find_shells(lat, sym, np.array(frac), "Ni", "H", tol=0.05)
    # 6 a 2.0 Å, 12 a 2.83 Å, 8 a 3.46 Å
    assert [s["g"] for s in sh] == [6, 12, 8], [s["g"] for s in sh]


def test_localization_entropy():
    A = np.zeros((1, 2, 4))
    A[0, 0] = [1, 0, 0, 0]                                    # localizado em 1 átomo
    A[0, 1] = [1, 1, 1, 1]                                    # espalhado em 4
    s, n_eff, ipr, _ = localization_entropy(A)
    assert abs(n_eff[0, 0] - 1) < 1e-12 and abs(n_eff[0, 1] - 4) < 1e-12
    assert abs(ipr[0, 1] - 0.25) < 1e-12


def test_weights_reader_roundtrip():
    txt = "# spin  k  band  atom  weight\nup 1 1 1 0.25\nup 1 1 2 0.75\ndown 2 2 1 1.0\n"
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "w.dat"
        p.write_text(txt)
        wt = read_weights(p, n_spins=2, n_kpoints=2, n_bands=2, base=1)
    assert wt.A.shape == (2, 2, 2, 2) and wt.has_atoms
    assert abs(wt.A[0, 0, 0].sum() - 1.0) < 1e-12 and wt.A[1, 1, 1, 0] == 1.0
    assert wt.covered[0, 0, 0] and not wt.covered[0, 1, 1]


def test_binary_entropy_limits():
    assert binary_entropy(0.0) == 0 or binary_entropy(0.0) < 1e-290
    assert abs(float(binary_entropy(0.5)) - np.log(2)) < 1e-12
    assert jensen_shannon(np.array([1.0, 0]), np.array([0, 1.0])) > 0.69   # ~ ln 2


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"OK   {name}")
            except Exception as exc:                              # noqa: BLE001
                fails += 1
                print(f"FAIL {name}: {exc!r}")
    raise SystemExit(1 if fails else 0)
