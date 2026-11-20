# BrUMaRe × ferramentas de referência (mesmos arquivos)

| Código | Sistema | Grandeza | BrUMaRe | Referência | Ferramenta | |Δ| |
|---|---|---|---|---|---|---|
| VASP | Si (PBE) | gap (eV) | 0.6143 | 0.6143 | pymatgen | 1.1e-16 |
| VASP | Si (PBE) | E_F (eV) | 5.8561 | 5.8561 | pymatgen | 1.3e-05 |
| VASP | Si (PBE) | CBM − VBM (eV) | 0.6143 | 0.6143 | pymatgen | 1.1e-16 |
| VASP | Si (PBE) | m* elétron, lado 1 (m_e) | 0.9410 | 0.6771 | sumo | 2.6e-01 |
| VASP | Si (PBE) | m* elétron, lado 2 (m_e) | 0.9410 | 1.5400 | sumo | 6.0e-01 |
| VASP | Si (PBE) | m* elétron (exp. 0,92 m_e) | 0.9410 | 0.9200 | experimento | 2.1e-02 |
| VASP | Fe2O3 (LDA+U) | momento (μB) | 20.0000 | 20.0000 | pymatgen | 1.6e-06 |
| VASP | Fe2O3 (LDA+U) | gap spin ↑ (eV) | 4.0199 | 4.0199 | pymatgen | 0.0e+00 |
| VASP | Fe2O3 (LDA+U) | gap spin ↓ (eV) | 2.1923 | 2.1923 | pymatgen | 0.0e+00 |
| CASTEP | MgH2:Fe (162 at.) | gap total (eV) | 1.3940 | 0.0000 | sumo | 1.4e+00 |
| CASTEP | MgH2:Fe (162 at.) | E_F (eV) | 1.8368 | 1.8368 | sumo | 4.2e-06 |
| CASTEP | MgH2:Fe (162 at.) | gap, contagem de elétrons nos autovalores do sumo (eV) | 1.3940 | 1.3941 | sumo | 7.4e-05 |
| SIESTA | SrTiO3 (LDA) | gap total (eV) | 2.2500 | 2.2500 | sisl | 0.0e+00 |
| SIESTA | Bi2Se3 | gap (eV) | 0.0418 | 0.0418 | sisl | 6.0e-16 |
| ORCA | divinylbenzene | HOMO (eV) | -4.0652 | -4.0652 | cclib | 1.8e-05 |
| ORCA | divinylbenzene | LUMO (eV) | 1.1197 | 1.1197 | cclib | 5.9e-06 |
| ORCA | divinylbenzene | gap HOMO–LUMO (eV) | 5.1849 | 5.1849 | cclib | 1.2e-05 |
| CHGCAR | MgH2:Fe CHGCAR (162 at.) | Bader on-grid: |Δq| médio (e) | 0.0000 | 0.0127 | pybader | 1.3e-02 |
| CHGCAR | MgH2:Fe CHGCAR (162 at.) | Bader: |Δq| médio (e) | 0.0000 | 0.0008 | pybader | 8.0e-04 |
| CHGCAR | MgH2:Fe CHGCAR (162 at.) | Bader: Σq (e) | 651.0000 | 651.0000 | pybader | 2.3e-13 |
| CHGCAR | MgH2:Fe CHGCAR (162 at.) | Bader: |Δq| máximo (e) | 0.0000 | 0.0257 | pybader | 2.6e-02 |

## Tempos de execução

| Caso | Tempo (s) |
|---|---|
| vasp_Si_bands | 3.68 |
| vasp_Fe2O3_ldau | 4.83 |
| mgh2_Fe | 4.49 |
| SrTiO3 | 7.55 |
| Bi2Se3 | 1.52 |
| dvb_sp.out | 1.63 |
| Bader on-grid (BrUMaRe) | 0.13 |
| Bader refinado (BrUMaRe) | 0.66 |
| Bader (pybader) | 1.18 |
