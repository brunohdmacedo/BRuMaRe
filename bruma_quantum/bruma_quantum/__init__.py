"""bruma_quantum: entropia, emaranhamento e descritores de informação quântica a partir de DFT (CASTEP).

API pública pensada para integração no BruMaRe: cada função recebe/retorna
estruturas simples (numpy/dict) e não faz I/O além dos leitores em io_castep.
"""
__version__ = "0.1.0"

from .io_castep import BandsData, DensityData, read_bands, read_cell, read_den_fmt
from .info_density import (analyze_density, fisher_information, jensen_shannon,
                           kl_divergence, shannon_report, sphere_mask)
from .entanglement import (binary_entropy, diagonal_entropy, entanglement_report,
                           localization_entropy, peschel_entropy, suggest_active_space,
                           thermal_entropy)
from .hsites import boltzmann_site_entropy, find_shells
from .weights import WeightTensors, read_weights
