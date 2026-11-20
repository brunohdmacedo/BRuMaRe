**Title:** CASTEP: spin-polarised insulator reported as metallic (zero gap) when the top valence band exceeds E_F by ~0.1 meV

### Summary

For a spin-polarised CASTEP band-structure calculation of an insulating 162-atom supercell (rutile MgH₂ with one Fe substituting Mg), `sumo.io.castep.band_structure()` returns a `BandStructure` for which `is_metal()` is `True`, so `get_band_gap()` reports a gap of 0 eV. Counting electrons on the same eigenvalues (band N vs N+1, with N taken from the `.bands` header) gives the expected gap of 1.394 eV.

### What happens

1. CASTEP writes one Fermi energy per spin channel in the `.bands` header (`0.067501` and `0.067499` Ha).
2. The `BandStructure` is built with a single Fermi energy, apparently the spin-up value (1.836796 eV here), which is the higher of the two.
3. Along the band path, the highest occupied band (band 323 of 323 electrons per spin) lies **0.044 meV (spin up) and 0.167 meV (spin down) above** that Fermi energy. This is not explained by the printed precision (half a unit of the last digit is ≈ 0.014 meV); my guess is that, with fixed occupancies, CASTEP determines E_F on the self-consistent k-point mesh, so path eigenvalues can exceed it by a sub-meV amount. Either way, the offset is far below any physical energy scale.
4. pymatgen's `BandStructure.is_metal(efermi_tol=1e-4)` treats any band crossing E_F by more than 0.1 meV as metallic, so the spin-down channel triggers the metallic classification and the gap becomes 0.

With `efermi_tol=5e-4` the system is correctly classified as an insulator (gap 1.394 eV, indirect M→R).

### Steps to reproduce

```python
from sumo.io.castep import band_structure

bs = band_structure("mgh2_Fe.bands")
print(bs.efermi)                      # 1.836796
print(bs.is_metal())                  # True
print(bs.get_band_gap()["energy"])    # 0.0
print(bs.is_metal(efermi_tol=5e-4))   # False
```

### Possible fixes (happy to help with a PR)

- For CASTEP `.bands` files, determine the occupied bands from the electron count in the header (number of electrons per spin) when deciding whether the system is metallic and when computing the gap; or
- use a tolerance larger than the default 0.1 meV (e.g. 1 meV) when checking for metallicity of CASTEP band structures.

Note that using the spin-resolved Fermi energies would not help here: relative to the spin-down Fermi level, the spin-down valence-band maximum is 0.22 meV higher. Since the tolerance lives in pymatgen's `is_metal`, please let me know if you think this is better reported there.

### Environment

sumo 3.0.0, pymatgen (latest release), Python 3.12. I can share the `.bands` file (≈ 450 kB) if useful.

This came up while cross-validating another post-processing tool against sumo on the same files; apart from this case, the band gaps, band edges and Fermi energies agreed to the printed precision. Thanks for maintaining sumo!
