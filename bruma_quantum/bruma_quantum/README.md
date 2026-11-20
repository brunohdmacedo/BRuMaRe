# bruma_quantum

Entropia, teoria da informação e descritores de emaranhamento a partir de cálculos
DFT do **CASTEP**. Só `numpy`. Pensado como módulo para o BruMaRe (cada função é pura e
devolve `dict`/`ndarray`; o I/O fica em `io_castep.py` e `weights.py`).

```
python -m bruma_quantum <comando> ...        # testes: python test_bruma_quantum.py
```

## Comandos

| Comando | Entrada | Saída |
|---|---|---|
| `density` | `.den_fmt` (+ `--cell --center --radius`, `--ref`) | Shannon, Rényi-2, Onicescu, Fisher por canal (total/up/down), JS(up,down), região esférica, KL/JS contra o puro |
| `shells` | `.cell`, `--center Ni`, `--target H` | cascas de H por distância ao dopante + multiplicidade g; gera CSV-modelo |
| `hsites` | CSV de cascas com `E_eV` | P_i de Boltzmann, S configuracional (famílias e sítios), levantamento da degenerescência, vs T |
| `orbital` | `.bands` + pesos de projeção `.dat` (+ `--region`) | S de emaranhamento (cota diagonal), entropia por banda, localização (S_loc, N_eff, IPR), espaço ativo/qubits |
| `fermi` | `.bands` | entropia térmica de ocupação (Fermi-Dirac) vs T |

Fluxo para o H no MgH₂:Ni: `shells` → calcular a energia de remoção de H de cada
casca (CASTEP/Bruma) → preencher `E_eV` → `hsites`.

## O que cada grandeza significa

- **Shannon da densidade** `S = -∫ρ̃ ln ρ̃`, `V_ef = exp(S)` (volume efetivo em Å³). Só
  faz sentido comparando a mesma célula/grade (defeito × puro).
- **Fisher** `I = 4∫|∇√ρ̃|²` (FFT, rede qualquer): mede o quão "pontuda" é a densidade.
- **Onicescu / Rényi-2**: desequilíbrio (concentração) da distribuição.
- **KL / Jensen-Shannon** (mesma grade): onde e quanto o defeito mudou a densidade.
- **Emaranhamento (von Neumann) KS**: `S = Σ h(ν)`, `h(x) = -x ln x - (1-x) ln(1-x)`
  (Peschel). Função `peschel_entropy(C)` é **exata** quando se tem a matriz de overlap
  `C_ij = <ψ_i|P_A|ψ_j>` dos orbitais ocupados. Com só pesos de projeção `w_nk` (diagonal
  de C) obtém-se uma **cota superior** (`diagonal_entropy`).
- **Localização**: `S_loc = -Σ_a p_a ln p_a`, `N_eff = exp(S_loc)`, IPR = `Σ p_a²`.
- **Espaço ativo**: bandas com |⟨E⟩-E_F| ≤ janela e peso na região ≥ `--wmin`;
  `spin-orbitais ≈ qubits` (Jordan-Wigner, 1 qubit por spin-orbital).
- **Boltzmann de sítios**: `P_i = g_i e^{-E_i/kT}/Z`; `S_sítios = -Σ P_i ln(P_i/g_i)`.

## Limitações (leia antes de interpretar)

1. **É emaranhamento do sistema auxiliar de Kohn-Sham**, não do sistema interagente.
   Use para tendências relativas (Ni × Cr × Zn × puro), não como valor absoluto.
2. **A cota diagonal é frouxa para estados estendidos** (num teste de cadeia
   tight-binding: exata 0,76 × diagonal 6,93). O número total é dominado por bandas de
   bulk; olhe a entropia **por banda** e a **diferença para o puro** com a mesma região.
3. Pesos de projeção LCAO têm *spilling* (não somam 1); a cota só vale com projeção
   completa. A mistura entre k-pontos é ignorada (média ponderada em k).
4. `.pdos_bin` **não** é lido aqui (binário Fortran não documentado). Use os pesos
   exportados pelo Bruma; se o layout do `.dat` for diferente do descrito em
   `weights.py`, ajuste só aquele módulo.
5. `fermi` usa autovalores congelados e μ = E_F do `.bands`. Se houver banda de defeito
   colada em E_F, o resultado depende fortemente de μ; confira `N(T)` (deve ≈ N_elétrons).
6. `shells` agrupa por distância (triagem): confirme equivalência por simetria.
7. Valores de `.den_fmt` são ρ·V; o pacote divide pelo volume. Ruído numérico negativo
   (~1e-6 da carga) é truncado em zero. Densidade de **diferença** não serve para Shannon.

## Integração no BruMaRe

Funções públicas (`from bruma_quantum import ...`): `read_bands`, `read_cell`,
`read_den_fmt`, `analyze_density`, `shannon_report`, `fisher_information`,
`kl_divergence`, `jensen_shannon`, `sphere_mask`, `read_weights`, `entanglement_report`,
`peschel_entropy`, `diagonal_entropy`, `localization_entropy`, `suggest_active_space`,
`thermal_entropy`, `find_shells`, `boltzmann_site_entropy`.
Saídas `.dat` usam TAB + linha de nomes (Origin) e os relatórios saem em `.json`.

## Próximos passos naturais

- Peschel **exato**: exige coeficientes de onda plana (ou funções de Wannier) para montar C.
- Informação mútua entre regiões A e B (precisa da mesma matriz completa).
- Entropia de orbital de uma partição de espaço ativo para DMRG/VQE.
