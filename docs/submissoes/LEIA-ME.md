# BrUMaRe: manuscritos prontos para seis revistas

> **Submeta a UMA revista por vez.** Enviar o mesmo trabalho a duas revistas ao mesmo tempo é má conduta.
> Só envie a outra se a primeira recusar ou se você retirar a submissão.

| Pasta | Revista | Tipo de artigo | Acesso aberto | Formato |
|---|---|---|---|---|
| `scipost-codebases/` | SciPost Physics Codebases | artigo completo + versão do código | diamante (grátis) | LaTeX; modelo oficial SciPost no Overleaf |
| `joss/` | Journal of Open Source Software | artigo curto (918 palavras) | diamante (grátis) | `paper.md` + `paper.bib` (o JOSS gera o PDF) |
| `computational-materials-science/` | Computational Materials Science (Elsevier) | artigo completo | híbrida: grátis pelo acordo CAPES–Elsevier* | LaTeX `elsarticle` (Overleaf) |
| `j-comput-chem/` | Journal of Computational Chemistry (Wiley) | Software News and Updates | híbrida: grátis pelo acordo CAPES–Wiley* | LaTeX livre (1ª submissão) |
| `jcim/` | J. Chem. Inf. Model. (ACS) | Application Note | híbrida: grátis pelo acordo CAPES–ACS* | LaTeX `achemso` (modelo oficial da ACS) |
| `braz-j-phys/` | Brazilian Journal of Physics (Springer) | artigo completo | híbrida: grátis pelo acordo CAPES–Springer Nature* | LaTeX livre (1ª submissão) |

\* Confirme com a biblioteca da UNILA que a instituição e a revista estão no acordo e que **você** é o autor correspondente.
A versão para a **Computer Physics Communications** continua em `../overleaf-brumare/` (acordo CAPES–Elsevier).

## Antes de submeter (qualquer revista)
1. **Revise a declaração de uso de IA** (presente em todos os manuscritos): ela precisa descrever exatamente o que foi feito. IA não pode ser autora; a responsabilidade é sua.
2. **Coautores:** decida com os orientadores; nas versões LaTeX eles estão indicados em comentários.
3. **Repositório público** no GitHub com a versão 1.2.0 **arquivada no Zenodo (DOI)**. O DOI entra na seção de disponibilidade.
4. **Testes automatizados** (por exemplo `pytest` com os arquivos de validação) e **guia de contribuição**: o JOSS exige; as outras revistas valorizam.
5. **Comparação com outras ferramentas** (sumo, pymatgen, Bader) nos mesmos dados: é o pedido mais provável dos revisores.
6. **Confira os números da aplicação (MgH₂:Fe)** com os seus cálculos finais.

## Particularidades
- **JOSS:** as seções obrigatórias já estão no `paper.md` (Summary, Statement of need, State of the field, Software design, Research impact statement, AI usage disclosure). A revisão acontece no GitHub, sobre o código.
- **SciPost Codebases:** além do artigo, publica uma "codebase release" (versão do código). Prepare o arquivamento da versão.
- **JCIM:** exige o resumo gráfico (TOC graphic, já incluído) e a seção "Data and Software Availability" (incluída).
- **Elsevier (CMS e CPC):** destaques (highlights) incluídos e a seção "Declaration of generative AI and AI-assisted technologies in the writing process".
- **Wiley e Springer:** aceitam formato livre na primeira submissão; após o aceite, use o modelo da editora.
