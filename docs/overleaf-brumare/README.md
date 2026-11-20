# BrUMaRe — artigo para a Computer Physics Communications (projeto Overleaf)

## Como abrir no Overleaf
1. Overleaf → **New Project → Upload Project** → envie `overleaf-brumare.zip`.
2. **Menu → Compiler: pdfLaTeX** e **Main document: main.tex**.
3. Clique em **Recompile**. O Overleaf já tem a classe oficial `elsarticle` e o estilo `elsarticle-num`.

## Onde está cada parte
| Pasta | Conteúdo |
|---|---|
| `main.tex` | arquivo principal: monta tudo na ordem da revista |
| `frontmatter/` | `title.tex`, `authors.tex`, `abstract.tex`, `highlights.tex`, `keywords.tex`, `program_summary.tex` |
| `sections/` | `01_introduction.tex` … `09_conclusions.tex` (uma seção por arquivo) |
| `backmatter/` | agradecimentos, disponibilidade de código e dados, CRediT, conflito de interesses |
| `appendix/` | `implementation_map.tex` (onde cada método está no código) |
| `tables/` | `readers.tex` (leitores e unidades), `validation.tex` (validação) |
| `figures/` | imagens `.png` e `architecture.tex` (diagrama em TikZ) |
| `bib/` | `references.bib` |
| `macros/` | `preamble.tex` (pacotes) e o comando `\BMR` (nome BrUMaRe com "Re" em âmbar) |

## Antes de submeter
- **Coautores:** os orientadores estão comentados em `frontmatter/authors.tex`; descomente após a concordância deles.
- **Highlights:** a Elsevier pede 3 a 5 itens com até 85 caracteres (já estão em `frontmatter/highlights.tex`).
- **Figuras:** a revista aceita PNG de alta resolução; para a versão final, gere as figuras em 300 dpi ou em PDF vetorial.
- **Números da aplicação (MgH₂:Fe):** confira com seus cálculos finais antes de submeter.
