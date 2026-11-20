# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025-2026 Bruno Henrique Dourado Macedo
# BrUMaRe: Brillouin-zone Utilities for Materials Research - https://github.com/brunohdmacedo/brumare
"""Interface do BrUMaRe para computador: a mesma do site (brumare_app/web), offline."""
import os
import sys


def web_dir():
    """Pasta com a interface (index.html, tutorial, bibliotecas e fontes embutidas)."""
    base = getattr(sys, "_MEIPASS", None)
    cands = ([os.path.join(base, "brumare_app", "web")] if base else []) + [os.path.join(os.path.dirname(os.path.abspath(__file__)), "web")]
    for c in cands:
        if os.path.isfile(os.path.join(c, "index.html")):
            return c
    raise FileNotFoundError("interface do BrUMaRe nao encontrada (brumare_app/web). Rode: python3 tools/montar_app.py")
