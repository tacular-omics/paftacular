# Copilot instructions for paftacular

Read [`CLAUDE.md`](../CLAUDE.md) at the repo root. It is the canonical guide: commands,
architecture, public API, scientific conventions, parser rules and gotchas. Library usage is
in [`llms-full.txt`](../llms-full.txt). This repo implements mzPAF 1.0.1. ProForma
interpretation comes from the optional peptacular dependency.

Key rules:

1. `mass()` is the charged-species mass. Without an embedded or resolved sequence, peptide,
   internal and precursor ions return only offsets and modifiers. Keep that behavior and use
   `resolve()` for complete masses.
2. Never split annotations on bare commas. Use `syntax.annotation_spans`. Keep the possessive
   `*+` quantifiers in `constants._ATOM_TOKEN` (ReDoS regression test).
3. The core must import and parse without any extra. peptacular, pysmiles and the MCP SDK are
   imported lazily and fail with an install hint. Keep SDK imports inside `mcp/`.
4. Before a commit: `just check` (ruff check and format on src tests benchmarks, `ty check src`,
   pytest). Never use em dashes or semicolons in authored text or code comments.
5. Never bump the version, tag, or publish. Only the tacular-omics overseer releases.
