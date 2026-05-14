# SCHEMA — how to maintain this wiki (humans + LLM)

## Layers

| Layer | What it is | Rule |
|-------|------------|------|
| **Raw** | Code, CSV outputs, configs in each GitHub repo | LLMs **read**; do not replace with prose here as if it were code. |
| **Wiki (this folder)** | Narrative, navigation, assumptions, cohort story | LLMs **edit** with human review; keep pages short; link out. |

## Ingest (when something new lands)

0. **Proactive rule:** substantive code/config/analysis changes made in chat should also touch this hub (at least `log.md` + one sentence on the right page) even if the user did not say “wiki” — unless they said **no wiki** / **위키 빼고**, or the change is trivial (typos/format/comments only).

1. Identify **which code repo** owns the fact (PR, commit, path).
2. Add or update **one** canonical paragraph here + link to GitHub file or README section.
3. If multiple repos are affected, update [cross/opioid_TRAP_genelist.md](cross/opioid_TRAP_genelist.md) or the relevant per-repo page — not three copies of the same paragraph.
4. Append a dated line to [log.md](log.md).

## Links

- Prefer stable URLs: **tag**, **release**, or **commit** for “this is exactly what we ran.” `main` is OK for “current entry point.”
- Cross-links between wiki pages use **relative** paths (e.g. `cross/opioid_TRAP_genelist.md`).

## Lint (occasional)

- Broken relative links, orphan pages not reachable from [index.md](index.md).
- Contradictions between two wiki pages → pick one canonical page, fix the other to point to it.

## Safety

- No secrets, tokens, or private identifiers in public wiki text.
- Local drive letters and internal-only paths: mark clearly as **local**; consider separate private notes if needed.

## Reference pattern

- Karpathy LLM Wiki gist: https://gist.github.com/karpathy/442a6bf555914893e9891c11519de94f
