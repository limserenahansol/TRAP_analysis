# Research wiki — index

**On GitHub:** [TRAP_analysis/wiki/research_hub](https://github.com/limserenahansol/TRAP_analysis/tree/main/wiki/research_hub) (this folder).

**Purpose:** One map for cross-repo context (cohort rules, paths, how pipelines relate). **Do not** treat this as executable truth for numbers — verify in code, manifests, and raw outputs.

## Code repositories (source of truth for behavior)

| Project | GitHub | Summary in this hub |
|---------|--------|---------------------|
| Opioid / behavior MATLAB | [opioidaddiction-matlab](https://github.com/limserenahansol/opioidaddiction-matlab) | [opioidaddiction-matlab.md](opioidaddiction-matlab.md) |
| TRAP | [TRAP_analysis](https://github.com/limserenahansol/TRAP_analysis) | [TRAP_analysis.md](TRAP_analysis.md) |
| Genelist / WMB | [Genelist_analysis_WMB](https://github.com/limserenahansol/Genelist_analysis_WMB) | [Genelist_analysis_WMB.md](Genelist_analysis_WMB.md) |

## Cross-project

- [Opioid ↔ TRAP ↔ Genelist — how they connect](cross/opioid_TRAP_genelist.md)

## Housekeeping

- [Schema / LLM ingest rules](SCHEMA.md)
- [Change log](log.md)

---

### Quick pointers for LLM sessions

1. Prefer **links to files/commits** in code repos over pasting large code blocks here.
2. When updating this wiki after a code change, add a line to [log.md](log.md).
3. Sensitive paths (e.g. local `K:\`) belong in **schema-approved** notes only; redact for public forks.
