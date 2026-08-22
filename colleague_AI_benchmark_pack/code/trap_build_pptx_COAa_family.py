"""
COAa-family deck: regions with same numeric pattern as COAa
  Passive > Active at Post AND Withdrawal
  Active > Passive at During AND Reinstatement

Figs + PPTX like behavior_match_12regions deck.
Out: forMark/COAa_family/
     forMark/TRAP_COAa_family.pptx
"""
from pathlib import Path
from datetime import datetime
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyBboxPatch
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from PIL import Image
import trap_region_zoomin as Z

ROOT = Path(r"C:\Users\hsollim\behavior_task\TRAP_analysis_sync")
ROSTER = ROOT / "TRAP_OUTPUT_calculated_mm3" / "13_universal_cluster_PCA_density" / \
    "forebrain_no_bs" / "z_within_phase" / "02_cluster_region_roster.csv"
WB = ROOT / "Hansol Lim 561 cell counts + densitynew.xlsx"
MANIFEST = ROOT / "TRAP_sample_manifest.csv"
OUT = ROOT / "TRAP_OUTPUT_calculated_mm3" / "forMark" / "COAa_family"
OUT.mkdir(parents=True, exist_ok=True)
PPTX = OUT.parent / "TRAP_COAa_family.pptx"

PH = ["During", "Post", "Withdrawal", "Reinstatement"]
PHS = ["During", "Post", "Withdr.", "Reinst."]
FULL = {
    "COAa": "Cortical amygdalar area, anterior",
    "LGv": "Ventral lateral geniculate complex",
    "FC": "Fasciola cinerea",
    "ASO": "Accessory supraoptic group",
    "SCH": "Suprachiasmatic nucleus",
    "RH": "Rhomboid nucleus",
    "ARH": "Arcuate hypothalamic nucleus",
    "FF": "Fields of Forel",
    "PIL": "Posterior intralaminar thalamic nucleus",
    "SubG": "Subgeniculate nucleus",
    "ME": "Median eminence",
    "COAp": "Cortical amygdalar area, posterior",
}

roster = pd.read_csv(ROSTER)
reg2clu = dict(zip(roster["acronym"], roster["cluster"]))
wb = pd.read_excel(WB, sheet_name="All Samples")
wb["base"] = wb["acronym"].astype(str).str.replace(r"-[LR]$", "", regex=True)
man = pd.read_csv(MANIFEST)
man = man[man["include"] == 1]
meta = man.drop_duplicates("mouse_id").set_index("mouse_id")[["delivery", "phase"]].to_dict("index")
dens_cols = [c for c in wb.columns if c.endswith("density (cells/mm^3)")]

rows = []
for reg, clu in reg2clu.items():
    rsub = wb[wb["base"] == reg]
    if rsub.empty:
        continue
    rec = {"region": reg, "cluster": int(clu)}
    ok = True
    for grp in ("Active", "Passive"):
        for ph in PH:
            vals = []
            for dcol in dens_cols:
                m = re.match(r"(HaLi_\d+_\d+)", dcol)
                if not m or m.group(1) not in meta:
                    continue
                if meta[m.group(1)]["delivery"] != grp or meta[m.group(1)]["phase"] != ph:
                    continue
                v = pd.to_numeric(rsub[dcol], errors="coerce").mean()
                if pd.notna(v):
                    vals.append(v)
            rec[f"{grp[0]}_{ph}"] = np.mean(vals) if vals else np.nan
            if not vals:
                ok = False
    if not ok:
        continue
    for ph in PH:
        rec[f"delta_{ph}"] = rec[f"A_{ph}"] - rec[f"P_{ph}"]
    # COAa-like strict
    rec["coa_like"] = (
        rec["P_Post"] > rec["A_Post"]
        and rec["P_Withdrawal"] > rec["A_Withdrawal"]
        and rec["A_During"] > rec["P_During"]
        and rec["A_Reinstatement"] > rec["P_Reinstatement"]
    )
    # softer: only P>A at Post+WD
    rec["p_gt_a_post_wd"] = rec["P_Post"] > rec["A_Post"] and rec["P_Withdrawal"] > rec["A_Withdrawal"]
    rec["delta_Post"] = rec["P_Post"] - rec["A_Post"]
    rec["delta_WD"] = rec["P_Withdrawal"] - rec["A_Withdrawal"]
    rec["score"] = rec["delta_Post"] + rec["delta_WD"]
    apeak = PH[int(np.argmax([rec[f"A_{p}"] for p in PH]))]
    ppeak = PH[int(np.argmax([rec[f"P_{p}"] for p in PH]))]
    rec["Active_peak"] = apeak
    rec["Passive_peak"] = ppeak
    rows.append(rec)

df = pd.DataFrame(rows)
df = df.sort_values(["coa_like", "score"], ascending=[False, False])
df.to_csv(OUT / "COAa_family_screen.csv", index=False)

strict = df[df["coa_like"]].copy()
soft = df[df["p_gt_a_post_wd"]].copy()
print("Strict COAa-like (A>P During+Rein & P>A Post+WD):", strict["region"].tolist())
print("Soft P>A Post+WD:", soft["region"].tolist())

# Use strict as core family; if few, also include soft top for deck
REGS = strict["region"].tolist()
# Always ensure COAa first
if "COAa" in REGS:
    REGS = ["COAa"] + [r for r in REGS if r != "COAa"]
CLUSTER = {r: int(df[df.region == r].iloc[0].cluster) for r in REGS}

# Also soft peers for heatmap context (not all in zoom)
SOFT_EXTRA = [r for r in soft["region"].tolist() if r not in REGS][:8]

# ========== FIGURES ==========
Z.REGIONS = REGS
Z.REGION_CLUSTER = CLUSTER
Z.OUT = OUT
long = Z.load_long()
from trap_region_zoomin import _withinphase_panel

x = np.arange(4)

# Heatmap mean A-P for strict family
fig, ax = plt.subplots(figsize=(9, max(3.8, 0.5 * len(REGS) + 2)))
M = np.array([[df[df.region == r].iloc[0][f"delta_{p}"] for p in PH] for r in REGS])
vmax = np.nanpercentile(np.abs(M), 98)
im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax.set_xticks(range(4)); ax.set_xticklabels(PHS, fontsize=11)
ax.set_yticks(range(len(REGS)))
ax.set_yticklabels([f"{r} (c{CLUSTER[r]})" for r in REGS], fontsize=11)
for i, r in enumerate(REGS):
    lw = 2.8 if r == "COAa" else 1.4
    ec = "#b30000" if r == "COAa" else "#333"
    ax.add_patch(plt.Rectangle((-0.5, i - 0.5), 4, 1, fill=False, edgecolor=ec, lw=lw))
    if r == "COAa":
        ax.get_yticklabels()[i].set_color("#b30000")
        ax.get_yticklabels()[i].set_fontweight("bold")
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
cb.set_label("Active − Passive (mean density, cells/mm³)", fontsize=9)
ax.set_title(
    "COAa family — RAW mean Active−Passive\n"
    "Blue at Post/WD = Passive > Active · Red at During/Rein = Active > Passive",
    fontsize=12, fontweight="bold", color="#1F3A5F",
)
fig.text(0.5, -0.02,
         "Strict pattern: Passive>Active at Post & Withdrawal; Active>Passive at During & Reinstatement · COAa = best Finding B anatomy",
         ha="center", fontsize=9,
         bbox=dict(boxstyle="round,pad=0.4", fc="#fff8e1", ec="#e0c060"))
fig.savefig(OUT / "heatmap_COAa_family_meanAP.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)

# Why COAa schematic
fig, ax = plt.subplots(figsize=(12, 5.2))
ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 5.5)
ax.text(6, 5.1, "COAa family — same numeric pattern; why pick COAa?",
        ha="center", fontsize=15, fontweight="bold", color="#1F3A5F")
boxes = [
    (0.4, 2.5, 3.5, 2.2, "#ffebee", "#c62828", "NUMERIC FAMILY",
     "\n".join(REGS) + f"\n\n({len(REGS)} regions)"),
    (4.25, 2.5, 3.5, 2.2, "#fff3e0", "#e65100", "PATTERN",
     "Passive > Active\nat Post & Withdrawal\n\nActive > Passive\nat During & Reinstatement"),
    (8.1, 2.5, 3.5, 2.2, "#e8f5e9", "#2e7d32", "CHOOSE COAa",
     "Cortical amygdala\n→ cue / context + affect\n\nBest fit for Finding B\nPIT-like schematic\n(LGv=visual; FC=hippocampal)"),
]
for x0, y0, w, h, fc, ec, t, body in boxes:
    ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=fc, ec=ec, lw=2))
    ax.text(x0 + w/2, y0 + h - 0.35, t, ha="center", fontsize=11, fontweight="bold", color=ec)
    ax.text(x0 + w/2, y0 + 0.2, body, ha="center", va="bottom", fontsize=9, color="#333")
ax.text(6, 1.5, "Candidate for cue/context (Post) + withdrawal-state → generalized motivation",
        ha="center", fontsize=11, fontweight="bold", color="#e65100",
        bbox=dict(boxstyle="round", fc="#fff8e1", ec="#e65100"))
ax.text(6, 0.55, "Not: Withdrawal-only peak  ·  Not: only region with this pattern  ·  Not: causal proof of PIT",
        ha="center", fontsize=9, color="#666", style="italic")
fig.savefig(OUT / "why_COAa_family.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)

# Contrast trajectories
ncols = min(3, len(REGS))
nrows = int(np.ceil(len(REGS) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.6 * nrows))
axes = np.atleast_2d(axes)
for ax, reg in zip(axes.flatten(), REGS):
    d = long[long.region == reg]
    for grp, col, mk in [("Active", "#d62728", "o"), ("Passive", "#1f77b4", "s")]:
        g = d[d.group == grp]
        mu, se = [], []
        for i, ph in enumerate(PH):
            vals = g[g.phase == ph].density.dropna().values
            if len(vals):
                rng = np.random.RandomState(abs(hash(reg+grp+ph)) % (2**31))
                jit = (rng.rand(len(vals)) - 0.5) * 0.18
                off = -0.08 if grp == "Active" else 0.08
                ax.scatter(np.full(len(vals), i)+off+jit, vals, s=24, color=col, alpha=0.65,
                           edgecolors="k", linewidths=0.3, zorder=3)
                mu.append(np.mean(vals))
                se.append(np.std(vals, ddof=1)/np.sqrt(len(vals)) if len(vals) > 1 else 0)
            else:
                mu.append(np.nan); se.append(0)
        ax.errorbar(x, mu, yerr=se, color=col, lw=2.2, marker=mk, ms=6, capsize=3, zorder=4)
    ax.axvspan(0.5, 2.5, color="#fff3e0", alpha=0.5, zorder=0)
    ax.set_xticks(x); ax.set_xticklabels(PHS, fontsize=8)
    title_c = "#b30000" if reg == "COAa" else "#333"
    ax.set_title(f"{reg} (c{CLUSTER[reg]})", fontsize=11, fontweight="bold", color=title_c)
    ax.grid(axis="y", color="#eee", lw=0.4)
for j in range(len(REGS), nrows * ncols):
    axes.flatten()[j].axis("off")
axes.flatten()[0].legend(["Active", "Passive"], fontsize=8)
fig.suptitle("COAa family — density by phase (orange = Post+WD where Passive > Active)",
             fontsize=12, fontweight="bold", color="#1F3A5F")
plt.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(OUT / "summary_density_COAa_family.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)

# Per region density/count + withinphase
for reg in REGS:
    Z.region_figure(long, reg)
    d = long[long.region == reg]
    fig, (axd, axc) = plt.subplots(2, 1, figsize=(12, 9))
    _withinphase_panel(axd, d, "density", "density (cells/mm³)",
                       f"{reg} · density — one bar = one mouse")
    _withinphase_panel(axc, d, "count", "cell count (L+R total)",
                       f"{reg} · cell count — one bar = one mouse")
    note = " ★ Finding B pick" if reg == "COAa" else " (numeric peer)"
    fig.suptitle(
        f"{reg} (Cluster {CLUSTER[reg]}){note}\n"
        "COAa-like: Passive>Active at Post & WD · Active>Passive at During & Rein",
        fontsize=12, fontweight="bold",
    )
    plt.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / f"{reg}_withinphase_ActiveVsPassive.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("plotted", reg)

# withinphase summary
fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.8 * nrows))
axes = np.atleast_2d(axes)
for ax, reg in zip(axes.flatten(), REGS):
    d = long[long.region == reg]
    _withinphase_panel(ax, d, "density", "density", reg)
    ax.set_title(f"{reg} (c{CLUSTER[reg]})", fontsize=10, fontweight="bold",
                 color="#b30000" if reg == "COAa" else "#333")
    if ax.get_legend():
        ax.get_legend().remove()
for j in range(len(REGS), nrows * ncols):
    axes.flatten()[j].axis("off")
axes.flatten()[0].legend(handles=[Patch(color="#d62728", label="Active"),
                                   Patch(color="#1f77b4", label="Passive")], fontsize=8)
fig.suptitle("COAa family — within-phase (one bar = one mouse)",
             fontsize=12, fontweight="bold", color="#1F3A5F")
plt.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(OUT / "summary_withinphase_COAa_family.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("Figs done:", OUT)

# ========== PPTX ==========
NAVY = RGBColor(0x1F, 0x3A, 0x5F)
GREY = RGBColor(0x55, 0x55, 0x55)
ORANGE = RGBColor(0xE6, 0x51, 0x00)
RED = RGBColor(0xB3, 0x00, 0x00)
BLUE = RGBColor(0x15, 0x65, 0xC0)
GREEN = RGBColor(0x2E, 0x7D, 0x32)

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]
SW, SH = 13.333, 7.5


def title(s, text, color=NAVY, size=24):
    tb = s.shapes.add_textbox(Inches(0.4), Inches(0.2), Inches(SW - 0.8), Inches(0.8))
    p = tb.text_frame.paragraphs[0]
    p.text = text
    p.font.size = Pt(size)
    p.font.bold = True
    p.font.color.rgb = color


def bullets(s, items, left, top, w, h, size=15):
    tb = s.shapes.add_textbox(Inches(left), Inches(top), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    for i, (txt, c) in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = "•  " + txt
        p.font.size = Pt(size)
        p.font.color.rgb = c
        p.space_after = Pt(7)


def img_fit(s, path, left, top, maxw, maxh):
    if not Path(path).exists():
        return
    iw, ih = Image.open(path).size
    r = min(maxw / iw, maxh / ih)
    w, h = iw * r, ih * r
    s.shapes.add_picture(str(path), Inches(left + (maxw - w) / 2), Inches(top + (maxh - h) / 2),
                         width=Inches(w), height=Inches(h))


def slide_big(t, img, caption, tcolor=ORANGE):
    s = prs.slides.add_slide(BLANK)
    title(s, t, tcolor)
    img_fit(s, img, 0.35, 1.1, SW - 0.7, SH - 1.95)
    tb = s.shapes.add_textbox(Inches(0.4), Inches(SH - 0.7), Inches(SW - 0.8), Inches(0.55))
    p = tb.text_frame.paragraphs[0]
    p.text = caption
    p.font.size = Pt(11)
    p.font.color.rgb = GREY


def slide_img_text(t, img, items, tcolor=ORANGE, img_w=9.0):
    s = prs.slides.add_slide(BLANK)
    title(s, t, tcolor)
    img_fit(s, img, 0.3, 1.1, img_w, SH - 1.5)
    bullets(s, items, img_w + 0.5, 1.4, SW - img_w - 0.9, SH - 2, size=13)


# title
s = prs.slides.add_slide(BLANK)
tb = s.shapes.add_textbox(Inches(0.8), Inches(2.2), Inches(SW - 1.6), Inches(3))
tf = tb.text_frame
tf.word_wrap = True
p = tf.paragraphs[0]
p.text = "COAa family — Finding B candidates"
p.font.size = Pt(32)
p.font.bold = True
p.font.color.rgb = NAVY
p2 = tf.add_paragraph()
p2.text = "Same numeric pattern as COAa · Passive > Active at Post & Withdrawal"
p2.font.size = Pt(16)
p2.font.color.rgb = ORANGE
p3 = tf.add_paragraph()
p3.text = f"{len(REGS)} regions · COAa = best anatomical fit for cue/context + withdrawal-state"
p3.font.size = Pt(14)
p3.font.color.rgb = GREY

# criteria
s = prs.slides.add_slide(BLANK)
title(s, "Selection rule — COAa-like pattern", ORANGE)
bullets(s, [
    ("PASSIVE > ACTIVE at Post AND Withdrawal.", BLUE),
    ("ACTIVE > PASSIVE at During AND Reinstatement.", RED),
    (f"Strict matches: {', '.join(REGS)}.", ORANGE),
    ("COAa chosen for Finding B schematic (cortical amygdala = cue/affect).", GREEN),
    ("LGv (visual thalamus) & FC (hippocampal) = numeric peers, weaker anatomy.", GREY),
    ("Separate from Finding A (ORBm/BMAp) and from behavior-match WD-peak set (SI/EPd/AVP).", GREY),
], 0.7, 1.4, SW - 1.4, 5.5, size=16)

# roster
s = prs.slides.add_slide(BLANK)
title(s, "COAa family roster", ORANGE)
lines = []
for r in REGS:
    name = FULL.get(r, r)
    tag = " ★ Finding B pick" if r == "COAa" else " (numeric peer)"
    lines.append((f"{r} (c{CLUSTER[r]}) — {name}{tag}", ORANGE if r == "COAa" else GREY))
bullets(s, lines, 0.7, 1.3, SW - 1.4, 5.5, size=15)

slide_img_text(
    "Why COAa among this family?",
    OUT / "why_COAa_family.png",
    [
        ("Same numeric pattern in all listed regions.", GREY),
        ("COAa = cortical amygdala → cue + affect.", ORANGE),
        ("Best fit for Finding B PIT schematic.", ORANGE),
        ("Candidate — not causal proof.", GREY),
    ],
)

slide_big(
    "Heatmap — RAW mean Active − Passive (cells/mm³)",
    OUT / "heatmap_COAa_family_meanAP.png",
    "NOT z-score · Blue at Post/WD = Passive higher · Red box = COAa.",
)

slide_big(
    "All family — density by phase",
    OUT / "summary_density_COAa_family.png",
    "Orange band = Post + Withdrawal (Passive > Active).",
)

slide_big(
    "All family — within-phase bars (one bar = one mouse)",
    OUT / "summary_withinphase_COAa_family.png",
    "Red = Active · Blue = Passive · horizontal line = mean.",
)

for reg in REGS:
    tcolor = ORANGE if reg == "COAa" else GREY
    tag = " ★ Finding B" if reg == "COAa" else " (peer)"
    slide_big(
        f"{reg} (c{CLUSTER[reg]}){tag} — density & cell count",
        OUT / f"{reg}_density_count_by_phase.png",
        f"{FULL.get(reg, reg)} · Passive>Active at Post & WD · Active>Passive at During & Rein.",
        tcolor,
    )
    slide_big(
        f"{reg} (c{CLUSTER[reg]}){tag} — within-phase (one bar = one mouse)",
        OUT / f"{reg}_withinphase_ActiveVsPassive.png",
        "One bar = one mouse · horizontal line = group mean.",
        tcolor,
    )

# summary
s = prs.slides.add_slide(BLANK)
title(s, "Summary — COAa family / Finding B")
bullets(s, [
    (f"Numeric family ({len(REGS)}): " + ", ".join(REGS) + ".", GREY),
    ("Pattern: Passive > Active at Post & Withdrawal; Active > Passive at During & Reinstatement.", ORANGE),
    ("COAa = best anatomical fit for cue/context + withdrawal-state → generalized motivation candidate.", ORANGE),
    ("Complementary to Finding A (ORBm/BMAp Active craving) and WD-peak behavior-match set.", NAVY),
], 0.8, 1.8, SW - 1.6, 4.5, size=16)

out = None
for cand in [PPTX] + [PPTX.with_name(f"{PPTX.stem}_alt{i}.pptx") for i in range(1, 4)]:
    try:
        prs.save(str(cand))
        out = cand
        break
    except PermissionError:
        continue
if out is None:
    out = PPTX.with_name(f"{PPTX.stem}_{datetime.now():%H%M%S}.pptx")
    prs.save(str(out))
print(f"Saved: {out}\nSlides: {len(prs.slides)}")
