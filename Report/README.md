# Report/ — results reports on production runs

A **report** here is a short PDF about one set of production runs, for the advisor and for
the project record. It is not the whitepaper (`Whitepaper/`), which describes the method.
`populations/populations_report.tex` is the finished example: when in doubt, do what it does.

## Layout

```
Report/
  README.md             this file
  TEMPLATE_report.tex   the skeleton every report starts from
  refs.bib              shared bibliography (reports cite it as ../refs)
  <topic>/              ONE directory per report: its .tex, its PDF, its build products
    <topic>_report.tex
```

Figures are never copied into a report directory; `\graphicspath` points at
`../../figures/<run>/`, where the analysis scripts wrote them.

## Starting a new report

```bash
mkdir Report/<topic>
cp Report/TEMPLATE_report.tex Report/<topic>/<topic>_report.tex
# point \graphicspath at the figure directories the analysis scripts wrote
# fill every \TODO, delete sections that do not apply
cd Report/<topic> && latexmk -pdf <topic>_report.tex
grep -c '^TODO:' <topic>_report.log        # must print 0 before the report is shared
```

`\TODO{...}` prints red in the PDF and logs itself, so an unfinished report cannot pass for a
finished one.

## Where the numbers come from

Every number in a report comes from a script's output file, **never typed in from a
terminal or from memory**. Cite the output file in the text or the reproduction block.

| Content | Script (`analysis/`) | Output |
|---|---|---|
| detection share, joint-fit gain, astrometry, resolution | `p6_synergy_resolution.py` | `figures/<run>/p6_*` |
| precision CDFs, mass–distance, sky maps | `p7_forecast_figures.py` | `figures/<run>/p7_*` |
| absolute yields at an F grid, plus validation | `y1_absolute_yield.py` | `y1_yields.{md,csv}` |
| yield per unit F, N(F) figure | `y3_yield_vs_F.py` | `y3/` |
| comparison with Sajadian & Sahu (2023) | `y2_ss23_compare.py` | `y2_ss23.md` |
| sample light curves | `s2_sample_lightcurves.py` | `figures/samples/` |

Every script that reads an event table goes through `analysis/romanlib.py`, which applies
the event-rate weight and the `-1.0` sentinel rules (CLAUDE.md, "Python side"). The p6, p7 and
y1 scripts take whole run directories (`--run name=dir`); `y3` reads only y1's CSV.

## Checklist before a report is shared

1. **No `TODO` left:** `grep -c '^TODO:' *.log` prints 0.
2. **It builds clean:** no `undefined` references or citations in the `.log`, and no BibTeX
   warnings in the `.blg`.
3. **Every figure was read** at page size after the build: no overlapping labels, legend
   clear of the data, floats next to their text.
4. **Every pooled number is weighted and has its N_eff**, and every denominator is named.
5. **Counts are labelled** as a Monte Carlo sample size or as a yield. A yield carries its F
   and its counting convention (per resolved object or all stars).
6. **Known-bug status is stated:** if the runs predate a fix (see `DEVIATIONS.md`), the
   abstract and the caveats say so and give the direction of the bias.
7. **Literature comparisons are like for like:** same mass function, abundance, area, seasons,
   detection cut and denominator (see DEVIATIONS 54a for how this went wrong once).
8. **The bibliography is checked:** each new `refs.bib` entry verified against arXiv or the
   publisher, and only papers actually cited. `refs.bib` is shared by all reports.
9. **The record is updated:** a `PROGRESS.md` line naming the report and its data, and a
   `DEVIATIONS.md` entry if writing it turned up anything new.

## Files

- `TEMPLATE_report.tex`: the skeleton; builds as is (from a `Report/<topic>/` copy), full
  of red TODOs.
- `populations/populations_report.tex`: the black-hole / neutron-star report (runs of 2026-09-17/18,
  **pre-extinction-fix**). The post-fix re-runs get a new report, not an edit of this one
  (the user's decision, 2026-09-22).
- `refs.bib`: shared bibliography, natbib/`plainnat`.
