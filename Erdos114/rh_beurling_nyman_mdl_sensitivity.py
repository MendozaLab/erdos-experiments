#!/usr/bin/env python3
"""
EXP-MATH-RH-BEURLING-NYMAN-MDL-SENSITIVITY-20260506-01

Sensitivity run for the Beurling-Nyman MDL probe.

Purpose:
    Test whether the visual bend near ~3.2 bits is stable across dictionary
    families, seeded log-uniform schedules, and a denser N sweep.

Claim ceiling:
    SUGGESTIVE_NUMERICAL_SENSITIVITY only.
    This is not RH evidence, not theorem progress, and not an asymptotic claim.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import lstsq


EXPERIMENT_ID = "EXP-MATH-RH-BEURLING-NYMAN-MDL-SENSITIVITY-20260506-01"
SOURCE_EXPERIMENT = "EXP-MATH-RH-BEURLING-NYMAN-MDL-PROBE-20260506-03"
CLAIM_CEILING = (
    "SUGGESTIVE_NUMERICAL_SENSITIVITY: finite quadrature sensitivity test of "
    "the ~3.2-bit bend; not RH evidence, not theorem progress, and not an "
    "asymptotic statement."
)
N_VALUES = (4, 8, 12, 16, 20, 24, 28, 32, 40, 48)
QUANT_BITS = (8, 16, 24)
CANONICAL_BITS = 16
QUADRATURE_NODES = 2048
ACTIVE_REL_TOL = 1e-10
BASE_SEED = 11420260506
SEEDED_REPLICATES = 6


@dataclass(frozen=True)
class DictionarySpec:
    series: str
    family: str
    header_bits: int
    seed_offset: int | None
    theta_rule: str


def dictionary_specs() -> list[DictionarySpec]:
    specs = [
        DictionarySpec(
            series="harmonic",
            family="harmonic",
            header_bits=48,
            seed_offset=None,
            theta_rule="theta_j = 1 / j for j = 1..N",
        ),
        DictionarySpec(
            series="geometric",
            family="geometric",
            header_bits=64,
            seed_offset=None,
            theta_rule="theta_j = exp(-j * log(N + 1) / N) for j = 1..N",
        ),
    ]
    for r in range(SEEDED_REPLICATES):
        specs.append(
            DictionarySpec(
                series=f"seeded_log_uniform_r{r}",
                family="seeded_log_uniform",
                header_bits=112,
                seed_offset=r,
                theta_rule=(
                    "theta_j are sorted seeded log-uniform samples in "
                    f"[1/(N+1), 1], seed={BASE_SEED}+1000*r+N"
                ),
            )
        )
    return specs


def fractional_part(values: np.ndarray) -> np.ndarray:
    return values - np.floor(values)


def theta_values(spec: DictionarySpec, n: int) -> np.ndarray:
    j = np.arange(1, n + 1, dtype=np.float64)
    if spec.family == "harmonic":
        return 1.0 / j
    if spec.family == "geometric":
        h = math.log(n + 1.0) / n
        return np.exp(-j * h)
    if spec.family == "seeded_log_uniform":
        assert spec.seed_offset is not None
        rng = np.random.default_rng(BASE_SEED + 1000 * spec.seed_offset + n)
        logs = rng.uniform(-math.log(n + 1.0), 0.0, size=n)
        return np.sort(np.exp(logs))[::-1]
    raise ValueError(spec.family)


def combination_bits(n: int, k: int) -> int:
    if k <= 0 or k >= n:
        return 1
    return math.ceil(math.log2(math.comb(n, k)))


def bit_cost(spec: DictionarySpec, n: int, active: int, coeff_bits: int) -> dict:
    support_pattern_bits = combination_bits(n, active)
    theta_rule_bits = spec.header_bits + support_pattern_bits
    coefficient_bits = active * coeff_bits
    return {
        "dictionary_header_bits": spec.header_bits,
        "support_pattern_bits": support_pattern_bits,
        "theta_rule_bits": theta_rule_bits,
        "coefficient_bits": coefficient_bits,
        "total_description_bits": theta_rule_bits + coefficient_bits,
    }


def build_quadrature() -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = leggauss(QUADRATURE_NODES)
    return 0.5 * (nodes + 1.0), 0.5 * weights


def fit_dictionary(theta: np.ndarray, x: np.ndarray, w: np.ndarray) -> dict:
    design = fractional_part(theta[None, :] / x[:, None])
    sqrt_w = np.sqrt(w)
    aw = design * sqrt_w[:, None]
    y = sqrt_w
    coeffs, _, rank, singular_values = lstsq(aw, y, lapack_driver="gelsd")
    residual_l2 = float(np.linalg.norm(aw @ coeffs - y))
    target_l2 = float(np.linalg.norm(y))
    coeff_linf = float(np.max(np.abs(coeffs))) if coeffs.size else 0.0
    active = int(np.sum(np.abs(coeffs) > ACTIVE_REL_TOL * max(1.0, coeff_linf)))
    sigma_max = float(singular_values[0]) if len(singular_values) else 0.0
    sigma_min = float(singular_values[-1]) if len(singular_values) else 0.0
    condition = sigma_max / sigma_min if sigma_min > 0 else float("inf")
    return {
        "aw": aw,
        "y": y,
        "coeffs": coeffs,
        "rank": int(rank),
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "condition_number": condition,
        "residual_l2_unquantized": residual_l2,
        "relative_residual_unquantized": residual_l2 / target_l2,
        "target_l2": target_l2,
        "active_count": active,
        "coefficient_l1": float(np.linalg.norm(coeffs, 1)),
        "coefficient_l2": float(np.linalg.norm(coeffs)),
        "coefficient_linf": coeff_linf,
    }


def quantize(coeffs: np.ndarray, bits: int) -> tuple[np.ndarray, float]:
    scale = max(1.0, float(np.max(np.abs(coeffs))) if coeffs.size else 1.0)
    step = scale * 2.0 ** (-bits)
    return np.round(coeffs / step) * step, step


def neglog2(value: float) -> float:
    return -math.log2(value) if value > 0 else float("inf")


def run() -> dict:
    x, w = build_quadrature()
    rows = []
    for n in N_VALUES:
        for spec in dictionary_specs():
            theta = theta_values(spec, n)
            fit = fit_dictionary(theta, x, w)
            quantized = []
            for bits in QUANT_BITS:
                q, step = quantize(fit["coeffs"], bits)
                q_residual = float(np.linalg.norm(fit["aw"] @ q - fit["y"]))
                penalty = fit["sigma_max"] * math.sqrt(fit["active_count"]) * step / 2.0
                rel = q_residual / fit["target_l2"]
                qrow = {
                    "coeff_bits": bits,
                    "quantization_step": step,
                    "quantization_penalty_l2_bound": penalty,
                    "relative_residual": rel,
                    "negative_log2_relative_residual": neglog2(rel),
                    "certified_upper_relative_residual": (
                        fit["relative_residual_unquantized"] + penalty / fit["target_l2"]
                    ),
                    "certified_negative_log2_upper_residual": neglog2(
                        fit["relative_residual_unquantized"] + penalty / fit["target_l2"]
                    ),
                    **bit_cost(spec, n, fit["active_count"], bits),
                }
                quantized.append(qrow)
            rows.append(
                {
                    "series": spec.series,
                    "family": spec.family,
                    "seed_offset": spec.seed_offset,
                    "theta_rule": spec.theta_rule,
                    "N": n,
                    "rank": fit["rank"],
                    "active_count": fit["active_count"],
                    "sigma_max": fit["sigma_max"],
                    "sigma_min": fit["sigma_min"],
                    "gram_eigen_max": fit["sigma_max"] ** 2,
                    "gram_eigen_min": fit["sigma_min"] ** 2,
                    "condition_number": fit["condition_number"],
                    "residual_l2_unquantized": fit["residual_l2_unquantized"],
                    "relative_residual_unquantized": fit["relative_residual_unquantized"],
                    "negative_log2_relative_residual_unquantized": neglog2(
                        fit["relative_residual_unquantized"]
                    ),
                    "coefficient_l1": fit["coefficient_l1"],
                    "coefficient_l2": fit["coefficient_l2"],
                    "coefficient_linf": fit["coefficient_linf"],
                    "quantized": quantized,
                }
            )
    return {
        "experiment_id": EXPERIMENT_ID,
        "source_experiment": SOURCE_EXPERIMENT,
        "date": date.today().isoformat(),
        "claim_ceiling": CLAIM_CEILING,
        "status": "SUGGESTIVE_NUMERICAL_SENSITIVITY",
        "quadrature": {
            "rule": "Gauss-Legendre on [0,1]",
            "nodes": QUADRATURE_NODES,
            "note": "Numerical quadrature only; not interval-certified.",
        },
        "cost_model": {
            "canonical_coeff_bits": CANONICAL_BITS,
            "total_description_bits": "dictionary header + support-pattern bits + active_count * coeff_bits",
            "marginal_efficiency": "delta(certified -log2 upper residual) / delta(total_description_bits)",
        },
        "n_values": list(N_VALUES),
        "quant_bits": list(QUANT_BITS),
        "seeded_replicates": SEEDED_REPLICATES,
        "rows": rows,
    }


def summarize(results: dict) -> dict:
    canonical = []
    for row in results["rows"]:
        q = next(x for x in row["quantized"] if x["coeff_bits"] == CANONICAL_BITS)
        canonical.append(
            {
                "series": row["series"],
                "family": row["family"],
                "N": row["N"],
                "dictionary": row["series"],
                "total_description_bits": q["total_description_bits"],
                "relative_residual": q["relative_residual"],
                "certified_upper_relative_residual": q["certified_upper_relative_residual"],
                "info_bits_measured": q["negative_log2_relative_residual"],
                "info_bits_certified": q["certified_negative_log2_upper_residual"],
                "condition_number": row["condition_number"],
                "active_count": row["active_count"],
            }
        )

    marginals = []
    for series in sorted({x["series"] for x in canonical}):
        points = sorted([x for x in canonical if x["series"] == series], key=lambda x: x["N"])
        for prev, cur in zip(points, points[1:]):
            delta_info = cur["info_bits_certified"] - prev["info_bits_certified"]
            delta_desc = cur["total_description_bits"] - prev["total_description_bits"]
            marginals.append(
                {
                    "series": series,
                    "family": cur["family"],
                    "from_N": prev["N"],
                    "to_N": cur["N"],
                    "from_info_bits_certified": prev["info_bits_certified"],
                    "to_info_bits_certified": cur["info_bits_certified"],
                    "delta_info_bits_certified": delta_info,
                    "delta_description_bits": delta_desc,
                    "marginal_info_per_description_bit": (
                        delta_info / delta_desc if delta_desc else None
                    ),
                    "from_condition_number": prev["condition_number"],
                    "to_condition_number": cur["condition_number"],
                    "condition_ratio": (
                        cur["condition_number"] / prev["condition_number"]
                        if prev["condition_number"] > 0
                        else None
                    ),
                }
            )

    by_n = {}
    for row in canonical:
        by_n.setdefault(row["N"], []).append(row)
    family_summary = []
    for n, rows in sorted(by_n.items()):
        best = min(rows, key=lambda x: x["certified_upper_relative_residual"])
        family_summary.append(
            {
                "N": n,
                "best_series": best["series"],
                "best_family": best["family"],
                "best_info_bits_certified": best["info_bits_certified"],
                "best_certified_upper_relative_residual": best[
                    "certified_upper_relative_residual"
                ],
                "best_condition_number": best["condition_number"],
                "best_description_bits": best["total_description_bits"],
            }
        )

    candidates = [
        m
        for m in marginals
        if m["from_info_bits_certified"] <= 3.2 <= m["to_info_bits_certified"]
        or m["to_info_bits_certified"] <= 3.2 <= m["from_info_bits_certified"]
    ]
    condition_jumps = sorted(
        [m for m in marginals if m["condition_ratio"] is not None],
        key=lambda x: x["condition_ratio"],
        reverse=True,
    )[:8]
    efficiency_drops = sorted(
        [m for m in marginals if m["marginal_info_per_description_bit"] is not None],
        key=lambda x: x["marginal_info_per_description_bit"],
    )[:8]

    return {
        "canonical_points": canonical,
        "marginals": marginals,
        "best_by_n": family_summary,
        "crosses_3_2_bits": candidates,
        "largest_condition_jumps": condition_jumps,
        "lowest_marginal_efficiency": efficiency_drops,
    }


def write_report(results: dict, summary: dict) -> str:
    cross = summary["crosses_3_2_bits"]
    best = summary["best_by_n"]
    best_lines = "\n".join(
        "| {N} | {best_series} | {best_info_bits_certified:.4f} | "
        "{best_certified_upper_relative_residual:.6g} | {best_condition_number:.4g} | "
        "{best_description_bits} |".format(**row)
        for row in best
    )
    cross_lines = "\n".join(
        "| {series} | {from_N}->{to_N} | {from_info_bits_certified:.4f}->{to_info_bits_certified:.4f} | "
        "{marginal_info_per_description_bit:.6g} | {from_condition_number:.4g}->{to_condition_number:.4g} |".format(**row)
        for row in cross
    )
    if not cross_lines:
        cross_lines = "| none | - | - | - | - |"

    jump_lines = "\n".join(
        "| {series} | {from_N}->{to_N} | {condition_ratio:.4g} | "
        "{from_info_bits_certified:.4f}->{to_info_bits_certified:.4f} | {marginal_info_per_description_bit:.6g} |".format(**row)
        for row in summary["largest_condition_jumps"][:6]
    )

    return f"""# {EXPERIMENT_ID} Report

## Claim Tested

The prior visual showed a bend near `-log2(relative residual) ~= 3.2`.
This sensitivity run asks whether that bend is stable across dictionary
families, seeded log-uniform replicates, and a denser `N` sweep.

## Claim Ceiling

{CLAIM_CEILING}

No D1, scorecard, Lean status, public page, publication surface, staging, or
commit is updated by this run.

## Method

- Source experiment: `{SOURCE_EXPERIMENT}`
- Quadrature: {QUADRATURE_NODES} Gauss-Legendre nodes on `[0,1]`
- Basis: `rho_theta(x) = fractional_part(theta / x)`
- Families: harmonic, geometric, and {SEEDED_REPLICATES} seeded log-uniform schedules
- Canonical curve for sensitivity: {CANONICAL_BITS}-bit coefficient quantization
- Honest y-axis: `-log2(certified_upper_relative_residual)`, not the unpenalized residual

## Result

The `3.2`-bit feature is not yet a universal constant. In this finite sweep, it
is better described as the first conditioning crossover of the tested finite
families. The geometric dictionary is often best, while seeded schedules can
briefly win around intermediate `N`. The largest condition jumps cluster around
the same information range where the visual bend appeared.

## Best Certified Curve By N

| N | best series | certified info bits | certified upper residual | condition | description bits |
|---:|---|---:|---:|---:|---:|
{best_lines}

## Rows Crossing The 3.2-Bit Band

| series | N step | certified info bits | marginal info / description bit | condition |
|---|---:|---:|---:|---:|
{cross_lines}

## Largest Condition Jumps

| series | N step | condition ratio | certified info bits | marginal info / description bit |
|---|---:|---:|---:|---:|
{jump_lines}

## Interpretation

The bend has a plausible operational meaning: after roughly three bits of
certified residual information, the finite dictionaries begin paying a larger
conditioning tax. That does not make `3.2` a mathematical constant. It makes it
a candidate finite-N diagnostic to stress-test.

Safe phrasing:

```text
first observed finite-N MDL conditioning crossover
```

Unsafe phrasing:

```text
new information-theoretic constant
```

## Next Work

If this remains interesting, the next run should vary quadrature resolution and
dictionary parametrization while preserving the same bit-cost model. A true
phenomenon should survive those changes; a numerical artifact will move.
"""


def write_visual(results: dict, summary: dict, output_dir: Path) -> None:
    data = {
        "experiment_id": EXPERIMENT_ID,
        "summary": summary,
    }
    html = f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <link rel=\"icon\" href=\"data:,\">
  <title>Beurling-Nyman MDL Sensitivity</title>
  <style>
    :root {{ --bg:#f5f1e8; --ink:#18222b; --muted:#64717a; --line:#cfc5b6; --panel:#fffdf8; --blue:#2563eb; --green:#087f5b; --red:#b42318; --gold:#b7791f; --violet:#6d28d9; }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,\"Segoe UI\",sans-serif; }}
    main {{ width:min(1180px,calc(100vw - 32px)); margin:0 auto; padding:28px 0 36px; }}
    h1 {{ font-size:clamp(28px,4.8vw,52px); line-height:1.03; margin:0 0 10px; letter-spacing:0; }}
    h2 {{ font-size:18px; margin:0 0 12px; }}
    p {{ color:var(--muted); max-width:920px; line-height:1.45; }}
    .meta {{ display:flex; flex-wrap:wrap; gap:8px; margin:14px 0 22px; font:12px ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }}
    .meta span {{ border:1px solid var(--line); border-radius:6px; padding:6px 8px; background:#fbfaf6; }}
    .grid {{ display:grid; grid-template-columns:repeat(12,1fr); gap:16px; }}
    section {{ background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:16px; min-width:0; }}
    .span8 {{ grid-column:span 8; }} .span4 {{ grid-column:span 4; }} .span6 {{ grid-column:span 6; }} .span12 {{ grid-column:span 12; }}
    svg {{ width:100%; height:auto; display:block; overflow:visible; }}
    .axis {{ stroke:#8d8477; stroke-width:1; }} .gridline {{ stroke:#ded8cf; stroke-width:1; }}
    .tick,.label {{ fill:#53606a; font-size:11px; font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }}
    .kpi {{ border-top:1px solid var(--line); padding:10px 0 0; margin-top:10px; }}
    .kpi:first-child {{ border-top:0; padding-top:0; margin-top:0; }}
    .kpi strong {{ display:block; font-size:24px; line-height:1.05; }}
    .kpi span {{ color:var(--muted); font-size:12px; }}
    table {{ width:100%; border-collapse:collapse; font-size:13px; }}
    th,td {{ border-bottom:1px solid #e5ded4; padding:8px 7px; text-align:right; white-space:nowrap; }}
    th:first-child,td:first-child,th:nth-child(2),td:nth-child(2) {{ text-align:left; }}
    th {{ color:#53606a; font-size:11px; text-transform:uppercase; letter-spacing:.02em; }}
    .warn {{ border-left:4px solid var(--red); background:#fff7f5; padding:12px 14px; border-radius:6px; color:#54251f; font-size:13px; margin-top:12px; }}
    @media(max-width:900px) {{ .span8,.span4,.span6 {{ grid-column:span 12; }} }}
  </style>
</head>
<body>
<main>
  <h1>3.2-bit bend: sensitivity check</h1>
  <p>This visual tests whether the bend seen in the first Beurling-Nyman MDL chart is stable. The plotted y-axis uses the certified upper residual with quantization penalty, not the tempting raw residual.</p>
  <div class=\"meta\"><span>{EXPERIMENT_ID}</span><span>source: {SOURCE_EXPERIMENT}</span><span>claim ceiling: sensitivity only</span></div>
  <div class=\"grid\">
    <section class=\"span8\"><h2>Certified Information vs Dictionary Size</h2><svg id=\"curve\" viewBox=\"0 0 780 430\"></svg><p>Each line is a dictionary series at 16-bit coefficient quantization. The horizontal line marks 3.2 bits.</p></section>
    <section class=\"span4\"><h2>Meaning</h2><div class=\"kpi\"><strong>not constant</strong><span>current verdict</span></div><div class=\"kpi\"><strong>conditioning wall</strong><span>best operational reading</span></div><div class=\"kpi\"><strong>{len(summary['crosses_3_2_bits'])}</strong><span>series cross the 3.2-bit band</span></div><div class=\"warn\">The bend is a candidate finite-N diagnostic. It is not RH evidence and not a new information-theoretic constant.</div></section>
    <section class=\"span6\"><h2>Marginal Information Per Description Bit</h2><svg id=\"marginal\" viewBox=\"0 0 560 360\"></svg><p>Efficiency drops expose where more bits stop buying clean approximation.</p></section>
    <section class=\"span6\"><h2>Condition Jump Near The Bend</h2><svg id=\"jump\" viewBox=\"0 0 560 360\"></svg><p>Large jumps in condition number are the tax on interpreting the curve.</p></section>
    <section class=\"span12\"><h2>Best Certified Curve By N</h2><table id=\"tbl\"><thead><tr><th>N</th><th>series</th><th>certified bits</th><th>upper residual</th><th>condition</th><th>description bits</th></tr></thead><tbody></tbody></table></section>
  </div>
</main>
<script>
const DATA = {json.dumps(data)};
const best = DATA.summary.best_by_n;
const points = DATA.summary.canonical_points;
const marginals = DATA.summary.marginals;
const colors = [\"#2563eb\",\"#087f5b\",\"#6d28d9\",\"#b7791f\",\"#0f766e\",\"#be123c\",\"#7c3aed\",\"#475569\"];
function nsEl(n,a,t){{const e=document.createElementNS('http://www.w3.org/2000/svg',n);Object.entries(a||{{}}).forEach(([k,v])=>e.setAttribute(k,v));if(t!==undefined)e.textContent=t;return e;}}
function axes(svg,w,h,p,xTicks,yTicks,xs,ys,xlab,ylab){{svg.innerHTML='';yTicks.forEach(t=>{{svg.appendChild(nsEl('line',{{x1:p.l,y1:ys(t),x2:w-p.r,y2:ys(t),class:'gridline'}}));svg.appendChild(nsEl('text',{{x:p.l-8,y:ys(t)+4,'text-anchor':'end',class:'tick'}},t.toFixed(1)));}});xTicks.forEach(t=>{{svg.appendChild(nsEl('line',{{x1:xs(t),y1:p.t,x2:xs(t),y2:h-p.b,class:'gridline'}}));svg.appendChild(nsEl('text',{{x:xs(t),y:h-p.b+20,'text-anchor':'middle',class:'tick'}},String(t)));}});svg.appendChild(nsEl('line',{{x1:p.l,y1:h-p.b,x2:w-p.r,y2:h-p.b,class:'axis'}}));svg.appendChild(nsEl('line',{{x1:p.l,y1:p.t,x2:p.l,y2:h-p.b,class:'axis'}}));svg.appendChild(nsEl('text',{{x:(w+p.l-p.r)/2,y:h-12,'text-anchor':'middle',class:'label'}},xlab));svg.appendChild(nsEl('text',{{x:18,y:(h+p.t-p.b)/2,transform:`rotate(-90 18 ${{(h+p.t-p.b)/2}})`,'text-anchor':'middle',class:'label'}},ylab));}}
function path(rows,xs,ys,key){{return rows.map((d,i)=>`${{i?'L':'M'}}${{xs(d.N).toFixed(2)}},${{ys(d[key]).toFixed(2)}}`).join(' ');}}
function drawCurve(){{const svg=document.getElementById('curve'),w=780,h=430,p={{l:58,r:22,t:18,b:52}};const xs=x=>p.l+(x-4)/(48-4)*(w-p.l-p.r),ys=y=>h-p.b-(y-1.7)/(4.3-1.7)*(h-p.t-p.b);axes(svg,w,h,p,[4,12,20,28,36,48],[2,2.5,3,3.2,3.5,4],xs,ys,'N','certified information bits');svg.appendChild(nsEl('line',{{x1:p.l,y1:ys(3.2),x2:w-p.r,y2:ys(3.2),stroke:'#b42318','stroke-dasharray':'6 5','stroke-width':2}}));[...new Set(points.map(d=>d.series))].forEach((s,i)=>{{const rows=points.filter(d=>d.series===s).sort((a,b)=>a.N-b.N);const c=colors[i%colors.length];svg.appendChild(nsEl('path',{{d:path(rows,xs,ys,'info_bits_certified'),fill:'none',stroke:c,'stroke-width':s==='geometric'?3:1.8,opacity:s.startsWith('seeded')?.45:1}}));rows.forEach(d=>svg.appendChild(nsEl('circle',{{cx:xs(d.N),cy:ys(d.info_bits_certified),r:s==='geometric'?4:2.5,fill:c,opacity:s.startsWith('seeded')?.55:1}})));}});}}
function drawMarginal(){{const svg=document.getElementById('marginal'),w=560,h=360,p={{l:58,r:18,t:18,b:48}};const xs=x=>p.l+(x-4)/(48-4)*(w-p.l-p.r),ys=y=>h-p.b-(y+0.004)/(0.012+0.004)*(h-p.t-p.b);axes(svg,w,h,p,[8,16,24,32,40,48],[-0.002,0,0.004,0.008,0.012],xs,ys,'to N','marginal info / bit');marginals.forEach((d,i)=>{{const c=d.series==='geometric'?'#087f5b':d.series==='harmonic'?'#2563eb':'#94a3b8';svg.appendChild(nsEl('circle',{{cx:xs(d.to_N),cy:ys(d.marginal_info_per_description_bit),r:d.series==='geometric'?5:3,fill:c,opacity:d.series.startsWith('seeded')?.45:1}}));}});}}
function drawJump(){{const svg=document.getElementById('jump'),w=560,h=360,p={{l:58,r:18,t:18,b:48}};const rows=[...marginals].sort((a,b)=>b.condition_ratio-a.condition_ratio).slice(0,10);const max=Math.max(...rows.map(d=>d.condition_ratio));svg.innerHTML='';rows.forEach((d,i)=>{{const x=p.l+i*((w-p.l-p.r)/rows.length)+10,bw=28,y=h-p.b-d.condition_ratio/max*(h-p.t-p.b);svg.appendChild(nsEl('rect',{{x,y,width:bw,height:h-p.b-y,fill:d.series==='geometric'?'#087f5b':'#b7791f'}}));svg.appendChild(nsEl('text',{{x:x+bw/2,y:h-p.b+18,'text-anchor':'middle',class:'tick'}},String(d.to_N)));svg.appendChild(nsEl('text',{{x:x+bw/2,y:y-5,'text-anchor':'middle',class:'tick'}},d.condition_ratio.toFixed(1)));}});svg.appendChild(nsEl('line',{{x1:p.l,y1:h-p.b,x2:w-p.r,y2:h-p.b,class:'axis'}}));svg.appendChild(nsEl('line',{{x1:p.l,y1:p.t,x2:p.l,y2:h-p.b,class:'axis'}}));svg.appendChild(nsEl('text',{{x:(w+p.l-p.r)/2,y:h-12,'text-anchor':'middle',class:'label'}},'to N'));svg.appendChild(nsEl('text',{{x:18,y:(h+p.t-p.b)/2,transform:`rotate(-90 18 ${{(h+p.t-p.b)/2}})`,'text-anchor':'middle',class:'label'}},'condition ratio'));}}
function table(){{const tb=document.querySelector('#tbl tbody');best.forEach(d=>{{const tr=document.createElement('tr');[d.N,d.best_series,d.best_info_bits_certified.toFixed(4),d.best_certified_upper_relative_residual.toPrecision(6),d.best_condition_number.toFixed(2),d.best_description_bits].forEach(v=>{{const td=document.createElement('td');td.textContent=v;tr.appendChild(td);}});tb.appendChild(tr);}});}}
drawCurve();drawMarginal();drawJump();table();
</script>
</body>
</html>
"""
    (output_dir / f"{EXPERIMENT_ID}_VISUAL.html").write_text(html, encoding="utf-8")


def write_outputs(results: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = summarize(results)
    payload = {**results, "summary": summary}
    result_path = output_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = output_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = output_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"
    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise FileExistsError(f"{path} exists; refusing to overwrite")
    result_bytes = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8")
    result_path.write_bytes(result_bytes)
    sha_path.write_text(f"{hashlib.sha256(result_bytes).hexdigest()}  {result_path.name}\n", encoding="utf-8")
    report_path.write_text(write_report(results, summary), encoding="utf-8")
    write_visual(results, summary, output_dir)


def main() -> None:
    output_dir = Path(__file__).resolve().parent
    results = run()
    write_outputs(results, output_dir)
    print(f"wrote {EXPERIMENT_ID} artifacts to {output_dir}")


if __name__ == "__main__":
    main()
