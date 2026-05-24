# EHP114 n=14 Interval Taylor / M14 Packet

Experiment: `EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-M14-PACKET-20260505-01`

## Meaning

This packet builds the first interval Taylor scaffold for the local
mixed-remainder target:

```text
M14(eps,s) <= 12 eps^(1/14) + 50000 ||s||^2.
```

It evaluates interval lemniscate lengths on a structured Taylor stencil around
radial bases `r(eps) * roots_of_unity`, with `eps` in `[0.0001, 0.001, 0.01, 0.1]`.

This is a shadow signature, not universal law. It is not a solution of
Erdos #114 and not the uniform mixed-remainder theorem.

## Verdict

- Status: `INTERVAL_TAYLOR_MATRIX_FAIL`
- Shape basis rank: `24`
- Oracle point count: `4804`
- Taylor step: `0.004`
- Local axis cap: `0.008`
- Global Gershgorin interval lower bound: `-131824.7210585011`
- Working `lambda14`: `100000.0`
- All matrix epsilon slices pass lambda: `False`
- Axis endpoints passing budget: `True`
- Axis minimum margin: `3.5420727715548637`
- All axis endpoint roots admissible: `True`
- All oracle points root-admissible: `False`
- Max oracle root radius: `1.003016572595277`
- Max axis root radius: `0.9999992856811205`

## Per-Epsilon Taylor Matrix Summary

| eps | lambda pass | Gershgorin lower | max interval width | worst row |
|---:|---:|---:|---:|---:|
| 1e-04 | no | -74273.966431 | 1.46658e-07 | 20 |
| 1e-03 | no | -131824.721059 | 1.08579e-07 | 7 |
| 1e-02 | no | -42486.401566 | 7.72743e-08 | 1 |
| 1e-01 | no | -589.086951476 | 4.10783e-08 | 1 |

## Worst Axis Endpoint

```json
{
  "admissible_at_endpoint": true,
  "deficit_interval": {
    "hi": 12.957442386633552,
    "lo": 12.95744238663232,
    "width": 1.2327916465437738e-12
  },
  "deficit_lower": 12.95744238663232,
  "eps": 0.0001,
  "kind": "axis_cap",
  "label": "eps:1e-04:axis-cap:21:-",
  "length_interval": {
    "hi": 17.895468454916212,
    "lo": 17.895468454914994,
    "width": 1.2185807918285718e-12
  },
  "local_cap_active": true,
  "margin": 3.5420727715548637,
  "max_root_radius": 0.9999972019044991,
  "pass": true,
  "rhs_mixed_absorption_budget": 9.415369615077456,
  "shape_index": 21,
  "shape_label": "m6_sin_tangent",
  "sign": "-",
  "t": 0.008,
  "tmax": 0.010257388440467812
}
```

## Interpretation

The packet gives a split verdict:

```text
axis endpoint budget: passes
uniform positive Taylor matrix at radial bases: fails
```

The strongest safe reading is:

```text
The conservative M14 budget survives the tested admissible signed axis
endpoints, but the naive assumption that the boundary shape cone remains
uniformly positive after radial contraction is false on this stencil.
```

The unsafe reading is that the local cone is fully certified. This packet does
not yet bound every multi-mode shape vector inside the 24-dimensional shape
ball. The central-difference Taylor matrix is an ambient diagnostic; some
off-axis stencil points may lie just outside root admissibility. The admissible
axis endpoints are reported separately.

## Next Theorem Target

Upgrade this packet from stencil evidence to a box theorem:

```text
For every eps in (0, 1/10] and every shape vector s with
||s|| <= min(eta14Boundary(eps), 0.008),
M14(eps,s) <= 12 eps^(1/14) + 50000 ||s||^2.
```

The next run should add derivative/remainder interval bounds over multi-mode
boxes, not just axis endpoints.

## Guardrails

- No scorecard update.
- No D1 update.
- No public claim.
- No Lean status change.
- No email.
- No git operation.
