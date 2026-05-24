//! Batch interval lemniscate-length oracle for EHP114 shape-cone packets.
//!
//! Input JSON:
//! {
//!   "degree": 14,
//!   "res": 220,
//!   "extent": 3.0,
//!   "points": [
//!     {"label": "...", "coeffs": [[re, im], ...]}  // ascending a_0..a_(n-1)
//!   ]
//! }
//!
//! Output JSON writes each point with `[inf, sup]` length enclosure. This is a
//! local verifier component, not a standalone proof of EHP114.

use inari::{interval, Interval};
use serde::{Deserialize, Serialize};
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Deserialize)]
struct InputPoint {
    label: String,
    coeffs: Vec<[f64; 2]>,
}

#[derive(Deserialize)]
struct BatchInput {
    degree: usize,
    res: usize,
    extent: f64,
    points: Vec<InputPoint>,
}

#[derive(Serialize)]
struct OutputPoint {
    label: String,
    length_lower: f64,
    length_upper: f64,
    elapsed_secs: f64,
}

#[derive(Serialize)]
struct BatchOutput {
    degree: usize,
    res: usize,
    extent: f64,
    point_count: usize,
    elapsed_secs: f64,
    points: Vec<OutputPoint>,
}

fn eval_poly_f64(z_re: f64, z_im: f64, degree: usize, coeffs: &[(f64, f64)]) -> (f64, f64) {
    let mut acc_re = 1.0;
    let mut acc_im = 0.0;
    for k in (0..degree).rev() {
        let t_re = acc_re * z_re - acc_im * z_im;
        let t_im = acc_re * z_im + acc_im * z_re;
        if k < coeffs.len() {
            acc_re = t_re + coeffs[k].0;
            acc_im = t_im + coeffs[k].1;
        } else {
            acc_re = t_re;
            acc_im = t_im;
        }
    }
    (acc_re, acc_im)
}

fn lemniscate_length_interval(
    degree: usize,
    coeffs: &[(f64, f64)],
    res: usize,
    extent: f64,
) -> Interval {
    let step_f = 2.0 * extent / res as f64;
    let step = interval!(step_f, step_f).unwrap();
    let mut total = interval!(0.0, 0.0).unwrap();

    for iy in 0..res {
        for ix in 0..res {
            let x0_f = -extent + ix as f64 * step_f;
            let y0_f = -extent + iy as f64 * step_f;
            let x1_f = x0_f + step_f;
            let y1_f = y0_f + step_f;

            let f = |x: f64, y: f64| -> f64 {
                let (pr, pi) = eval_poly_f64(x, y, degree, coeffs);
                pr * pr + pi * pi - 1.0
            };

            let fsw = f(x0_f, y0_f);
            let fse = f(x1_f, y0_f);
            let fne = f(x1_f, y1_f);
            let fnw = f(x0_f, y1_f);

            let case = ((fsw > 0.0) as u8)
                | (((fse > 0.0) as u8) << 1)
                | (((fne > 0.0) as u8) << 2)
                | (((fnw > 0.0) as u8) << 3);
            if case == 0 || case == 15 {
                continue;
            }

            let interp = |fa: f64, fb: f64| -> Interval {
                let d = fa - fb;
                if d.abs() < 1e-30 {
                    interval!(0.5, 0.5).unwrap()
                } else {
                    interval!((fa / d).clamp(0.0, 1.0), (fa / d).clamp(0.0, 1.0)).unwrap()
                }
            };

            let seg_interval =
                |t1x: Interval, t1y: Interval, t2x: Interval, t2y: Interval| -> Interval {
                    let dx = t1x - t2x;
                    let dy = t1y - t2y;
                    (dx * dx + dy * dy).sqrt()
                };

            let x0 = interval!(x0_f, x0_f).unwrap();
            let y0 = interval!(y0_f, y0_f).unwrap();
            let x1 = interval!(x1_f, x1_f).unwrap();
            let y1 = interval!(y1_f, y1_f).unwrap();

            let s = (x0 + interp(fsw, fse) * step, y0);
            let e = (x1, y0 + interp(fse, fne) * step);
            let n = (x0 + interp(fnw, fne) * step, y1);
            let w = (x0, y0 + interp(fsw, fnw) * step);

            let seg = |a: (Interval, Interval), b: (Interval, Interval)| -> Interval {
                seg_interval(a.0, a.1, b.0, b.1)
            };

            let cell_len = match case {
                1 | 14 => seg(s, w),
                2 | 13 => seg(s, e),
                3 | 12 => seg(w, e),
                4 | 11 => seg(e, n),
                5 => {
                    let avg = (fsw + fse + fne + fnw) / 4.0;
                    if avg > 0.0 {
                        seg(s, w) + seg(e, n)
                    } else {
                        seg(s, e) + seg(w, n)
                    }
                }
                6 | 9 => seg(s, n),
                7 | 8 => seg(w, n),
                10 => {
                    let avg = (fsw + fse + fne + fnw) / 4.0;
                    if avg > 0.0 {
                        seg(s, e) + seg(w, n)
                    } else {
                        seg(s, w) + seg(e, n)
                    }
                }
                _ => interval!(0.0, 0.0).unwrap(),
            };

            total = total + cell_len;
        }
    }

    total
}

fn parse_args() -> (PathBuf, PathBuf) {
    let args: Vec<String> = env::args().collect();
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--input" if i + 1 < args.len() => {
                input = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--output" if i + 1 < args.len() => {
                output = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            _ => {
                i += 1;
            }
        }
    }
    (
        input.expect("missing --input"),
        output.expect("missing --output"),
    )
}

fn quiet_requested() -> bool {
    env::args().any(|arg| arg == "--quiet")
}

fn main() {
    let (input_path, output_path) = parse_args();
    let raw = fs::read_to_string(&input_path).expect("failed to read input JSON");
    let input: BatchInput = serde_json::from_str(&raw).expect("failed to parse input JSON");
    let started = Instant::now();
    let mut points = Vec::with_capacity(input.points.len());

    for point in input.points {
        let coeffs: Vec<(f64, f64)> = point.coeffs.iter().map(|pair| (pair[0], pair[1])).collect();
        let point_started = Instant::now();
        let interval = lemniscate_length_interval(input.degree, &coeffs, input.res, input.extent);
        points.push(OutputPoint {
            label: point.label,
            length_lower: interval.inf(),
            length_upper: interval.sup(),
            elapsed_secs: point_started.elapsed().as_secs_f64(),
        });
    }

    let output = BatchOutput {
        degree: input.degree,
        res: input.res,
        extent: input.extent,
        point_count: points.len(),
        elapsed_secs: started.elapsed().as_secs_f64(),
        points,
    };
    let text = serde_json::to_string_pretty(&output).unwrap();
    fs::write(&output_path, text.clone() + "\n").expect("failed to write output JSON");
    if quiet_requested() {
        println!(
            "{{\"degree\":{},\"res\":{},\"point_count\":{},\"elapsed_secs\":{}}}",
            output.degree, output.res, output.point_count, output.elapsed_secs
        );
    } else {
        println!("{}", text);
    }
}
