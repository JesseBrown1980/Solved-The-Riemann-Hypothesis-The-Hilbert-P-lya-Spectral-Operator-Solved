//! Experiment builder: certified ordinates of the frozen operator at size K.
//!
//! Reads no zeta data. Writes an HBP receipt to stdout:
//!
//! ```text
//! cargo run --release --bin build_spectrum -- K [INERTIA=full|spot8|none]
//! ```

use riemann_hilbert_polya_spectral::family::{
    build, gram_first_row, inertia_of_shifted_gram, ordinate_enclosure, pow_mod_u64,
    sorted_positive_ordinates, FrozenOperator,
};
use riemann_hilbert_polya_spectral::interval::SCALE;
use std::sync::{Arc, Mutex};

type InertiaResult = (usize, i64, Result<usize, String>);

fn algebra_checks(op: &FrozenOperator) -> Vec<String> {
    let (k, p, w) = (op.k as u64, op.p, op.w);
    let mut out = Vec::new();
    let order_ok =
        pow_mod_u64(w, k, p) == 1 && (1..k).all(|e| k % e != 0 || pow_mod_u64(w, e, p) != 1);
    out.push(format!("ALGEBRA w_order_is_K={}", order_ok as u8));
    let skew = (0..op.k).all(|d| op.c[d] == -op.c[(op.k - d) % op.k]);
    out.push(format!("ALGEBRA integer_skew={}", skew as u8));
    let row_sum: i128 = op.c.iter().map(|&x| x as i128).sum();
    out.push(format!("ALGEBRA row_sum={row_sum}"));
    // DC of prism(Tv) mod P for a fixed test vector v_i = i + 1
    let mut dc: i128 = 0;
    for i in 0..op.k {
        for j in 0..op.k {
            dc += op.c[(j + op.k - i) % op.k] as i128 * (j as i128 + 1);
        }
    }
    out.push(format!(
        "ALGEBRA prism_dc_of_Tv={}",
        dc.rem_euclid(p as i128)
    ));
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let k: usize = args.get(1).expect("K").parse().expect("K integer");
    let mode = args.get(2).map(String::as_str).unwrap_or("full");
    let op = build(k);
    let m = op.m();
    println!("HBP FROZEN-SPECTRUM v1");
    println!("K={k} P={} G={} W={} M={m} SCALE={SCALE}", op.p, op.g, op.w);
    for line in algebra_checks(&op) {
        println!("{line}");
    }
    for d in 1..=m {
        println!("C d={d} c={}", op.c[d]);
    }
    println!("GRAM a0={}", gram_first_row(&op)[0]);
    for j in 1..=m {
        let (lo, hi) = ordinate_enclosure(&op, j);
        println!("GAMMA j={j} lo={lo} hi={hi}");
    }
    let ords = sorted_positive_ordinates(&op).unwrap_or_else(|e| panic!("{e}"));
    let mut disjoint = 0usize;
    for (r, o) in ords.iter().enumerate() {
        let sep = if r + 1 < ords.len() {
            o.hi < ords[r + 1].lo
        } else {
            true
        };
        disjoint += sep as usize;
        println!(
            "ORD rank={} j={} lo={} hi={} separated={}",
            r + 1,
            o.j,
            o.lo,
            o.hi,
            sep as u8
        );
    }
    println!("ORD_SUMMARY count={m} separated={disjoint}");

    // Inertia cross-check: τ strictly between ordinate r and r+1 (squared, integer units).
    let mut ranks: Vec<usize> = match mode {
        "full" => (0..=m).collect(),
        "spot8" => (1..=8).map(|i| (m * i).div_ceil(8)).collect(),
        _ => Vec::new(),
    };
    if mode == "spot8" {
        ranks.insert(0, 0);
    }
    let mut jobs: Vec<(usize, i64)> = Vec::new();
    let mut skipped = 0usize;
    for &r in &ranks {
        let tau: i128 = if r == 0 {
            if ords[0].lo / SCALE >= 2 {
                1
            } else {
                skipped += 1;
                continue;
            }
        } else {
            let a = ords[r - 1].hi / SCALE + 1;
            let b = if r < m { ords[r].lo / SCALE } else { a + 2 };
            if a >= b {
                skipped += 1;
                continue;
            }
            a * b
        };
        jobs.push((r, i64::try_from(tau).expect("tau fits i64")));
    }
    let results: Arc<Mutex<Vec<InertiaResult>>> = Arc::new(Mutex::new(Vec::new()));
    let queue = Arc::new(Mutex::new(jobs.clone()));
    let workers: usize = std::env::var("WORKERS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2);
    let mut handles = Vec::new();
    for _ in 0..workers {
        let (q, res, op) = (Arc::clone(&queue), Arc::clone(&results), op.clone());
        handles.push(std::thread::spawn(move || loop {
            let job = q.lock().unwrap().pop();
            let Some((r, tau)) = job else { break };
            let n = inertia_of_shifted_gram(&op, tau);
            eprintln!("K={} rank {r} done", op.k);
            res.lock().unwrap().push((r, tau, n));
        }));
    }
    for h in handles {
        h.join().unwrap();
    }
    let mut res = results.lock().unwrap().clone();
    res.sort_by_key(|x| x.0);
    let mut agree = 0usize;
    for (r, tau, n) in &res {
        match n {
            Ok(n) => {
                let ok = *n == 1 + 2 * r;
                agree += ok as usize;
                println!(
                    "INERTIA rank={r} tau={tau} neg={n} expected={} agree={}",
                    1 + 2 * r,
                    ok as u8
                );
            }
            Err(e) => println!("INERTIA rank={r} tau={tau} error={}", e.replace(' ', "_")),
        }
    }
    println!(
        "INERTIA_SUMMARY mode={mode} checked={} agree={agree} unseparated_skipped={skipped}",
        res.len()
    );
}
