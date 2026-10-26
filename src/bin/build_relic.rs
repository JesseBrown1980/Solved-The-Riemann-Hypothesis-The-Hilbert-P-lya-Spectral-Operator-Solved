//! Relic-landing builder: certified ordinates of the nested operator.
//!
//! Reads no zeta data. Writes an HBP receipt to stdout:
//!
//! ```text
//! cargo run --release --bin build_relic -- LEVEL [INERTIA=full|none]
//! ```

use riemann_hilbert_polya_spectral::interval::SCALE;
use riemann_hilbert_polya_spectral::nested::{
    build_with_prime_offset, dense_generator, factors, inertia_dense_gram, landing, landing_size,
    positive_ordinates,
};
use std::sync::{Arc, Mutex};

type Job = (usize, i64);
type Done = (usize, i64, Result<usize, String>);

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let level: usize = args.get(1).expect("LEVEL").parse().expect("LEVEL integer");
    let mode = args.get(2).map(String::as_str).unwrap_or("none");
    let l = landing(level);
    assert_eq!(l.size, landing_size(level));
    let m = (l.size - 1) / 2;
    println!("HBP RELIC-LANDING-SPECTRUM v1");
    println!("K={} M={m} LEVEL={level} SCALE={SCALE}", l.size);
    for (i, &(q, offset)) in factors(level).iter().enumerate() {
        let op = build_with_prime_offset(q, offset);
        let role = if i == 0 { "outer" } else { "inner" };
        let cs: Vec<String> = op.c[1..=op.m()].iter().map(|c| c.to_string()).collect();
        println!(
            "FACTOR index={i} role={role} K={q} P={} G={} W={} prime_offset={offset} c={}",
            op.p,
            op.g,
            op.w,
            cs.join(",")
        );
    }
    let pos = positive_ordinates(&l).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(pos.len(), m, "one free centre plus symmetric pairs");
    println!(
        "CENTRE zero_modes=1 positive={} negative={}",
        pos.len(),
        pos.len()
    );
    let mut separated = 0usize;
    for (r, &(lo, hi, idx)) in pos.iter().enumerate() {
        let sep = r + 1 == pos.len() || hi < pos[r + 1].0;
        separated += sep as usize;
        println!(
            "ORD rank={} j={idx} lo={lo} hi={hi} separated={}",
            r + 1,
            sep as u8
        );
    }
    println!("ORD_SUMMARY count={m} separated={separated}");
    if mode != "full" {
        println!("INERTIA_SUMMARY mode=structure checked=0 agree=0 note=kronecker_sum_of_inertia_validated_factors");
        return;
    }
    let (n, t) = dense_generator(level);
    assert_eq!(n, l.size);
    for i in 0..n {
        for j in 0..n {
            assert_eq!(t[i * n + j], -t[j * n + i], "dense generator is skew");
        }
    }
    println!("DENSE size={n} skew=1");
    let mut jobs: Vec<Job> = Vec::new();
    let mut skipped = 0usize;
    for r in 0..=m {
        let tau: i128 = if r == 0 {
            if pos[0].0 / SCALE >= 2 {
                1
            } else {
                skipped += 1;
                continue;
            }
        } else {
            let a = pos[r - 1].1 / SCALE + 1;
            let b = if r < m { pos[r].0 / SCALE } else { a + 2 };
            if a >= b {
                skipped += 1;
                continue;
            }
            a * b
        };
        jobs.push((r, i64::try_from(tau).expect("tau fits")));
    }
    let queue = Arc::new(Mutex::new(jobs));
    let done: Arc<Mutex<Vec<Done>>> = Arc::new(Mutex::new(Vec::new()));
    let t = Arc::new(t);
    let workers: usize = std::env::var("WORKERS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2);
    let handles: Vec<_> = (0..workers)
        .map(|_| {
            let (q, d, t) = (Arc::clone(&queue), Arc::clone(&done), Arc::clone(&t));
            std::thread::spawn(move || loop {
                let job = q.lock().unwrap().pop();
                let Some((r, tau)) = job else { break };
                let res = inertia_dense_gram(n, &t, tau);
                d.lock().unwrap().push((r, tau, res));
            })
        })
        .collect();
    for h in handles {
        h.join().unwrap();
    }
    let mut res = done.lock().unwrap().clone();
    res.sort_by_key(|x| x.0);
    let mut agree = 0usize;
    for (r, tau, out) in &res {
        match out {
            Ok(neg) => {
                let ok = *neg == 1 + 2 * r;
                agree += ok as usize;
                println!(
                    "INERTIA rank={r} tau={tau} neg={neg} expected={} agree={}",
                    1 + 2 * r,
                    ok as u8
                );
            }
            Err(e) => println!("INERTIA rank={r} tau={tau} error={}", e.replace(' ', "_")),
        }
    }
    println!(
        "INERTIA_SUMMARY mode=full checked={} agree={agree} unseparated_skipped={skipped}",
        res.len()
    );
}
