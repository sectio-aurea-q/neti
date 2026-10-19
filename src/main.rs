//! eckart-ref command line
//!   eckart-ref reduced <v1> <v2> <u>             barriers in k_B T, u = h nu / k_B T
//!   eckart-ref si <V1 kJ/mol> <V2 kJ/mol> <nu cm^-1> <T K>
//!   eckart-ref grid <csv>                        rows: v1,v2,u  -> appends kappa
use eckart_ref::{kappa, kappa_si, Kappa};
use std::io::{BufRead, BufReader};

fn print(k: &Kappa) {
    println!(
        "{{\"kappa\": {:.17e}, \"quad_err\": {:.3e}, \"tail_bound\": {:.3e}, \"cos_branch\": {}, \"evals\": {}}}",
        k.value, k.quad_err, k.tail_bound, k.cos_branch, k.evals
    );
}

fn num(s: &str) -> f64 {
    s.parse().unwrap_or_else(|_| { eprintln!("not a number: {s}"); std::process::exit(2) })
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let r = match a.get(1).map(|s| s.as_str()) {
        Some("reduced") if a.len() == 5 => kappa(num(&a[2]), num(&a[3]), num(&a[4])),
        Some("si") if a.len() == 6 => kappa_si(num(&a[2]) * 1e3, num(&a[3]) * 1e3, num(&a[4]), num(&a[5])),
        Some("grid") if a.len() == 3 => {
            let f = std::fs::File::open(&a[2]).unwrap_or_else(|e| { eprintln!("{e}"); std::process::exit(2) });
            for line in BufReader::new(f).lines().map_while(Result::ok) {
                let p: Vec<&str> = line.split(',').collect();
                if p.len() < 3 || p[0].trim().parse::<f64>().is_err() { continue; }
                match kappa(num(p[0].trim()), num(p[1].trim()), num(p[2].trim())) {
                    Ok(k) => println!("{},{},{},{:.17e},{:.3e}", p[0].trim(), p[1].trim(), p[2].trim(), k.value, k.quad_err),
                    Err(e) => println!("{},{},{},ERR {}", p[0].trim(), p[1].trim(), p[2].trim(), e),
                }
            }
            return;
        }
        _ => {
            eprintln!("usage:\n  eckart-ref reduced <v1> <v2> <u>\n  eckart-ref si <V1 kJ/mol> <V2 kJ/mol> <nu cm^-1> <T K>\n  eckart-ref grid <csv: v1,v2,u>");
            std::process::exit(2)
        }
    };
    match r {
        Ok(k) => print(&k),
        Err(e) => { eprintln!("error: {e}"); std::process::exit(1) }
    }
}
