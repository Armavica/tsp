extern crate tsp;
extern crate rand;
use rand::Rng;

fn main() {
    let mut avg2 = 0.;
    let mut avg3 = 0.;
    let mut avg2_3 = 0.;
    let n = 10000;
    for _ in 0..n {
        let mut v = Vec::new();
        let mut rng = rand::weak_rng();
        for _ in 0..100 {
            v.push((rng.gen::<f64>(), rng.gen::<f64>()));
        }
        let tsp = tsp::TSP::new_euc2d(&v);
        let (dist2, _) = tsp.do_2opt(None);
        let (dist3, _) = tsp.do_3opt(None);
        avg2 += dist2;
        avg3 += dist3;
        avg2_3 += (dist3 - dist2) / dist2;
    }
    println!("avg 2-opt distance: {}", avg2 / n as f64);
    println!("avg 3-opt distance: {}", avg3 / n as f64);
    println!(
        "avg 3-opt improvement % vs 2-opt: {}",
        100. * avg2_3 / n as f64
    );
}
