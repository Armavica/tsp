extern crate num;

#[cfg(test)]
#[macro_use]
extern crate quickcheck;

use num::Float;
use std::cmp;
use std::fmt;

/// Structure to represent a Travelling Salesman Problem.
pub struct TSP<N: Float> {
    distances: Vec<Vec<N>>,
}

impl<N: Float + fmt::Debug> TSP<N> {
    /// Returns the euclidian distance between two points in 2D.
    fn euc2d(a: (N, N), b: (N, N)) -> N {
        ((a.0 - b.0) * (a.0 - b.0) + (a.1 - b.1) * (a.1 - b.1)).sqrt()
    }
    /// Returns the euclidian distance between two points in 3D.
    fn euc3d(a: (N, N, N), b: (N, N, N)) -> N {
        ((a.0 - b.0) * (a.0 - b.0) + (a.1 - b.1) * (a.1 - b.1) + (a.2 - b.2) * (a.2 - b.2)).sqrt()
    }
    /// Exchanges two consecutive segments.
    fn exchange(path: &mut [usize], a: usize, b: usize, c: usize) {
        let ab = b - a;
        let bc = c - b;
        if ab == 0 || bc == 0 {
            return;
        }
        for i in 0..cmp::min(ab, bc) {
            path.swap(a + i + 1, b + i + 1);
        }
        if ab < bc {
            Self::exchange(path, b, b + ab, c);
        } else if ab > bc {
            Self::exchange(path, a + bc, b, c);
        }
    }
    /// Creates a problem from a list of 2D real coordinates.
    pub fn new_euc2d(vertices: &[(N, N)]) -> Self {
        let mut distances = Vec::new();
        for a in vertices {
            distances.push(vertices.iter().map(|&b| TSP::euc2d(*a, b)).collect());
        }
        TSP { distances: distances }
    }
    /// Creates a problem from a list of 3D real coordinates.
    pub fn new_euc3d(vertices: &[(N, N, N)]) -> Self {
        let mut distances = Vec::new();
        for a in vertices {
            distances.push(vertices.iter().map(|&b| TSP::euc3d(*a, b)).collect());
        }
        TSP { distances: distances }
    }
    /// Returns a 2-opt tour and its score, from an optionally specified
    /// starting tour.
    pub fn do_2opt(&self, start: Option<&[usize]>) -> (N, Vec<usize>) {
        let mut path = match start {
            None => (0..self.distances.len()).collect(),
            Some(path) => path.to_vec(),
        };
        let n = path.len();
        let mut stable = false;
        while !stable {
            stable = true;
            for a in 0..n {
                for b in a + 2..n {
                    let delta = self.distances[path[a]][path[b]] +
                        self.distances[path[a + 1]][path[(b + 1) % n]] -
                        self.distances[path[a]][path[a + 1]] -
                        self.distances[path[b]][path[(b + 1) % n]];
                    if delta < N::from(-1e-10).unwrap() {
                        stable = false;
                        path[a + 1..b + 1].reverse();
                    }
                }
            }
        }
        let mut distance = N::zero();
        for i in 0..path.len() {
            distance = distance + self.distances[path[i]][path[(i + 1) % n]];
        }
        (distance, path)
    }
    /// Returns a 3-opt tour and its score, from an optionally specified
    /// starting tour.
    pub fn do_3opt(&self, start: Option<&[usize]>) -> (N, Vec<usize>) {
        let mut path = match start {
            None => (0..self.distances.len()).collect(),
            Some(path) => path.to_vec(),
        };
        let n = path.len();
        let mut stable = false;
        while !stable {
            stable = true;
            for a in 0..n {
                for b in a + 1..n {
                    // 2-opt here?
                    for c in b + 1..n {
                        let current = self.distances[path[a]][path[a + 1]] +
                            self.distances[path[b]][path[b + 1]] +
                            self.distances[path[c]][path[(c + 1) % n]];
                        let swaps = vec![
                            // (a - a+1, b - c, b+1 - c+1) // 2-opt (a)
                            self.distances[path[a]][path[a + 1]] +
                                self.distances[path[b]][path[c]] +
                                self.distances[path[b + 1]][path[(c + 1) % n]],
                            // (a - c, b+1 - b, a+1 - c+1) // 2-opt (b)
                            self.distances[path[a]][path[c]] +
                                self.distances[path[b + 1]][path[b]] +
                                self.distances[path[a + 1]][path[(c + 1) % n]],
                            // (a - b, a+1 - b+1, c - c+1) // 2-opt (c)
                            self.distances[path[a]][path[b]] +
                                self.distances[path[a + 1]][path[b + 1]] +
                                self.distances[path[c]][path[(c + 1) % n]],
                            // (a - b+1, c - b, a+1 - c+1) // 3-opt (a - b+1)
                            self.distances[path[a]][path[b + 1]] +
                                self.distances[path[c]][path[b]] +
                                self.distances[path[a + 1]][path[(c + 1) % n]],
                            // (a - b, a+1 - c, b+1 - c+1) // 3-opt (c - a+1)
                            self.distances[path[a]][path[b]] +
                                self.distances[path[a + 1]][path[c]] +
                                self.distances[path[b + 1]][path[(c + 1) % n]],
                            // (a - c, b+1 - a+1, b - c+1) // 3-opt (b - c+1)
                            self.distances[path[a]][path[c]] +
                                self.distances[path[b + 1]][path[a + 1]] +
                                self.distances[path[b]][path[(c + 1) % n]],
                            // (a - b+1, c - a+1, b - c+1) // 3-opt (star)
                            self.distances[path[a]][path[b + 1]] +
                                self.distances[path[c]][path[a + 1]] +
                                self.distances[path[b]][path[(c + 1) % n]],
                        ];
                        let mut min = 0;
                        for i in 1..swaps.len() {
                            if swaps[i] < swaps[min] {
                                min = i;
                            }
                        }
                        if swaps[min] - current < N::from(-1e-10).unwrap() {
                            stable = false;
                            match min {
                                0 => path[b + 1..c + 1].reverse(),
                                1 => path[a + 1..c + 1].reverse(),
                                2 => path[a + 1..b + 1].reverse(),
                                3 => {
                                    path[a + 1..b + 1].reverse();
                                    Self::exchange(&mut path, a, b, c);
                                }
                                4 => {
                                    path[a + 1..b + 1].reverse();
                                    path[b + 1..c + 1].reverse();
                                }
                                5 => {
                                    path[b + 1..c + 1].reverse();
                                    Self::exchange(&mut path, a, b, c);
                                }
                                6 => Self::exchange(&mut path, a, b, c),
                                _ => unreachable!(),
                            }
                        }
                    }
                }
            }
        }
        let mut distance = N::zero();
        for i in 0..path.len() {
            distance = distance + self.distances[path[i]][path[(i + 1) % n]];
        }
        (distance, path)
    }
}

#[cfg(test)]
mod tests {
    use super::TSP;
    quickcheck! {
        fn exchange(v1: Vec<usize>, v2: Vec<usize>, v3: Vec<usize>, v4: Vec<usize>) -> bool {
            if v1.is_empty() {
                return true
            }
            let mut v1234 = Vec::new();
            v1234.extend(v1.clone());
            v1234.extend(v2.clone());
            v1234.extend(v3.clone());
            v1234.extend(v4.clone());
            let mut v1324 = Vec::new();
            v1324.extend(v1.clone());
            v1324.extend(v3.clone());
            v1324.extend(v2.clone());
            v1324.extend(v4.clone());
            let a = v1.len() - 1;
            let b = a + v2.len();
            let c = b + v3.len();
            TSP::<f64>::exchange(&mut v1234, a, b, c);
            v1234 == v1324
        }
    }
    #[test]
    fn pair() {
        let tsp = TSP::new_euc2d(&[(0., 0.), (4., 2.)]);
        assert_eq!(tsp.do_2opt(None), (2. * (20f64).sqrt(), vec![0, 1]));
    }
    #[test]
    fn square() {
        let tsp = TSP::new_euc2d(&[(0., 0.), (1., 1.), (0., 1.), (1., 0.)]);
        assert_eq!(tsp.do_2opt(None), (4., vec![0, 2, 1, 3]));
    }
    quickcheck! {
        fn no_point_missed_2opt(vs: Vec<(f64, f64)>) -> bool {
            let tsp = TSP::new_euc2d(&vs);
            let (_, mut path) = tsp.do_2opt(None);
            path.sort();
            path == (0..vs.len()).collect::<Vec<usize>>()
        }
    }
    quickcheck! {
        fn path_starts_at_base_2opt(vs: Vec<(f64, f64)>) -> bool {
            let tsp = TSP::new_euc2d(&vs);
            let (_, path) = tsp.do_2opt(None);
            match path.first() {
                None | Some(&0) => true,
                Some(_) => false,
            }
        }
    }
    quickcheck! {
        fn no_point_missed_3opt(vs: Vec<(f64, f64)>) -> bool {
            let tsp = TSP::new_euc2d(&vs);
            let (_, mut path) = tsp.do_3opt(None);
            path.sort();
            path == (0..vs.len()).collect::<Vec<usize>>()
        }
    }
    quickcheck! {
        fn path_starts_at_base_3opt(vs: Vec<(f64, f64)>) -> bool {
            let tsp = TSP::new_euc2d(&vs);
            let (_, path) = tsp.do_3opt(None);
            match path.first() {
                None | Some(&0) => true,
                Some(_) => false,
            }
        }
    }
}
