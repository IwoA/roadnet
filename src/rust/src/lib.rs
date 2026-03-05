use extendr_api::prelude::*;
use petgraph::algo::dijkstra;
use petgraph::graph::{NodeIndex, UnGraph};
use std::collections::HashMap;

#[extendr]
fn ridwgraph(
    segment_ids: Integers,
    centroids_x: Doubles,
    centroids_y: Doubles,
    start_x: Doubles,
    start_y: Doubles,
    end_x: Doubles,
    end_y: Doubles,
    initial_values: Doubles, 
    power: f64
) -> Doubles {
    let n = segment_ids.len();
    let mut endpoint_to_segments: HashMap<(u64, u64), Vec<usize>> = HashMap::new();

    // Coordinate hashing for spatial joining
    let hash_coord = |x: f64, y: f64| -> (u64, u64) {
        ((x * 1_000_000.0) as u64, (y * 1_000_000.0) as u64)
    };

    for i in 0..n {
        // Use .inner() to extract the f64 value from Rfloat
        let sx = start_x[i].inner();
        let sy = start_y[i].inner();
        let ex = end_x[i].inner();
        let ey = end_y[i].inner();

        for &(x, y) in &[(sx, sy), (ex, ey)] {
            endpoint_to_segments.entry(hash_coord(x, y)).or_default().push(i);
        }
    }

    // Build Dual Graph
    let mut graph = UnGraph::<usize, f64>::new_undirected();
    let nodes: Vec<NodeIndex> = (0..n).map(|i| graph.add_node(i)).collect();

    for (_, connected_indices) in endpoint_to_segments {
        if connected_indices.len() > 1 {
            for i in 0..connected_indices.len() {
                for j in i + 1..connected_indices.len() {
                    let idx1 = connected_indices[i];
                    let idx2 = connected_indices[j];
                    
                    // Use .inner() for mathematical operations
                    let dx = centroids_x[idx1].inner() - centroids_x[idx2].inner();
                    let dy = centroids_y[idx1].inner() - centroids_y[idx2].inner();
                    let dist = (dx.powi(2) + dy.powi(2)).sqrt();
                    
                    graph.add_edge(nodes[idx1], nodes[idx2], dist);
                }
            }
        }
    }

    // Source identification
    let sources: Vec<(NodeIndex, f64)> = (0..n)
        .filter(|&i| !initial_values[i].is_na())
        .map(|i| (nodes[i], initial_values[i].inner())) 
        .collect();

    // Pre-calculate Dijkstra paths from all sources
    let mut distances_from_sources = Vec::new();
    for &(source_node, value) in &sources {
        let dists = dijkstra(&graph, source_node, None, |e| *e.weight());
        distances_from_sources.push((dists, value));
    }

    // IDW Calculation
    let results: Vec<f64> = (0..n).map(|i| {
        if !initial_values[i].is_na() { 
            return initial_values[i].inner(); 
        }

        let (mut sum_w, mut sum_wv) = (0.0, 0.0);
        for (dists, val) in &distances_from_sources {
            if let Some(&d) = dists.get(&nodes[i]) {
                if d > 0.0 {
                    let w = 1.0 / d.powf(power);
                    sum_w += w;
                    sum_wv += w * val;
                } else {
                    return *val; // Direct hit on source
                }
            }
        }
        if sum_w > 0.0 { sum_wv / sum_w } else { f64::NAN }
    }).collect();

    Doubles::from_values(results)
}

extendr_module! {
    mod roadnet; // Changed from 'networkidw' to 'roadnet'
    fn ridwgraph;
}