use extendr_api::prelude::*;

#[extendr]
fn ridwgraph(power: f64) -> Doubles {
    let _ = power;
    Doubles::from_values([1.0])
}

extendr_module! {
    mod roadnet;
    fn ridwgraph;
}