extern crate ndarray;
extern crate itertools;

mod util;

pub type Grid = ndarray::Array2<f64>;

pub enum CycleType {
    VCycle,
    WCycle,
    FCycle
}

pub mod geometric;
