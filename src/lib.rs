extern crate ndarray;
extern crate itertools;

/// The grid is represented using a two-dimensional array of 64-bit floating point values.
pub type Grid = ndarray::Array2<f64>;

/// The two types of cycles commonly used in multigrid algorithms.
#[derive(Debug, Copy, Clone)]
pub enum CycleType {
    /// The V-Cycle performs one recursive call to the next coarser level.
    V,
    /// The W-Cycle performs two recursive calls to the next coarser level.
    W,
    /// The F-Cycle performs one F-Cycle recursive call and one V-Cycle recursive call to the next coarser level.
    F
}

pub mod geometric;
