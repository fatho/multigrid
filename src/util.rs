use ndarray::{Array2, Ix};

// Extracts the shape of a two-dimensional array.
pub fn shape2<A>(arr : &Array2<A>) -> (Ix, Ix) {
    let s = arr.shape();
    (s[0], s[1])
}
