use ndarray::Array2;

/// Interpolates from a `(m + 1) * (n + 1)` coarse grid to a `(2m + 1) * (2n + 1)` fine grid.
pub fn interpolate(coarse : &Array2<f64>, fine : &mut Array2<f64>) {
    let (cnr, cnc) = coarse.dim();
    let (fnr, fnc) = fine.dim();

    assert_eq!((fnr - 1) / 2, cnr - 1);
    assert_eq!((fnc - 1) / 2, cnc - 1);

    // copy exact data points
    for i in 0..cnr {
        for j in 0..cnc {
            fine[(i*2,j*2)] = coarse[(i,j)]
        }
    }

    // interpolate everything but right/bottom boundary
    for i in 0..cnr - 1 {
        for j in 0..cnc - 1 {
            fine[(i*2+1,j*2)] = 0.5 * (coarse[(i,j)] + coarse[(i+1,j)]);
            fine[(i*2,j*2+1)] = 0.5 * (coarse[(i,j)] + coarse[(i,j+1)]);
            fine[(i*2+1,j*2+1)] = 0.25 * (coarse[(i,j)] + coarse[(i,j+1)] + coarse[(i+1,j)] + coarse[(i+1,j+1)]);
        }
    }

    // interpolate bottom row
    for j in 0..cnc-1 {
        fine[(fnr-1,j*2+1)] = 0.5 * (coarse[(cnr-1,j)] + coarse[(cnr-1,j+1)]);
    }
    // interpolate right column
    for i in 0..cnr-1 {
        fine[(i*2+1,fnc-1)] = 0.5 * (coarse[(i,cnc-1)] + coarse[(i+1,cnc-1)]);
    }
}
