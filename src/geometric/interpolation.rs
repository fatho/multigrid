use ndarray::Array2;
use util::shape2;

pub fn interpolate(coarse : &Array2<f64>, fine : &mut Array2<f64>) {
    let (cnr, cnc) = shape2(coarse);
    let (fnr, fnc) = shape2(fine);

    assert_eq!((fnr - 1) / 2, cnr - 1);
    assert_eq!((fnc - 1) / 2, cnc - 1);

    // copy exact data points
    fine.slice_with_steps_mut((0,0),(cnr,cnc),(2,2)).copy_from(coarse);

    // interpolate everything but right/bottom boundary
    for j in 0..cnc - 1 {
        for i in 0..cnr - 1 {
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
