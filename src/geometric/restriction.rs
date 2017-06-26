use ndarray::Array2;

/// Halfweight-restriction from a `(2m+1) * (2n+1)` fine grid to a `(m+1) * (n+1)` coarse grid.
pub fn halfweight(fine : &Array2<f64>, coarse : &mut Array2<f64>) {
    let (cnr, cnc) = coarse.dim();
    let (fnr, fnc) = fine.dim();

    assert_eq!((fnr - 1) / 2, cnr - 1);
    assert_eq!((fnc - 1) / 2, cnc - 1);

    // restrict center of matrix
    for i in 2..fnr - 2 {
        for j in 2..fnc - 2 {
            coarse[(i/2,j/2)] = (4. * fine[(i,j)] + fine[(i-1,j)] + fine[(i+1,j)] + fine[(i,j-1)] + fine[(i,j+1)]) / 8.
        }
    }

    // restrict boundaries of matrix
    for i in 2..fnr - 2 {
        coarse[(i/2,0)] = (4. * fine[(i,0)] + fine[(i-1,0)] + fine[(i+1,0)] + fine[(i,1)]) / 7.;
        coarse[(i/2,cnc - 1)] = (4. * fine[(i,fnc-1)] + fine[(i-1,fnc-1)] + fine[(i+1,fnc-1)] + fine[(i,fnc-2)]) / 7.
    }
    for j in 2..fnc - 2 {
        coarse[(0,j/2)] = (4. * fine[(0,j)] + fine[(1,j)] + fine[(0,j-1)] + fine[(0,j+1)]) / 7.;
        coarse[(cnr-1,j/2)] = (4. * fine[(fnr-1,j)] + fine[(fnr-2,j)] + fine[(fnr-1,j-1)] + fine[(fnr-1,j+1)]) / 7.
    }

    // restrict corners of matrix
    coarse[(0,0)] = (4. * fine[(0,0)] + fine[(0,1)] + fine[(1,0)]) / 6.;
    coarse[(cnr-1,0)] = (4. * fine[(fnr-1,0)] + fine[(fnr-1,1)] + fine[(fnr-2,0)]) / 6.;
    coarse[(0,cnc-1)] = (4. * fine[(0,fnc-1)] + fine[(0,fnc-2)] + fine[(1,fnc-1)]) / 6.;
    coarse[(cnr-1,cnc-1)] = (4. * fine[(fnr-1,fnc-1)] + fine[(fnr-1,fnc-2)] + fine[(fnr-2,fnc-1)]) / 6.;
}

