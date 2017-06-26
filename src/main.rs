extern crate nalgebra;
extern crate image;

use std::fs::File;
use std::path::Path;

/*
struct Grids<N> where N : std::fmt::Debug + std::cmp::PartialEq + std::marker::Copy + 'static {
    /// grids on the various levels
    levels : std::vec::Vec<nalgebra::DMatrix<N>>,
    /// Columns of the coarsest grid
    cx : usize,
    /// Rows of the coarsest grid
    cy : usize,
    /// Width of the discretized space
    sx : N,
    /// Height of the discretized space
    sy : N
}

impl<N> Grids<N> where N : nalgebra::Scalar {
    pub fn new(sx : N, sy : N, cx : usize, cy : usize, nlevels : usize) -> Grids<N> {
        Grids {
            levels = std::
        }
    }
}*/

/// Right-hand side of equations
fn f(x : f64, y : f64) -> f64 {
    //12. * x * x + 24. * y * y
    if (x - 0.425).abs() < 0.025 && (y - 0.5).abs() < 0.15 {
        100.
    } else if (x - 0.575).abs() < 0.025 && (y - 0.5).abs() < 0.15 {
        -100.
    } else if (x - 0.45).abs() < 0.05 && (y - 0.675).abs() < 0.025 {
        100.
    } else if (x - 0.55).abs() < 0.05 && (y - 0.675).abs() < 0.025 {
        -100.
    } else {
        0.
    }
}

/// The unknown function "g" (just used for boundary conditions)
fn g(x : f64, y : f64) -> f64 {
    //x.powf(4.) + 2. * y.powf(4.)
    //f64::sin(x * 2. * std::f64::consts::PI)
    0.
}


fn restrict_halfweight(fine : &nalgebra::DMatrix<f64>, coarse : &mut nalgebra::DMatrix<f64>) {
    let (cnr, cnc) = coarse.shape();
    let (fnr, fnc) = fine.shape();

    assert_eq!((fnr - 1) / 2, cnr - 1);
    assert_eq!((fnc - 1) / 2, cnc - 1);

    // restrict center of matrix
    for j in 2..fnc - 2 {
        for i in 2..fnr - 2 {
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

fn interpolate(coarse : &nalgebra::DMatrix<f64>, fine : &mut nalgebra::DMatrix<f64>) {
    let (cnr, cnc) = coarse.shape();
    let (fnr, fnc) = fine.shape();

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

fn poisson_relax(u : &mut nalgebra::DMatrix<f64>, f : &nalgebra::DMatrix<f64>, h : f64) {
    assert_eq!(u.shape(), f.shape());

    let (nr, nc) = u.shape();
    let h2 = h*h;

    for j in 1..nc - 1 {
        for i in 1..nr-1 {
            u[(i,j)] = -0.25 * (h2 * f[(i,j)] - u[(i-1,j)] - u[(i+1,j)] - u[(i,j-1)] - u[(i,j+1)])
        }
    }
}

fn poisson_residual(u : & nalgebra::DMatrix<f64>, f : & nalgebra::DMatrix<f64>, h : f64) -> nalgebra::DMatrix<f64> {
    assert_eq!(u.shape(), f.shape());

    let (nr, nc) = u.shape();

    let mut r = nalgebra::DMatrix::from_element(nr, nc, 0.);

    let h2inv = 1. / (h * h);

    for j in 1..nc-1 {
        for i in 1..nr-1 {
            r[(i,j)] = f[(i,j)] - h2inv * (u[(i-1,j)] + u[(i+1,j)] + u[(i,j-1)] + u[(i,j+1)] - 4. * u[(i,j)])
        }
    }
    r
}

fn cycle(u : &mut nalgebra::DMatrix<f64>, f : & nalgebra::DMatrix<f64>, h : f64) {
    assert_eq!(u.shape(), f.shape());

    let (nr, nc) = u.shape();

    let (pre_iter, post_iter) = (10, 10);

    if nr == 3 || nc == 3 {
        // solve exact (probably only works on square matrices for now)
        poisson_relax(u, f, h)
    } else {
        for _ in 0..pre_iter {
            poisson_relax(u, f, h)
        }

        let r = poisson_residual(u, f, h);
        let mut rc = nalgebra::DMatrix::from_element(nr / 2 + 1, nc / 2 + 1, 0.);
        restrict_halfweight(& r, &mut rc);

        let mut uc = nalgebra::DMatrix::from_element(nr / 2 + 1, nc / 2 + 1, 0.);
        cycle(&mut uc, & rc, h * 2.);

        let mut uf = nalgebra::DMatrix::from_element(nr, nc, 0.);
        interpolate(&uc, &mut uf);

        for i in 0..nr {
            for j in 0..nc {
                u[(i,j)] += uf[(i,j)]
            }
        }

        for _ in 0..post_iter {
            poisson_relax(u, f, h)
        }
    }
}

fn poisson_err(u : & nalgebra::DMatrix<f64>, h : f64) -> f64 {
    let (nr, nc) = u.shape();

    let mut max_err = 0.;
    for i in 0..nr {
        for j in 0..nc {
            let err = (u[(i,j)] - g(j as f64 * h, i as f64 * h)).abs();
            max_err = err.max(max_err);
        }
    }
    max_err
}

fn save_plot(data : &nalgebra::DMatrix<f64>, filename : &str) {
    let max_val = data.iter().fold(std::f64::MIN, |acc, x| x.max(acc));
    let min_val = data.iter().fold(std::f64::MAX, |acc, x| x.min(acc));

    if max_val > min_val {
        let imgbuf = image::ImageBuffer::from_fn(data.ncols() as u32, data.nrows() as u32, |x,y| {
//            image::Luma([(255. * (u_arr[(y as usize,x as usize)] - min_val) / (max_val - min_val)).trunc() as u8])
            let norm = (data[(y as usize,x as usize)] - min_val) / (max_val - min_val);
            image::Rgb([(255. * norm).trunc() as u8, (255. * (1. - norm)).trunc() as u8, (255. * (0.5 - (norm - 0.5).abs())).trunc() as u8])
        });

        let ref mut fout = File::create(&Path::new(filename)).unwrap();

        // We must indicate the image’s color type and what format to save as
        let _ = image::ImageRgb8(imgbuf).save(fout, image::PNG);
    }
}

fn main() {
    let n = 1024;
    let h = 1. / n as f64;

    let f_arr = nalgebra::DMatrix::from_fn(n+1,n+1, |i : usize, j : usize| { f(h * j as f64, h * i as f64) });
    let mut u_arr = nalgebra::DMatrix::from_element(n+1, n+1, 0.);

    for i in 0..n+1 {
        u_arr[(0,i)] = g(i as f64 * h, 0.);
        u_arr[(n,i)] = g(i as f64 * h, n as f64 * h);
        u_arr[(i,0)] = g(0., i as f64 * h);
        u_arr[(i,n)] = g(n as f64 * h, i as f64 * h);
    }

    println!("err: {}", poisson_err(& u_arr, h));

    let niter = 20;

    for _ in 0..niter {
        cycle(&mut u_arr, & f_arr, h);
        println!("err: {}", poisson_err(& u_arr, h));
    }

    save_plot(& u_arr, "/tmp/multigrid_solution.png");
    save_plot(& f_arr, "/tmp/multigrid_rhs.png");
}
