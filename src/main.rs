extern crate ndarray;
extern crate image;
extern crate multigrid;

use multigrid::CycleType;
use multigrid::geometric::problems::Problem2D;

use ndarray::Array2;
use std::fs::File;
use std::path::Path;

/// Right-hand side of equations
fn f(x : f64, y : f64) -> f64 {
    12. * x * x + 24. * y * y

    // if (x - 0.425).abs() < 0.025 && (y - 0.5).abs() < 0.15 {
    //     100.
    // } else if (x - 0.575).abs() < 0.025 && (y - 0.5).abs() < 0.15 {
    //     -100.
    // } else if (x - 0.45).abs() < 0.05 && (y - 0.675).abs() < 0.025 {
    //     100.
    // } else if (x - 0.55).abs() < 0.05 && (y - 0.675).abs() < 0.025 {
    //     -100.
    // } else {
    //     0.
    // }
}

/// The unknown function "g" (just used for boundary conditions)
fn g(x : f64, y : f64) -> f64 {
    x.powi(4) + 2. * y.powi(4)
    //f64::sin(x * 2. * std::f64::consts::PI)
    //0.
}

fn save_plot(data : &Array2<f64>, filename : &str) {
    let max_val = data.iter().fold(std::f64::MIN, |acc, x| x.max(acc));
    let min_val = data.iter().fold(std::f64::MAX, |acc, x| x.min(acc));

    if max_val > min_val {
        let (nrows, ncols) = data.dim();
        let imgbuf = image::ImageBuffer::from_fn(ncols as u32, nrows as u32, |x,y| {
            let norm = (data[(y as usize,x as usize)] - min_val) / (max_val - min_val);
//            image::Rgb([(255. * norm).trunc() as u8, (255. * (1. - norm)).trunc() as u8, (255. * (0.5 - (norm - 0.5).abs())).trunc() as u8])
            image::Luma([(255. * norm).trunc() as u8])
        });

        let _ = image::DynamicImage::ImageLuma8(imgbuf).save(filename);
    }
}

fn main() {
    let n = 1024;
    let h = 1. / n as f64;

    let f_arr = Array2::from_shape_fn((n+1,n+1), |(i,j)| { f(h * j as f64, h * i as f64) });
    let mut u_arr = Array2::zeros((n+1, n+1));

    for i in 0..n+1 {
        u_arr[(0,i)] = g(i as f64 * h, 0.);
        u_arr[(n,i)] = g(i as f64 * h, n as f64 * h);
        u_arr[(i,0)] = g(0., i as f64 * h);
        u_arr[(i,n)] = g(n as f64 * h, i as f64 * h);
    }


    // maximum number of iterations before stopping
    let max_iter = 40;

    // stop iteration when maximal residual is below that value
    let conv_stop = 0.0000005;

    let cfg = multigrid::geometric::Settings
    {
        pre_relax: 3,
        post_relax: 3,
    };

    // solve Poisson's equation on unit square
    let prob = multigrid::geometric::problems::PoissonProblem::unit_square();

    println!("res: {}", prob.max_residual(& u_arr, & f_arr));

    for cur_iter in 0..max_iter {
        multigrid::geometric::cycle2d(CycleType::F, &cfg, &prob, &mut u_arr, &f_arr);
        let res = prob.max_residual(& u_arr, & f_arr);
        println!("res: {}", res);
        if res < conv_stop {
            println!("converged after {} iterations", cur_iter + 1);
            break;
        }
    }

    // save solution as image
    // save_plot(& u_arr, "/tmp/multigrid_solution.png");
    // save_plot(& f_arr, "/tmp/multigrid_rhs.png");
}
