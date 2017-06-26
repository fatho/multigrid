use Grid;

/// Describes a two-dimensional multigrid problem.
pub trait Problem2D {
    /// Directly solves the given system of equations on the coarsest grid.
    fn solve(&self, solution : &mut Grid, rhs : & Grid);
    /// Performs a single relaxation step on the system of equations.
    fn relax(&self, solution : &mut Grid, rhs : & Grid);
    /// Computes the residual of the system of equations.
    fn residual(&self, solution : &Grid, rhs : &Grid, residual : &mut Grid);
    /// Computes the maximum residual of the system of equations. Used as a convergence criterion.
    fn max_residual(&self, solution : &Grid, rhs : &Grid) -> f64;
}

/// Describes the discretization of Poisson's equation with Dirichlet boundary conditions as a multigrid problem.
#[derive(Debug, Copy, Clone)]
pub struct PoissonProblem {
    /// size of the source space along x-axis
    pub width : f64,
    /// size of the source space along y-axis
    pub height : f64
}

impl PoissonProblem {
    pub fn new(width : f64, height : f64) -> PoissonProblem {
        PoissonProblem { width: width, height: height }
    }

    pub fn unit_square() -> PoissonProblem {
        Self::new(1., 1.)
    }

    fn h(&self, nx : usize, ny : usize) -> (f64, f64) {
        (self.width / (nx - 1) as f64, self.height / (ny - 1) as f64)
    }
}

/// FIXME: de-duplicate code

/// This implementation assumes a five-point stencil discretization of the 2D space
impl Problem2D for PoissonProblem {
    fn solve(&self, solution : &mut Grid, rhs : & Grid) {
        self.relax(solution, rhs);
    }

    fn relax(&self, solution : &mut Grid, rhs : & Grid) {
        assert_eq!(solution.shape(), rhs.shape());

        let (nr, nc) = solution.dim();
        let (hx, hy) = self.h(nc, nr);
        let hxx = hx.powi(2);
        let hyy = hy.powi(2);
        let hxyi = - 1. / (2. / hxx + 2. / hyy);
        let ihxx = 1. / hxx;
        let ihyy = 1. / hyy;

        // red-black iteration
        // by processing the data in a checkerboard pattern we can avoid data-dependencies between cells of the same color
        for i in 1..nr-1 {
            let mut j = 2 - (i & 1);
            while j < nc - 1 {
                solution[(i,j)] = hxyi * (rhs[(i,j)] - ihxx * (solution[(i,j-1)] + solution[(i,j+1)]) - ihyy * (solution[(i-1,j)] + solution[(i+1,j)]));
                j += 2;
            }
        }
        for i in 1..nr-1 {
            let mut j = 1 + (i & 1);
            while j < nc - 1 {
                solution[(i,j)] = hxyi * (rhs[(i,j)] - ihxx * (solution[(i,j-1)] + solution[(i,j+1)]) - ihyy * (solution[(i-1,j)] + solution[(i+1,j)]));
                j += 2;
            }
        }
    }

    fn residual(&self, solution : &Grid, rhs : &Grid, residual : &mut Grid) {
        assert_eq!(solution.shape(), rhs.shape());
        assert_eq!(solution.shape(), residual.shape());

        let (nr, nc) = solution.dim();
        let (hx, hy) = self.h(nc, nr);
        let hxx = hx.powi(2);
        let hyy = hy.powi(2);
        let hxy = 2. / hxx + 2. / hyy;
        let ihxx = 1. / hxx;
        let ihyy = 1. / hyy;

        for i in 1..nr-1 {
            for j in 1..nc-1 {
                residual[(i,j)] = rhs[(i,j)] - ihxx * (solution[(i,j-1)] + solution[(i,j+1)]) + hxy * solution[(i,j)] - ihyy * (solution[(i-1,j)] + solution[(i+1,j)]);
            }
        }
    }

    fn max_residual(&self, solution : &Grid, rhs : &Grid) -> f64 {
        assert_eq!(solution.shape(), rhs.shape());

        let (nr, nc) = solution.dim();
        let (hx, hy) = self.h(nc, nr);
        let hxx = hx.powi(2);
        let hyy = hy.powi(2);
        let hxy = 2. / hxx + 2. / hyy;
        let ihxx = 1. / hxx;
        let ihyy = 1. / hyy;

        let mut rmax : f64 = 0.;

        for i in 1..nr-1 {
            for j in 1..nc-1 {
                let r = rhs[(i,j)] - ihxx * (solution[(i,j-1)] + solution[(i,j+1)]) + hxy * solution[(i,j)] - ihyy * (solution[(i-1,j)] + solution[(i+1,j)]);
                rmax = rmax.max(r.abs());
            }
        }
        rmax
    }
}
