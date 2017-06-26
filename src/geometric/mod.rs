// ==================== NESTED MODULES ====================
pub mod restriction;
pub mod interpolation;

// ====================    IMPORTS     ====================
use itertools::{iterate};
use ndarray::{Ix};
use std::vec::Vec;

use util;
use Grid;
use CycleType;

pub trait Problem2D {
    fn relax(&self, solution : &mut Grid, rhs : & Grid);
    fn residual(&self, &Grid, &Grid) -> Grid;
    fn restrict(&self, fine : &Grid, coarse : &mut Grid);
    fn interpolate(&self, coarse : &Grid, fine : &mut Grid);
}

pub struct PoissonProblem {
    /// number of discretization steps along x-axis
    nx     : usize,
    /// number of discretization steps along y-axis
    ny     : usize,
    /// size of the source space along the x-axis
    sx  : f64,
    /// size of the source space along the y-axis
    sy : f64
}

impl PoissonProblem {
    pub fn new(sx : f64, sy : f64) -> PoissonProblem {
        
    }
}

impl Problem2D for PoissonProblem {
    fn relax(&self, solution : &mut Grid, rhs : & Grid) {
        assert_eq!(solution.shape(), rhs.shape());

        let (nr, nc) = util::shape2(solution);
        let h2 = self.hx*self.hy;

        for j in 1..nc - 1 {
            for i in 1..nr-1 {
                u[(i,j)] = -0.25 * (h2 * f[(i,j)] - u[(i-1,j)] - u[(i+1,j)] - u[(i,j-1)] - u[(i,j+1)])
            }
        }
    }

    fn restrict(&self, fine : &Grid, coarse : &mut Grid) {
        restriction::halfweight(fine, coarse);
    }

    fn interpolate(&self, coarse : &Grid, fine : &mut Grid) {
        interpolation::interpolate(coarse, fine);
    }
}

pub struct Settings {
    pub pre_relax : usize,
    pub post_relax : usize,
    pub cycle_type : CycleType,
    pub num_cycles : usize,
}

pub fn cycle2d<P>(settings : &Settings, problem : &P) where P : Problem2D {

    for _ in 0..settings.pre_relax {
        
    }

}

/*
/// Data structure for managing the state of a multigrid solver.
pub struct Multigrid<P>
    where Restriction : Fn(&Grid, &mut Grid),
          Interpolation : Fn(&Grid, &mut Grid)
{
    /// stores the variable assignments
    solution : Vec<Grid>,
    /// stores the right hand side of the problem
    rhs : Vec<Grid>,
    /// temporary storage used internally
    tmp : Vec<Grid>,
    /// finest level of the multigrid
    finest_level : usize,
    /// width of the problem domain (in actual units)
    pub width : f64,
    /// height of the problem domain (in actual units)
    pub height : f64,
    pub restriction_op : Restriction
}

impl<R, I> Multigrid<R, I> {
    pub fn new(width : f64, height : f64, nxcoarse : usize, nycoarse : usize, levels : usize) -> Multigrid<R, I> {
        assert!(levels > 0, "multigrid needs at least one level");

        // note that the order for the shape is (y,x) as the indices are interpreted as [rows, columns]
        let empty_grids : Vec<_> =
            iterate((nycoarse, nxcoarse), |&(rows,cols)| (2 * (rows - 1) + 1, 2 * (cols - 1) + 1))
            .take(levels)
            .map(Array2::zeros)
            .collect();

        Multigrid {
            solution: empty_grids.clone(),
            rhs: empty_grids.clone(),
            tmp: empty_grids,
            finest_level: levels - 1,
            width: width,
            height: height
        }
    }

    /// Provides a mutable reference for modifying the right-hand side of the problem at the finest level.
    pub fn finest_rhs_mut(&mut self) -> &mut Grid {
        let lvl = self.finest_level;
        &mut self.rhs[lvl]
    }

    /// Returns the right-hand side of the problem at the finest level.
    pub fn finest_rhs(& self) -> & Grid {
        let lvl = self.finest_level;
        & self.rhs[lvl]
    }

    /// Provides a mutable reference for modifying the right-hand side of the problem at the finest level.
    pub fn finest_solution_mut(&mut self) -> &mut Grid {
        let lvl = self.finest_level;
        &mut self.solution[lvl]
    }

    /// Returns the right-hand side of the problem at the finest level.
    pub fn finest_solution(& self) -> & Grid {
        let lvl = self.finest_level;
        & self.solution[lvl]
    }

    fn level_shape(&self, level : usize) -> (usize, usize) {
        let s = self.solution[level].shape();
        (s[0], s[1])
    }

    fn real_coords(& self, level : usize, row : Ix, col : Ix) -> (f64,f64) {
        let (rows, cols) = self.level_shape(level);
        ( (col as f64) / ((cols - 1) as f64) * self.width,
          (row as f64) / ((rows - 1) as f64) * self.height)
    }

    /// Initializes the boundary solutions with a known function.
    pub fn set_dirichlet_boundary<F>(&mut self, f : F) where F : Fn((f64, f64)) -> f64 {
        let (nrows, ncols) = self.level_shape(self.finest_level);
        for i in 0..nrows {
            self.finest_solution_mut()[[i, 0]] = f(self.real_coords(self.finest_level, i, 0));
            self.finest_solution_mut()[[i, ncols - 1]] = f(self.real_coords(self.finest_level, i, ncols - 1));
        }
        // corners have already been handled by previous loop
        for j in 1..ncols-1 {
            self.finest_solution_mut()[[0, j]] = f(self.real_coords(self.finest_level, 0, j));
            self.finest_solution_mut()[[nrows - 1, j]] = f(self.real_coords(self.finest_level, nrows - 1, j));
        }
    }
}
*/
