// ==================== NESTED MODULES ====================
pub mod restriction;
pub mod interpolation;
pub mod problems;

// ====================    IMPORTS     ====================
use Grid;
use CycleType;
use self::problems::Problem2D;

/* NOTE: The discretization grid resulting from ``n x n` intervals consists of `(n+1) x (n+1)` variables.
*/

/// Settings controlling the behavior of the multigrid solver.
#[derive(Debug, Copy, Clone)]
pub struct Settings {
    /// Number of relaxation steps prior to cycle
    pub pre_relax : usize,
    /// Number of relaxation steps after cycle
    pub post_relax : usize,
}

/// Performs a single multigrid cycle for a given problem and right-hand side, starting with an initial guess.
pub fn cycle2d<P>(cycle_type : CycleType, settings : &Settings, problem : &P, solution : &mut Grid, rhs : &Grid) where P : Problem2D {
    let (nr, nc) = solution.dim();

    if nr <= 3 || nc <= 3 || (nr - 1) % 2 != 0 || (nc - 1) % 2 != 0 {
        // we cannot coarsen the grid further, solve directly
        problem.solve(solution, rhs);
    } else {
        // 1. pre-smoothing
        for _ in 0..settings.pre_relax {
            problem.relax(solution, rhs);
        }

        // 2. compute restriction of residual to coarser grid
        let mut tmp = Grid::zeros(solution.dim());
        problem.residual(solution, rhs, &mut tmp);
        let mut coarse_rhs = Grid::zeros((nr / 2 + 1, nc / 2 + 1));
        // FIXME: make restriction operator configurable
        restriction::halfweight(& tmp, &mut coarse_rhs);

        // 3. compute coarse-grid corrections

        let mut coarse_solution = Grid::zeros((nr / 2 + 1, nc / 2 + 1));
        match cycle_type {
            CycleType::V => {
                cycle2d(cycle_type, settings, problem, &mut coarse_solution, & coarse_rhs)
            }
            CycleType::W => {
                cycle2d(cycle_type, settings, problem, &mut coarse_solution, & coarse_rhs);
                cycle2d(cycle_type, settings, problem, &mut coarse_solution, & coarse_rhs)
            }
            CycleType::F => {
                cycle2d(cycle_type, settings, problem, &mut coarse_solution, & coarse_rhs);
                cycle2d(CycleType::V, settings, problem, &mut coarse_solution, & coarse_rhs);
            }
        }

        // 4. interpolate and apply corrections to finer grid

        // FIXME: make interpolation operator configurable
        interpolation::interpolate(&coarse_solution, &mut tmp);
        *solution += &tmp;

        // 5. post-smoothing
        for _ in 0..settings.post_relax {
            problem.relax(solution, rhs);
        }
    }
}
