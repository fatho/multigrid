# multigrid

[![Build Status](https://travis-ci.org/fatho/multigrid.svg?branch=master)](https://travis-ci.org/fatho/multigrid)

A multigrid solver written in Rust. So far, it only supports geometric multigrid, but support for algebraic multigrid is planned.

The multigrid cycle is parameterized over the problem that is to be solved. An example for solving a discretized Poisson equation is included.
Every problem must implement the `Problem2D` trait providing methods for relaxation and computing the residual.
