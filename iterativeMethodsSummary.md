---
math: true
permalink: /iterativeMethodsSummary
title: Iterative Methods Summary
layout: page
---

The iterative methods used in this assignment are all really interesting examples of applying dynamic programming to solve initial value problems. There is no clear "best method", as each one has different convergence conditions depending on the context of the initial value problem, but some have better overall performance.

The Explicit Euler technique is by far the simplest technique of those examined in this assignment. It is taught early on in calculus classes, and represents a single step of integration at a point. The Implicit Euler technique instead uses the unknown next value of the iteration in the same setup. This results in an equation that must be solved using some other method such as the Newton method for finding zeroes. Implicit Euler is more stable in general that Explicit Euler, and can be useful in a wider context.

Runge Kutta (a favorite of engineers) is favored due to its stability and ease of use. Runge Kutta techniques exist for arbitrary order, and thus are very useful when you need an easy way to get a high-order approximation of an initial value problem. This was probably the easiest to code (save for Explicit Euler) of the problems in this assignment.

My favorite method of this assignment, however, is the Predictor Corrector method using Adams Bashforth for an initial guess, then Adams Moulton for a correction. I think that it is fascinating how you can turn an equation that was implicit into an explicit one just by taking a good guess beforehand. This can save you so much work, as solving implicit equations can often be computationally costly and have problems with convergence (especially with initial value problems whose characteristic equations have roots with algebraic multiplicity \\( > 1 \\) ).
