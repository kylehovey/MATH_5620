---
math: true
permalink: /vonNeumann
title: Von Neumann Stability Analysis
layout: page
---

## Lax Wendroff Method (calculations from p. 213 in textbook)

# \\[ g(\xi) = 1 - \frac{1}{2}v(e^{i\xi h} - e^{-i\xi h}) + \frac{1}{2}v^2(e^{i\xi h} - 2 + e^{-i\xi h}) \\]

# \\[ = 1 - iv\sin(\xi h) + v^2(\cos(\xi h) - 1) \\]

# \\[ = 1 - iv \Big( 2\sin(\frac{\xi h}{2})cos(\frac{\xi h}{2} \Big) + v^2 \Big( 2sin^2(\frac{\xi h}{2}) \Big) \\]

# \\[ \implies |g(\xi)|^2 = \Big( 1 - 2sin^2(\frac{\xi h}{2}) \Big)^2 + 4\sin^2(\frac{\xi h}{2})cos^2(\frac{\xi h}{2} \\]

# \\[ = 1 - 4v^2(1 - v^2)sin^4(\frac{\xi h}{2}) \\]

# \\[ \implies \text{stable for } |v| \leq 1 \\]

## Beam Warming Method

# \\[ g(\xi) = 1 - \frac{v}{2}(3 - 4e^{-i\xi h} + e^{-i\xi 2h}) + \frac{v^2}{2}(1 - 2e^{-i\xi h} + e^{-i\xi 2h}) \\]

# \\[ \implies \text{stable for } 0 \leq v \leq 2 \\]
