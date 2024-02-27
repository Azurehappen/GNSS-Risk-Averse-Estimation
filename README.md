# Optimization-Based Risk-Averse State Estimation for GNSS outlier accommodation in the urban environment

Current targeted signal:
GPS: L1
GLO:
GAL:
BDS: B1I (Before customizing the satellite position calculation, please read BDS ICD; different signal channels may have different computations.)

The Risk-Averse Optimization problem is an multi-convex problem. The convex reformulation of the binary case is evaluated for achieving a globally optimal solution and used as a base comparison. The convex form can be solved by CVX but its runtime is not suitable for real-time applications. The Block Coordinate Method for both binary and non-binary cases results in a locally optimal solution with real-time feasibility and achieves comparable results with the globally optimal solution.

This repo is for the experiment of the PhD dissertation and DGNSS Risk-Aver paper. More information will be added later.
