# MSF Analysis of Delayed Systems
Codes for MSF analysis for delayed coupled systems with examples in Stuart-Landau oscillators and the Lang-Kobayashi model for coupled lasers.
See the references below for more details.

# Usage 
- `StuartLandau` : This folder contains the routines to generate the MSF surface landscape (`LK_NonDelayed.m`, `LK_SemiDiffusive.m`, `LK_DelayedDiffusive.m`  ) and the study of the dependence of the region of stability in the coupling strength (`SL_DepthSize.m`) for the Stuart-Landau with identical indegrees.

- `LangKobayashi` : This folder contains the routines to generate the MSF surface landscape (`LK_NonDiffusive.m`) and the study of the dependence of the region of stability in the coupling strength (`LK_DepthSize.m`) for the Lang-Kobayashi with identical indegrees.

- `LyapExpAnalysis` : This folder contains the functions needed for the stability analysis, including the transcendental equations for Lang-Kobayashi homogeneous systems (`sync.m`) used to calculate the solutions of the transcendental equations for identical parameters with regular networks. The codes to solve the transcendental equations for heterogeneous systems (`sync_het_net.m`) used when the network is heterogeneous. This folder also includes the codes to calculate the Jacobian matrices for each model (`Jac_SL_hetnet.m `, `Jac_LK_hetnet.m`) used for the Stuart-Landau and Lang-Kobayashi models, respectively. Moreover, we also include the functions needed to calculate the rightmost Lyapunov exponent, `dde_rightmost_eig.m`, and its dependence `stst_stabil_mod`. The latter calculation requires the installation of the DDE-BIFTOOL toolbox. These codes were tested on version 3.0. See Ref. 2 and their website (https://sourceforge.net/projects/ddebiftool/) for more details.

- `NetworkOpt` : This folder contains optimization routines (`SL_NetOpt.m`, `LK_NetOpt.m` ) used to minimize the Lyapunov exponent of the Stuart-Landau and Lang-Kobayashi model, respectively associated with the corresponding constraints. We also include the code to sparsify the network (`Sparsify.m`) by removing edges with weight below a given cutoff and renormalizing the remaining edges.

- All codes were tested and run in MATLAB 2024b. To run the codes, download all files in this repository to a folder of your choice and run one of the scripts of the main folder. All codes generate/include the required data to run the simulations and optimization; simulations can take a few minutes on a standard laptop.

# License
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".

# References
1. AED Barioni, AN Montanari, AE Motter. Master stability with delays favors heterogeneous optimal networks. (2026)
2. K Engelborghs, T Luzyanina, D Roose. Numerical bifurcation analysis of delay differential equations using DDE-BIFTOOL. ACM Transactions on Mathematical Software, 28:1-21 (2002).
