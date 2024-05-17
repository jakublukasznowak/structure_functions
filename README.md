The repository contains MATLAB code developed for the purpose of the analysis presented in the manuscript

*Scale-by-scale budget of turbulence kinetic energy in the convective atmospheric boundary layer: analysis of structure functions*

by Jakub L. Nowak, Marta Wacławczyk, John C. Vassilicos, S. Król and Szymon P. Malinowski


## Abstract

We investigate scale-by-scale budget of turbulence kinetic energy in convective atmospheric boundary layers using airborne measurements performed in shallow trade-wind regime over ocean by the ATR aircraft during EUREC4A experiment. The simple and repeatable flight pattern sampled four altitude levels: near-surface, the middle and top of the subcloud layer, and the cloud base. This approach provides acceptable level of convergence for the statistics, including the mixed velocity-temperature structure function and the third order structure functions of velocity.
 
The derived scale-by-scale budget is approximately closed up to the length-scales of about 200 m at all four levels. At the near-surface and mid-subcloud levels, the buoyancy forcing supplies energy to the largest scales and becomes negligible at smaller scales. As a result, Kolmogorov equilibrium is observed over a wide range of scales, from a few up to hundred meters. On the other hand, at the top-subcloud and cloud-base levels, the budget at scales above 10 m is far from Kolmogorov equilibrium, with a significant contribution of buoyancy forcing and turbulent transport. 

The contribution to the buoyancy forcing related to humidity variations is substantial and strictly positive at all four levels, even at the cloud-base where it balances the negative contribution related to temperature only.


## Data

The measurements analysed in this study were performed within the field experiment EUREC4A (Elucidating the role of cloud–circulation coupling in climate) in Jan - Feb 2020 in trade-wind cumulus regime in northwestern Atlantic off the coast of Barbados. For the overview of EUREC4A, see [Stevens et al. (2021)](https://doi.org/10.5194/ESSD-13-4067-2021). The meteorological conditions and the structure of the boundary layer are analysed in detail by [Albright et al. (2022)](https://doi.org/10.1175/JAS-D-21-0337.1).

The turbulence measurements in the boundary layer were obtained with a the French research aircraft SAFIRE ATR42. Its instrumentation and sampling strategy are described by [Bony et al. (2022)](https://doi.org/10.5194/ESSD-14-2021-2022). The relevant datasets can be downloaded from the public repositories:
- [Lothon, M. & Brilouet, P. (2020). SAFIRE ATR42: Turbulence Data 25 Hz.](https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/SAFIRE-TURB/PROCESSED/) doi:10.25326/128,
- [Coutris, P. (2021). SAFIRE ATR42: PMA/Cloud composite dataset.](https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/PMA/PROCESSED/CloudComposite/) doi:10.25326/237.


## Code

The code was developed in MATLAB R2019b. The functionality in other versions of this environment was not tested.

The script `main.m` performs the analysis described in the main part of the manuscript.

One external package is used:
- [YAML 1.1 parser and emitter for MATLAB](https://www.mathworks.com/matlabcentral/fileexchange/106765-yaml) by Martin Koch
