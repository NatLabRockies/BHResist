# BHResist

A pure Python library for calculating the effective thermal resistance of grouted single U-tube, parallel double U-tube, and coaxial borehole heat exchangers. The single- and double-U-tube models use first-order closed-form multipole approximations. The coaxial model uses a one-dimensional resistance network. Uniform heat-flux and uniform borehole-wall-temperature boundary conditions are supported.

The public calculation methods accept total borehole mass flow rate. For a parallel double U-tube, BHResist divides that flow equally between the two U-tubes. Thermal resistances are reported in K/(W/m), equivalently m-K/W.

BHResist is intended to be a lightweight library that can be imported into other Python applications without bulky dependencies.

## Documentation

Documentation for BHResist can be found at https://bhresist.readthedocs.io.

## Citation

Mitchell, Matt, Adams, Sonja, Lee, Edwin, and Swindler, Alexander. BHResist [SWR-25-57]. Computer Software. https://github.com/NatLabRockies/BHResist. USDOE Office of Energy Efficiency and Renewable Energy (EERE), Renewable Power Office. Geothermal Technologies Office. 04 Apr. 2025. Web. doi:10.11578/dc.20250421.3.

## References

Hellström, G. 1991. "Ground Heat Storage: Thermal Analyses of Duct Storage Systems." PhD dissertation. Department of Mathematical Physics, University of Lund, Sweden.

Grundmann, R.M. 2016. "Improved design methods for ground heat exchangers." Master’s thesis, Oklahoma State University.

Javed, S. and J.D. Spitler. 2016. "Calculation of borehole thermal resistance." In _Advances in Ground-Source Heat Pump Systems_. Ed. S.J. Rees. Woodhead Publishing. https://doi.org/10.1016/B978-0-08-100311-4.00003-0

Javed, S., and J.D. Spitler. 2017. "Accuracy of borehole thermal resistance calculation methods for grouted single U-tube ground heat exchangers." _Applied Energy,_ 187:790-806. https://doi.org/10.1016/j.apenergy.2016.11.079

Claesson, J., and S. Javed. 2019. "Explicit multipole formulas and thermal network models for calculating thermal resistances of double U-pipe borehole heat exchangers." _Science and Technology for the Built Environment,_ 25(8) pp. 980–992. https://doi.org/10.1080/23744731.2019.1620565
