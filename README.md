# BHResist

A pure python library for computing thermal resistance within single-u, double-u, and coaxial grouted borehole heat exchangers. For single and double u-tube configurations, the methods use the 1st-order closed-form multipole approximations, which typically produces results with less than 1% error when compared to the 10th-order multipole method. Coaxial borehole methods apply a simple 1D resistance network method.

This is intended to be a lightweight library that can be easily imported into any other Python tool, with no bulky dependencies.

## Documentation

Documentation for BHResist can be found at https://bhresist.readthedocs.io.

## Citation

Mitchell, Matt, Adams, Sonja, Lee, Edwin, and Swindler, Alexander. BHResist [SWR-25-57]. Computer Software. https://github.com/NatLabRockies/BHResist. USDOE Office of Energy Efficiency and Renewable Energy (EERE), Renewable Power Office. Geothermal Technologies Office. 04 Apr. 2025. Web. doi:10.11578/dc.20250421.3.

## References

Hellström, G. 1991. "Ground Heat Storage: Thermal Analyses of Duct Storage Systems." PhD dissertation. Department of Mathematical Physics, University of Lund, Sweden.

Grundmann, R.M. 2016 "Improved design methods for ground heat exchangers." Master’s thesis, Oklahoma State University.

Javed, S. and J.D. Spitler. 2016. "Calculation of borehole thermal resistance." In _Advances in Ground-Source Heat Pump Systems_. Ed. S.J. Rees. Woodhead Publishing. https://doi.org/10.1016/B978-0-08-100311-4.00003-0

Javed, S., and J.D. Spitler. 2017. "Accuracy of borehole thermal resistance calculation methods for grouted single u-tube ground heat exchangers." _Applied Energy,_ 187:790-806. https://doi.org/10.1016/j.apenergy.2016.11.079

Claesson, J., and S. Javed. 2019. "Explicit multipole formulas and thermal network models for calculating thermal resistances of double U-pipe borehole heat exchangers." _Science and Technology for the Built Environment,_ 25(8) pp. 980–992. https://doi.org/10.1080/23744731.2019.1620565
