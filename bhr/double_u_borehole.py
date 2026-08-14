from math import log as ln
from math import pi, sqrt

from bhr.enums import DoubleUPipeInletArrangement
from bhr.u_tube import UTube
from bhr.utilities import coth


class DoubleUTube(UTube):
    def __init__(
        self,
        borehole_diameter: float,
        pipe_outer_diameter: float,
        pipe_dimension_ratio: float,
        length: float,
        shank_space: float,
        pipe_conductivity: float,
        pipe_inlet_arrangement: str,
        grout_conductivity: float,
        soil_conductivity: float,
        fluid_type: str,
        fluid_concentration: float = 0,
    ):
        super().__init__(
            pipe_outer_diameter,
            pipe_dimension_ratio,
            length,
            shank_space,
            pipe_conductivity,
            fluid_type,
            fluid_concentration,
        )

        """
        Implementation for computing borehole thermal resistance for grouted boreholes with parallel double u-tubes.

        Relies primarily on the following references:

        Javed, S. & Spitler, J.D. 2017. 'Accuracy of Borehole Thermal Resistance Calculation Methods
        for Grouted Single U-tube Ground Heat Exchangers.' Applied Energy.187:790-806.

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas and Thermal Network Models
        for Calculating Thermal Resistances of Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980-92. doi:10.1080/23744731.2019.1620565.

        :param borehole_diameter: borehole diameter, in m.
        :param pipe_outer_diameter: outer diameter of the pipe, in m.
        :param pipe_dimension_ratio: non-dimensional ratio of pipe diameter to pipe thickness.
        :param length: length of borehole from top to bottom, in m.
        :param shank_space: radial distance from the borehole center to the pipe center, in m.
        :param pipe_conductivity: pipe thermal conductivity, in W/m-K.
        :param pipe_inlet_arrangement: arrangement of the pipe inlets. "ADJACENT", or "DIAGONAL"
        :param grout_conductivity: grout thermal conductivity, in W/m-K.
        :param soil_conductivity: soil thermal conductivity, in W/m-K.
        :param fluid_type: fluid type. "ETHYLALCOHOL", "ETHYLENEGLYCOL", "METHYLALCOHOL",  "PROPYLENEGLYCOL", or "WATER"
        :param fluid_concentration: fractional concentration of antifreeze mixture, from 0-0.6.
        """

        # static parameters
        self.borehole_diameter = borehole_diameter
        self.grout_conductivity = grout_conductivity
        self.borehole_radius = borehole_diameter / 2  # radius of borehole (m)
        self.pipe_radius = pipe_outer_diameter / 2  # pipe outer radius (m)

        if pipe_inlet_arrangement == DoubleUPipeInletArrangement.ADJACENT.name:
            self.pipe_inlet_arrangement = DoubleUPipeInletArrangement.ADJACENT
        elif pipe_inlet_arrangement == DoubleUPipeInletArrangement.DIAGONAL.name:
            self.pipe_inlet_arrangement = DoubleUPipeInletArrangement.DIAGONAL
        else:
            msg = (
                f"Invalid pipe_inlet_arrangement. Use one of the allowed values: "
                f"{', '.join(map(str, DoubleUPipeInletArrangement._member_names_))}"
            )
            raise AssertionError(msg)

        self.bh_length = length  # length of borehole (m)
        self.grout_conductivity = grout_conductivity  # W/(m-K)
        self.soil_conductivity = soil_conductivity  # W/(m-K)

        # (m) radial distance between centers of symmetrically placed pipes and borehole center (rc)
        self.pipe_centers_radius = shank_space

        # thermal conductivity ratio, dimensionless
        self.sigma = (self.grout_conductivity - self.soil_conductivity) / (
            self.grout_conductivity + self.soil_conductivity
        )

        # Check if shank spacing realistic
        lower_shank_space_limit = sqrt(pipe_outer_diameter**2 / 2)
        upper_shank_space_limit = 0.5 * (borehole_diameter - pipe_outer_diameter)
        if shank_space < lower_shank_space_limit:
            msg = (
                "Shank spacing is too small and must be greater than the 2 pipe radii to prevent "
                "pipes from overlapping."
            )
            raise AssertionError(msg)
        elif shank_space > upper_shank_space_limit:
            msg = (
                "Shank spacing is too large and must be less than the borehole radius minus the pipe radius "
                "to prevent pipes from extending beyond the borehole wall."
            )
            raise AssertionError(msg)

        # static parameters - calc_bh_resist_local
        self.p_pc = self.pipe_radius**2 / (4 * self.pipe_centers_radius**2)
        self.p_c = self.pipe_centers_radius**2 / (self.borehole_radius**8 - self.pipe_centers_radius**8) ** 0.25
        self.p_b = self.borehole_radius**2 / (self.borehole_radius**8 - self.pipe_centers_radius**8) ** 0.25
        self.eight_pi_kg = 8 * pi * self.grout_conductivity
        self.b_2 = ln(self.borehole_radius**4 / (4 * self.pipe_radius * self.pipe_centers_radius**3))
        self.b_3 = ln(self.borehole_radius**8 / (self.borehole_radius**8 - self.pipe_centers_radius**8))

        # static parameter - calc_internal_resist
        self.two_pi_kg = 2 * pi * self.grout_conductivity
        self.c_1 = self.pipe_centers_radius / self.pipe_radius
        c_2 = self.borehole_radius**4 + self.pipe_centers_radius**4
        c_3 = self.borehole_radius**4 - self.pipe_centers_radius**4
        self.ln_c2_c3 = ln(c_2 / c_3)
        self.c_4 = self.p_c**2 * self.p_b**2
        self.c_5 = self.p_c**2 * self.p_b**6 + self.p_c**6 * self.p_b**2
        d_2 = self.borehole_radius**2 + self.pipe_centers_radius**2
        d_3 = self.borehole_radius**2 - self.pipe_centers_radius**2
        self.ln_d2_d3 = ln(d_2 / d_3)
        self.d_4 = 3 * self.p_c**3 * self.p_b**5 + self.p_c**7 * self.p_b
        self.d_5 = self.p_c * self.p_b**7 + 3 * self.p_c**5 * self.p_b**3

        # non-static parameters
        self.pipe_resist: float | None = None

    def update_b1(self, m_dot_per_u_tube: float, temperature: float) -> float:
        """
        Updates b1 coefficient.

        Javed, S. & Spitler, J.D. Calculation of Borehole Thermal Resistance. In 'Advances in
        Ground-Source Heat Pump Systems,' pp. 84. Rees, S.J. ed. Cambridge, MA. Elsevier Ltd. 2016.

        Eq: 3-47

        Javed, S. & Spitler, J.D. 2017. 'Accuracy of Borehole Thermal Resistance Calculation Methods
        for Grouted Single U-tube Ground Heat Exchangers.' Applied Energy.187:790-806.

        Eq: 14

        :param m_dot_per_u_tube: mass flow rate in each u-tube, kg/s
        :param temperature: temperature, Celsius

        :return: b1: a ratio of (1-beta)/(1+beta), dependent on pipe resistance
                    & grout conductivity, dimensionless
        """

        pipe_resist = self.calc_fluid_pipe_resist(m_dot_per_u_tube, temperature)
        self.pipe_resist = pipe_resist
        beta = 2 * pi * self.grout_conductivity * pipe_resist
        b1 = (1 - beta) / (1 + beta)  # dimensionless parameter

        return b1

    def calc_bh_resist_local(self, m_dot_per_u_tube: float, temperature: float) -> float:
        """
        Calculates tube-to-borehole resistance (aka local borehole resistance).

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980-92. doi:10.1080/23744731.2019.1620565.

        Eq: 13 & 14

        :param m_dot_per_u_tube: mass flow rate in each u-tube, kg/s
        :param temperature: temperature, Celsius

        :return: borehole_resist_local: local borehole resistance, K/(W/m)
        """

        b1 = self.update_b1(m_dot_per_u_tube, temperature)

        if self.pipe_resist is None:
            raise ValueError("Pipe resistance has not been calculated yet.")

        # --Borehole resistance, 0th order [K/(W/m)]--
        rb0 = self.pipe_resist / 4 + 1 / self.eight_pi_kg * (self.b_2 + self.sigma * self.b_3)

        # --Borehole resistance, 1st order [K/(W/m)]--
        borehole_resist_local = rb0 - 1 / self.eight_pi_kg * (
            b1 * self.p_pc * (3 - 8 * self.sigma * self.p_c**4) ** 2
        ) / (1 + b1 * self.p_pc * (5 + 64 * self.sigma * self.p_c**4 * self.p_b**4))

        return borehole_resist_local

    def calc_internal_resist(self, m_dot_per_u_tube: float, temperature: float) -> float:
        """
        Calculates tube-to-tube resistance (aka internal resistance).

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980-92. doi:10.1080/23744731.2019.1620565.

        Eq: 18, 19, 22, 23

        :param m_dot_per_u_tube: mass flow rate in each u-tube, kg/s
        :param temperature: temperature, Celsius

        :return: internal_resist: local internal resistance, K/(W/m)
        """

        b1 = self.update_b1(m_dot_per_u_tube, temperature)

        if self.pipe_resist is None:
            raise ValueError("Pipe resistance has not been calculated yet.")

        if self.pipe_inlet_arrangement == DoubleUPipeInletArrangement.DIAGONAL:
            # 0th order
            ra0 = 2 * self.pipe_resist + 2 / self.two_pi_kg * (ln(self.c_1) + self.sigma * self.ln_c2_c3)

            # 1st order
            internal_resist = ra0 - 2 / self.two_pi_kg * (b1 * self.p_pc * (1 + 8 * self.sigma * self.c_4) ** 2) / (
                1 - b1 * self.p_pc * (3 - 32 * self.sigma * self.c_5)
            )

            return internal_resist

        elif self.pipe_inlet_arrangement == DoubleUPipeInletArrangement.ADJACENT:
            # 0th order
            ra0 = 2 * self.pipe_resist + 2 / self.two_pi_kg * (ln(2 * self.c_1) + self.sigma * self.ln_d2_d3)

            # 1st order
            matrix_element_11 = 1 + 16 * b1 * self.sigma * self.p_pc * self.d_4
            matrix_element_22 = -1 - 16 * b1 * self.sigma * self.p_pc * self.d_5
            matrix_element_21 = b1 * self.p_pc
            vector_1 = 1 - 8 * self.sigma * self.p_c**3 * self.p_b
            vector_2 = 3 + 8 * self.sigma * self.p_c * self.p_b**3

            internal_resist = ra0 + 2 / self.two_pi_kg * b1 * self.p_pc / 2 * (
                vector_2**2 * matrix_element_11
                - 2 * vector_1 * vector_2 * matrix_element_21
                - vector_1**2 * matrix_element_22
            ) / (matrix_element_11 * matrix_element_22 + matrix_element_21**2)

            return internal_resist

        raise AssertionError("Developer error. Invalid pipe inlet arrangement.")

    def calc_effective_bh_resistance_uhf(self, m_dot: float, temp: float) -> float:
        """
        Calculates effective borehole resistance for uniform heat flux along the borehole.

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980-92. doi:10.1080/23744731.2019.1620565.

        Equation 44. The paper defines ``V_f`` as flow in one of the two
        parallel U-tubes; ``m_dot`` is total borehole flow.

        :param m_dot: total borehole mass flow rate, kg/s
        :param temp: temperature, Celsius

        :return: effective_bhr_uhf: effective borehole resistance under uniform heat flux boundary conditions [K/(W/m)]
        """

        # Claesson and Javed (2019), Eqs. 42-44: the symmetric double U-tube is
        # two identical U-tubes in parallel; V_f is the flow in one U-tube.
        m_dot_per_u_tube = m_dot / 2
        internal_resist = self.calc_internal_resist(m_dot_per_u_tube, temp)
        borehole_resist_local = self.calc_bh_resist_local(m_dot_per_u_tube, temp)
        rv = self.bh_length / (self.fluid.cp(temp) * m_dot_per_u_tube)  # (K/(w/m)) thermal resistance factor

        # Claesson and Javed (2019), Eq. 44: R_b* = R_b + R_v^2 / (6 R_a).
        effective_bhr_uhf = borehole_resist_local + rv**2 / (6 * internal_resist)

        return effective_bhr_uhf

    def calc_effective_bh_resistance_ubwt(self, m_dot: float, temp: float) -> float:
        """
        Calculates effective borehole resistance for uniform borehole wall temperature.

        Claesson, Johan, and Saqib Javed. 2019. “Explicit Multipole Formulas
        and Thermal Network Models for Calculating Thermal Resistances of
        Double U-Pipe Borehole Heat Exchangers.” Science and Technology for
        the Built Environment 25 (8): 980-92. doi:10.1080/23744731.2019.1620565.

        Equation 46. The paper defines ``V_f`` as flow in one of the two
        parallel U-tubes; ``m_dot`` is total borehole flow.

        :param m_dot: total borehole mass flow rate, kg/s
        :param temp: temperature, Celsius

        :return: effective_bhr_ubwt: effective borehole resistance for uniform borehole wall temperature
                                    boundary condition [K/(W/m)]
        """

        # Claesson and Javed (2019), Eqs. 42 and 46: the symmetric double U-tube
        # is two identical U-tubes in parallel; V_f is the flow in one U-tube.
        m_dot_per_u_tube = m_dot / 2
        internal_resist = self.calc_internal_resist(m_dot_per_u_tube, temp)
        borehole_resist_local = self.calc_bh_resist_local(m_dot_per_u_tube, temp)

        rv = self.bh_length / (self.fluid.cp(temp) * m_dot_per_u_tube)  # (K/(w/m)) thermal resistance factor
        # Claesson and Javed (2019), Eq. 46: eta = R_v / sqrt(2 R_b R_a).
        n = rv / (2 * borehole_resist_local * internal_resist) ** 0.5

        # Claesson and Javed (2019), Eq. 46: R_b* = R_b eta coth(eta).
        effective_bhr_ubwt = borehole_resist_local * n * coth(n)

        return effective_bhr_ubwt
