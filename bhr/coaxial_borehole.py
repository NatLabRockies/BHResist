from math import log, pi

from bhr.fluid import get_fluid
from bhr.pipe import Pipe
from bhr.utilities import coth, smoothing_function


class Coaxial:
    def __init__(
        self,
        borehole_diameter: float,
        outer_pipe_outer_diameter: float,
        outer_pipe_dimension_ratio: float,
        outer_pipe_conductivity: float,
        inner_pipe_outer_diameter: float,
        inner_pipe_dimension_ratio: float,
        inner_pipe_conductivity: float,
        length: float,
        grout_conductivity: float,
        soil_conductivity: float,
        fluid_type: str,
        fluid_concentration: float,
    ):
        """
        Implementation for computing borehole thermal resistance for a grouted coaxial borehole.

        :param borehole_diameter: borehole diameter, in m.
        :param outer_pipe_outer_diameter: outer diameter of outer pipe, in m.
        :param outer_pipe_dimension_ratio: non-dimensional ratio of outer pipe diameter to thickness.
        :param outer_pipe_conductivity: outer pipe thermal conductivity, in W/m-K.
        :param inner_pipe_outer_diameter: inner diameter of outer pipe, in m.
        :param inner_pipe_dimension_ratio: non-dimensional ratio of inner pipe diameter to thickness.
        :param inner_pipe_conductivity: inner pipe thermal conductivity, in W/m-K.
        :param length: length of borehole from top to bottom, in m.
        :param grout_conductivity: grout thermal conductivity, in W/m-K.
        :param soil_conductivity: pipe thermal conductivity, in W/m-K.
        :param fluid_type: fluid type. "ETHYLALCOHOL", "ETHYLENEGLYCOL", "METHYLALCOHOL",  "PROPYLENEGLYCOL", or "WATER"
        :param fluid_concentration: fractional concentration of antifreeze mixture, from 0-0.6.
        """

        self.borehole_diameter = borehole_diameter
        self.grout_conductivity = grout_conductivity
        self.soil_conductivity = soil_conductivity
        self.fluid = get_fluid(fluid_type, fluid_concentration)
        self.length = length

        self.outer_pipe = Pipe(
            outer_pipe_outer_diameter,
            outer_pipe_dimension_ratio,
            length,
            outer_pipe_conductivity,
            fluid_type,
            fluid_concentration,
        )
        self.inner_pipe = Pipe(
            inner_pipe_outer_diameter,
            inner_pipe_dimension_ratio,
            length,
            inner_pipe_conductivity,
            fluid_type,
            fluid_concentration,
        )

        self.annular_hydraulic_diameter = self.outer_pipe.pipe_inner_diameter - self.inner_pipe.pipe_outer_diameter
        self.annular_wetted_perimeter = pi * (self.outer_pipe.pipe_inner_diameter + self.inner_pipe.pipe_outer_diameter)

    def re_annulus(self, m_dot, temp):
        """
        Reynolds number for annulus flow

        :param m_dot: mass flow rate, kg/s
        :param temp: temperature, C
        :return: Reynolds number
        """

        return 4 * m_dot / (self.fluid.mu(temp) * self.annular_wetted_perimeter)

    def laminar_nusselt_annulus(self):
        """
        Laminar Nusselt numbers for annulus flow

        Hellström, G. 1991. Ground Heat Storage: Thermal Analyses of Duct Storage Systems.
        Department of Mathematical Physics, University of Lund, Sweden. pp 67-71

        :return: nu_ii: Laminar Nusselt number for inner surface of annulus pipe
        :return: nu_oo: Laminar Nusselt number for outer annulus pipe surface
        """
        nu_ii = 3.66 + 1.2 * (self.inner_pipe.pipe_outer_diameter / self.outer_pipe.pipe_inner_diameter) ** -0.8
        nu_oo = 3.66 + 1.2 * (self.inner_pipe.pipe_outer_diameter / self.outer_pipe.pipe_inner_diameter) ** 0.5

        return nu_ii, nu_oo

    def turbulent_nusselt_annulus(self, re, temp):
        """
        Turbulent Nusselt numbers for annulus flow

        Grundmann, Rachel Marie. "Improved design methods for ground heat exchangers."
        Master's thesis, Oklahoma State University, 2016.

        Eqns 4.10 and 4.11 based on the Dittus-Boelter equation

        :param re: Reynolds number
        :param temp: temperature, C
        :return: nu_ii: Turbulent Nusselt number for inner surface of annulus pipe
        :return: nu_oo: Turbulent Nusselt number for outer annulus pipe surface
        """

        pr = self.fluid.prandtl(temp)

        nu_ii = 0.023 * re**0.8 * pr**0.35
        nu_oo = nu_ii

        return nu_ii, nu_oo

    def calc_conv_resist_annulus(self, m_dot, temp):
        """
        Grundmann, Rachel Marie. "Improved design methods for ground heat exchangers."
        Master's thesis, Oklahoma State University, 2016.

        Eqns 4.4 - 4.11

        :param m_dot: mass flow rate, kg/s
        :param temp: temperature, C
        :return: r_conv_outside_inner_pipe: convective resistances along the outer wall of the inner pipe, K/(W/m)
        :return: r_conv_inside_outer_pipe: convective resistance along the inside wall of the outer pipe, K/(W/m)
        """

        # limit determined from Hellström, G. 1991. Ground Heat Storage: Thermal Analyses of
        # Duct Storage Systems. Department of Mathematical Physics, University of Lund, Sweden.

        low_reynolds = 2300

        # limit based on Dittus-Boelter equation
        high_reynolds = 10000

        re = self.re_annulus(m_dot, temp)

        if re < low_reynolds:
            # use this Nusselt number when the flow is laminar
            nu_ii, nu_oo = self.laminar_nusselt_annulus()

        elif low_reynolds <= re < high_reynolds:
            # in between
            nu_ii_low, nu_oo_low = self.laminar_nusselt_annulus()
            nu_ii_high, nu_oo_high = self.turbulent_nusselt_annulus(high_reynolds, temp)
            nu_ii = smoothing_function(re, low_reynolds, high_reynolds, nu_ii_low, nu_ii_high)
            nu_oo = smoothing_function(re, low_reynolds, high_reynolds, nu_oo_low, nu_oo_high)

        else:
            # use this Nusselt number when the flow is fully turbulent
            nu_ii, nu_oo = self.turbulent_nusselt_annulus(re, temp)

        r_conv_outside_inner_pipe = self.annular_hydraulic_diameter / (
            nu_ii * self.fluid.k(temp) * self.inner_pipe.pipe_outer_diameter * pi
        )

        r_conv_inside_outer_pipe = self.annular_hydraulic_diameter / (
            nu_oo * self.fluid.k(temp) * self.outer_pipe.pipe_inner_diameter * pi
        )

        return r_conv_outside_inner_pipe, r_conv_inside_outer_pipe

    def calc_local_bh_resistance(self, m_dot, temp):
        """
        Calculates the local resistance branches for a concentric coaxial borehole.

        Grundmann, Rachel Marie. "Improved design methods for ground heat exchangers."
        Master's thesis, Oklahoma State University, 2016.

        Equation 4.4 defines the fluid-to-fluid short-circuit resistance ``R_12 = R_a``.
        Equation 4.5 defines the annular-fluid-to-borehole-wall resistance ``R_b``.

        :param m_dot: mass flow rate, kg/s
        :param temp: temperature, C
        :return: total_series_resist: sum of ``R_a`` and ``R_b``, K/(W/m)
        :return: r_internal_resist: local internal borehole resistance ``R_a``, K/(W/m)
        :return: r_borehole_resist: local borehole resistance ``R_b``, K/(W/m)

        """
        # resistances progressing from inside to outside
        r_conv_inner_pipe = self.inner_pipe.calc_conv_resist(m_dot, temp)
        r_cond_inner_pipe, r_cond_outer_pipe = self.calc_cond_resist()
        r_conv_outside_inner_pipe = self.calc_conv_resist_annulus(m_dot, temp)[0]
        r_conv_inside_outer_pipe = self.calc_conv_resist_annulus(m_dot, temp)[1]
        r_cond_grout = log(self.borehole_diameter / self.outer_pipe.pipe_outer_diameter) / (
            2 * pi * self.grout_conductivity
        )

        # Grundmann (2016), Eq. 4.4: R_12 = R_a = R_c,dti + R_p,dt + R_c,ai.
        r_internal_resist = sum([r_conv_inner_pipe, r_cond_inner_pipe, r_conv_outside_inner_pipe])

        # Grundmann (2016), Eq. 4.5: R_b = R_c,ao + R_p,a + R_g.
        r_borehole_resist = sum([r_conv_inside_outer_pipe, r_cond_outer_pipe, r_cond_grout])
        total_series_resist = r_internal_resist + r_borehole_resist

        return [total_series_resist, r_internal_resist, r_borehole_resist]

    def calc_effective_bh_resistance_uhf(self, m_dot, temp):
        """
        Grundmann, Rachel Marie. "Improved design methods for ground heat exchangers."
        Master's thesis, Oklahoma State University, 2016.

        Eqn 4.33

        :param m_dot: mass flow rate, kg/s
        :param temp: temperature, C
        :return: effective_bhr_uhf: effective borehole resistance for
                 uniform heat flux boundary condition, K/(W/m)
        """

        _, r_a, r_b = self.calc_local_bh_resistance(m_dot, temp)
        r_v = self.length / (m_dot * self.fluid.cp(temp))  # K/(W/m)
        # Grundmann (2016), Eq. 4.33: R_b* = R_b + R_v^2 / (3 R_a).
        effective_bhr_uhf = r_b + 1 / (3 * r_a) * r_v**2

        return effective_bhr_uhf

    def calc_effective_bh_resistance_ubwt(self, m_dot, temp):
        """
        Grundmann, Rachel Marie. "Improved design methods for ground heat exchangers."
        Master's thesis, Oklahoma State University, 2016.

        Eqns 4.28 & 4.29

        :param m_dot: mass flow rate, kg/s
        :param temp: temperature, C
        :return: effective_bhr_ubwt: effective borehole resistance for
                 uniform borehole wall temperature boundary condition, K/(W/m)
        """

        _, r_a, r_b = self.calc_local_bh_resistance(m_dot, temp)
        r_v = self.length / (m_dot * self.fluid.cp(temp))  # K/(W/m)
        # Grundmann (2016), Eq. 4.29: eta accounts for axial fluid-temperature variation.
        eta = r_v / (2 * r_b) * (1 + 4 * r_b / r_a) ** (1 / 2)

        # Grundmann (2016), Eq. 4.28: R_b* = R_b eta coth(eta).
        effective_bhr_ubwt = r_b * eta * coth(eta)

        return effective_bhr_ubwt

    def calc_cond_resist(self) -> tuple[float, float]:
        """
        Computes the pipe conduction resistance for the inner and outer pipes.
        :return: pipe conduction resistance, K/(W/m)
        """
        return self.inner_pipe.calc_cond_resist(), self.outer_pipe.calc_cond_resist()

    def calc_conv_resist(self, m_dot, temp) -> tuple[float, float]:
        """
        Computes the convection resistance for the inner pipe and annular space between inner and outer pipes.
        :param m_dot: mass flow rate, kg/s
        :param temp: temperature, C
        :return: convection resistance, K/(W/m)
        """
        return self.inner_pipe.calc_conv_resist(m_dot, temp), sum(self.calc_conv_resist_annulus(m_dot, temp))

    def calc_fluid_pipe_resist(self, m_dot, temp):
        """
        Calculates the convection resistance at the inner surface of the outer pipe
        and the conduction resistance of the outer pipe.

        This is ``R_c,ao + R_p,a``, the fluid-to-outer-pipe-wall portion of
        Grundmann (2016), Equation 4.5. It excludes grout resistance and the
        inner-annulus convection resistance assigned to ``R_a`` by Equation 4.4.

        :param m_dot: mass flow rate, kg/s
        :param temp: temperature, C
        :return: outer annular convection resistance and outer pipe conduction resistance, K/(W/m)
        """

        _, r_cond_outer_pipe = self.calc_cond_resist()
        r_conv_inside_outer_pipe = self.calc_conv_resist_annulus(m_dot, temp)[1]

        # Grundmann (2016), Eq. 4.5, excluding the grout term R_g.
        return r_cond_outer_pipe + r_conv_inside_outer_pipe
