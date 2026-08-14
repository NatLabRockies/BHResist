from math import log, pi

from scp.base_fluid import BaseFluid

from bhr.fluid import resolve_fluid
from bhr.utilities import inch_to_m, smoothing_function


class Pipe:
    def __init__(
        self,
        pipe_outer_diameter: float,
        pipe_dimension_ratio: float,
        pipe_length: float,
        pipe_conductivity: float,
        fluid_type: str | None = None,
        fluid_concentration: float = 0,
        *,
        fluid: BaseFluid | None = None,
    ):
        self.fluid = resolve_fluid(fluid_type, fluid_concentration, fluid=fluid)

        # ratio of outer diameter to wall thickness
        self.dimension_ratio = pipe_dimension_ratio

        # set diameters and thickness
        self.pipe_outer_diameter = pipe_outer_diameter
        self.pipe_inner_diameter = self.pipe_outer_diameter * (1 - 2 / self.dimension_ratio)
        self.thickness = self.pipe_outer_diameter / self.dimension_ratio

        # set length
        self.pipe_length = pipe_length

        # set physical properties
        self.pipe_conductivity = pipe_conductivity

        # compute cross-sectional areas
        self.area_cr_inner = pi / 4 * self.pipe_inner_diameter**2
        self.area_cr_outer = pi / 4 * self.pipe_outer_diameter**2
        self.area_cr_pipe = self.area_cr_outer - self.area_cr_inner

        # compute surface areas
        self.area_s_inner = pi * self.pipe_inner_diameter * self.pipe_length
        self.area_s_outer = pi * self.pipe_outer_diameter * self.pipe_length

        # compute volumes
        self.total_vol = self.area_cr_outer * self.pipe_length
        self.fluid_vol = self.area_cr_inner * self.pipe_length
        self.pipe_wall_vol = self.area_cr_pipe * self.pipe_length

    def set_fluid(self, fluid: BaseFluid) -> None:
        """
        Adopt an existing SecondaryCoolantProps fluid instance.

        This supports user-defined fluids created by scp.get_fluid.

        :param fluid: SecondaryCoolantProps fluid instance
        """
        self.fluid = resolve_fluid(fluid=fluid)

    @staticmethod
    def get_inner_dia(outer_dia: float, dimension_ratio: float) -> float:
        return outer_dia * (1 - 2 / dimension_ratio)

    def get_pipe_diameters_imperial(self, nominal_pipe_size_inches: float, dimension_ratio: float):
        if nominal_pipe_size_inches == 0.75:
            outer_dia = 1.05
        elif nominal_pipe_size_inches == 1.0:
            outer_dia = 1.315
        elif nominal_pipe_size_inches == 1.25:
            outer_dia = 1.66
        elif nominal_pipe_size_inches == 1.5:
            outer_dia = 1.9
        elif nominal_pipe_size_inches == 2.0:
            outer_dia = 2.375
        elif nominal_pipe_size_inches == 3.0:
            outer_dia = 3.5
        elif nominal_pipe_size_inches == 4.0:
            outer_dia = 4.5
        elif nominal_pipe_size_inches == 6.0:
            outer_dia = 6.625
        elif nominal_pipe_size_inches == 8.0:
            outer_dia = 8.625
        else:
            raise ValueError("Unsupported pipe size")

        return inch_to_m(self.get_inner_dia(outer_dia, dimension_ratio)), inch_to_m(outer_dia)

    def mdot_to_vdot(self, m_dot: float, temp: float) -> float:
        """
        Computes volumetric flow rate based on mass flow rate.

        :param m_dot: mass flow rate, in kg/s
        :param temp: temperature, in C
        :return: volumetric flow rate, in m3/s
        """

        return m_dot / self.fluid.density(temp)

    def mdot_to_re(self, m_dot: float, temp: float) -> float:
        """
        Computes Reynolds number based on mass flow rate.

        :param m_dot: mass flow rate, in kg/s
        :param temp: temperature, in C
        :return: Reynolds number, dimensionless
        """

        return 4 * m_dot / (self.fluid.mu(temp) * pi * self.pipe_inner_diameter)

    def mdot_to_velocity(self, m_dot: float, temp: float) -> float:
        """
        Computes velocity based on mass flow rate.

        :param m_dot: mass flow rate in, kg/s
        :param temp: temperature, in C
        :return: velocity, in m/s
        """

        return m_dot / (self.area_cr_inner * self.fluid.density(temp))

    def friction_factor(self, re: float) -> float:
        """
        Calculates the friction factor in smooth tubes

        Petukhov, B.S. 1970. 'Heat transfer and friction in turbulent pipe flow with variable physical properties.'
        In Advances in Heat Transfer, ed. T.F. Irvine and J.P. Hartnett, Vol. 6. New York Academic Press.

        :param re: Reynolds number, dimensionless
        """

        # limits picked be within about 1% of actual values
        low_reynolds = 2000
        high_reynolds = 4000

        if re < low_reynolds:
            return self.laminar_friction_factor(re)
        if re > high_reynolds:
            return self.turbulent_friction_factor(re)

        # pure laminar flow
        f_low = self.laminar_friction_factor(re)

        # pure turbulent flow
        f_high = self.turbulent_friction_factor(re)

        return smoothing_function(re, low_reynolds, high_reynolds, f_low, f_high)

    @staticmethod
    def laminar_friction_factor(re: float):
        """
        Laminar friction factor

        :param re: Reynolds number
        :return: friction factor
        """

        return 64.0 / re

    @staticmethod
    def turbulent_friction_factor(re: float):
        """
        Turbulent friction factor

        Petukhov, B. S. (1970). Advances in Heat Transfer, volume 6, Heat transfer and
        friction in turbulent pipe flow with variable physical properties, pages 503-564.
        Academic Press, Inc., New York, NY.

        :param re: Reynolds number
        :return: friction factor
        """

        return (0.79 * log(re) - 1.64) ** (-2.0)

    def pressure_loss(self, m_dot: float, temp: float) -> float:
        """
        Pressure loss in straight pipe

        :param m_dot: mass flow rate, kg/s
        :param temp: temperature, C
        :return: pressure loss, Pa
        """

        if m_dot <= 0:
            return 0

        re = self.mdot_to_re(m_dot, temp)
        term_1 = self.friction_factor(re) * self.pipe_length / self.pipe_inner_diameter
        term_2 = (self.fluid.density(temp) * self.mdot_to_velocity(m_dot, temp) ** 2) / 2

        return term_1 * term_2

    def pressure_loss_v_dot(self, v_dot: float, temp: float) -> float:
        """
        Pressure loss in a straight pipe

        :param v_dot: volume flow rate, m3/s
        :param temp: temperature, C
        :return: pressure loss, Pa
        """

        m_dot = self.fluid.density(temp) * v_dot
        return self.pressure_loss(m_dot, temp)

    @staticmethod
    def laminar_nusselt():
        """
        Laminar Nusselt number for smooth pipes

        mean(4.36, 3.66)
        :return: Nusselt number
        """
        return 4.01

    def turbulent_nusselt(self, re: float, temp: float):
        """
        Turbulent Nusselt number for smooth pipes

        Gnielinski, V. 1976. 'New equations for heat and mass transfer in turbulent pipe and channel flow.'
        International Chemical Engineering 16(1976), pp. 359-368.

        :param re: Reynolds number
        :param temp: temperature, C
        :return: Nusselt number
        """

        f = self.friction_factor(re)
        pr = self.fluid.prandtl(temp)
        return (f / 8) * (re - 1000) * pr / (1 + 12.7 * (f / 8) ** 0.5 * (pr ** (2 / 3) - 1))

    def calc_cond_resist(self) -> float:
        """
        Calculates the pipe radial conduction thermal resistance, in [K/(W/m)].

        Javed, S. and Spitler, J.D. 2017. 'Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers.' Applied Energy. 187: 790-806.

        :return: conduction resistance, K/(W/m)
        """

        return log(self.pipe_outer_diameter / self.pipe_inner_diameter) / (2 * pi * self.pipe_conductivity)

    def calc_conv_resist(self, m_dot: float, temp: float) -> float:
        """
        Calculates the pipe internal convection thermal resistance, in [k/(W/m)]

        Gnielinski, V. 1976. 'New equations for heat and mass transfer in turbulent pipe and channel flow.'
        International Chemical Engineering 16(1976), pp. 359-368.

        :param m_dot: mass flow rate, kg/s
        :param temp: temperature, C
        :return: convection resistance, K/(W/m)
        """

        low_reynolds = 2000
        high_reynolds = 4000

        re = self.mdot_to_re(m_dot, temp)

        if re < low_reynolds:
            nu = self.laminar_nusselt()
        elif low_reynolds <= re < high_reynolds:
            nu_low = self.laminar_nusselt()
            nu_high = self.turbulent_nusselt(high_reynolds, temp)
            nu = smoothing_function(re, low_reynolds, high_reynolds, nu_low, nu_high)
        else:
            nu = self.turbulent_nusselt(re, temp)

        return 1 / (nu * pi * self.fluid.k(temp))

    def calc_fluid_pipe_resist(self, m_dot: float, temp: float):
        """
        Calculates the combined convection and conduction pipe resistance

        Javed, S. and Spitler, J.D. 2017. 'Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers.' Applied Energy. 187: 790-806.

        Equation 3

        :param m_dot: mass flow rate, kg/s
        :param temp: temperature, C
        :return: combined convection and conduction pipe resistance, K/(W/m)
        """

        return self.calc_conv_resist(m_dot, temp) + self.calc_cond_resist()
