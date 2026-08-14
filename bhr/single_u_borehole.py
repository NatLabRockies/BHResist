from math import log, pi

from scp.base_fluid import BaseFluid

from bhr.u_tube import UTube
from bhr.utilities import coth


class SingleUBorehole(UTube):
    """First-order multipole model for a symmetric, grouted single U-tube borehole."""

    def __init__(
        self,
        borehole_diameter: float,
        pipe_outer_diameter: float,
        pipe_dimension_ratio: float,
        length: float,
        shank_space: float,
        pipe_conductivity: float,
        grout_conductivity: float,
        soil_conductivity: float,
        fluid_type: str | None = None,
        fluid_concentration: float = 0,
        *,
        fluid: BaseFluid | None = None,
    ) -> None:
        """
        Constructs a grouted single U-tube borehole model.

        The resistance equations follow Javed and Spitler (2017), "Accuracy of
        Borehole Thermal Resistance Calculation Methods for Grouted Single U-tube
        Ground Heat Exchangers," Applied Energy 187:790-806.

        :param borehole_diameter: borehole diameter, m
        :param pipe_outer_diameter: outer diameter of one U-tube leg, m
        :param pipe_dimension_ratio: ratio of pipe outer diameter to wall thickness
        :param length: active borehole length, m
        :param shank_space: radial distance from the borehole center to a pipe center, m;
                            this is one-half of the center-to-center spacing ``s`` in the source
        :param pipe_conductivity: pipe thermal conductivity, W/(m-K)
        :param grout_conductivity: grout thermal conductivity, W/(m-K)
        :param soil_conductivity: ground thermal conductivity, W/(m-K)
        :param fluid_type: built-in fluid key accepted by scp.get_fluid;
                           omit when passing fluid
        :param fluid_concentration: antifreeze fraction from 0 to 0.6; ignored for water
        :param fluid: existing SecondaryCoolantProps fluid instance, including a
                      user-defined fluid created by scp.get_fluid
        """
        super().__init__(
            pipe_outer_diameter,
            pipe_dimension_ratio,
            length,
            shank_space,
            pipe_conductivity,
            fluid_type,
            fluid_concentration,
            fluid=fluid,
        )

        self.borehole_diameter = borehole_diameter
        self.grout_conductivity = grout_conductivity
        self.soil_conductivity = soil_conductivity
        # Javed and Spitler (2017), Eq. 14: dimensionless geometry and conductivity parameters.
        self.theta_1 = 2 * self.shank_space / self.borehole_diameter
        self.theta_2 = self.borehole_diameter / self.pipe_outer_diameter
        self.theta_3 = 1 / (2 * self.theta_1 * self.theta_2)
        self.sigma = (self.grout_conductivity - self.soil_conductivity) / (
            self.grout_conductivity + self.soil_conductivity
        )
        self.bh_length = length
        self.two_pi_kg = 2 * pi * self.grout_conductivity

        # Cached fluid-to-outer-pipe-wall resistance R_p for one U-tube leg.
        self.pipe_resist = None

    def update_beta(self, m_dot: float, temp: float) -> float:
        """
        Calculates and caches the dimensionless pipe-resistance parameter ``beta``.

        Javed and Spitler (2017), Equation 14, defines
        ``beta = 2 pi k_g R_p``. The equivalent definition is Equation 3.47 in
        Javed and Spitler (2016), in Rees (ed.).

        :param m_dot: total mass flow rate through the single U-tube, kg/s
        :param temp: mean fluid temperature, Celsius
        :return: dimensionless pipe-resistance parameter ``beta``
        """

        pipe_resist = self.calc_fluid_pipe_resist(m_dot, temp)
        self.pipe_resist = pipe_resist

        # Javed and Spitler (2017), Eq. 14: beta = 2 pi k_g R_p.
        beta = self.two_pi_kg * pipe_resist

        return beta

    def calc_direct_coupling_resistance(self, m_dot: float, temp: float) -> tuple[float, float]:
        """
        Calculates the direct fluid-to-fluid branch in the symmetric delta network.

        Javed and Spitler (2016), in Rees (ed.), Equations 3.12-3.14, relate
        ``R_1-b = 2 R_b`` and ``R_1-2 = 4 R_a R_b / (4 R_b - R_a)``.
        ``R_1-2`` is network-specific and may be negative even though the physical
        total internal resistance ``R_a`` remains positive.

        :param m_dot: total mass flow rate through the single U-tube, kg/s
        :param temp: mean fluid temperature, Celsius
        :return: direct coupling resistance ``R_1-2`` and local borehole resistance
                 ``R_b``, both in K/(W/m)
        """
        r_a = self.calc_total_internal_bh_resistance(m_dot, temp)
        r_b = self.calc_local_bh_resistance(m_dot, temp)

        # Javed and Spitler (2016), Eqs. 3.12-3.14: R_1-2 = 4 R_a R_b / (4 R_b - R_a).
        r_12 = (4 * r_a * r_b) / (4 * r_b - r_a)

        resist_bh_direct_coupling = r_12
        return resist_bh_direct_coupling, r_b

    def calc_local_bh_resistance(self, m_dot: float, temp: float) -> float:
        """
        Calculates local borehole resistance using the first-order multipole method.

        Javed and Spitler (2017), Equation 13, gives the resistance ``R_b``
        between the mean U-tube fluid temperature and the borehole wall.

        :param m_dot: total mass flow rate through the single U-tube, kg/s
        :param temp: mean fluid temperature, Celsius
        :return: local borehole resistance ``R_b``, K/(W/m)
        """
        beta = self.update_beta(m_dot, temp)

        final_term_1 = log(self.theta_2 / (2 * self.theta_1 * (1 - self.theta_1**4) ** self.sigma))

        term_2_num = self.theta_3**2 * (1 - (4 * self.sigma * self.theta_1**4) / (1 - self.theta_1**4)) ** 2
        term_2_den_pt_1 = (1 + beta) / (1 - beta)
        term_2_den_pt_2 = self.theta_3**2 * (1 + (16 * self.sigma * self.theta_1**4) / (1 - self.theta_1**4) ** 2)
        term_2_den = term_2_den_pt_1 + term_2_den_pt_2
        final_term_2 = term_2_num / term_2_den

        # Javed and Spitler (2017), Eq. 13: first-order local resistance R_b.
        resist_bh_ave = (1 / (4 * pi * self.grout_conductivity)) * (beta + final_term_1 - final_term_2)
        return resist_bh_ave

    def calc_total_internal_bh_resistance(self, m_dot: float, temp: float) -> float:
        """
        Calculates total internal resistance using the first-order multipole method.

        Javed and Spitler (2017), Equation 26, gives the resistance ``R_a``
        between the downward- and upward-flowing U-tube legs.

        :param m_dot: total mass flow rate through the single U-tube, kg/s
        :param temp: mean fluid temperature, Celsius
        :return: total internal resistance ``R_a``, K/(W/m)
        """
        beta = self.update_beta(m_dot, temp)

        term_1_num = (1 + self.theta_1**2) ** self.sigma
        term_1_den = self.theta_3 * (1 - self.theta_1**2) ** self.sigma
        final_term_1 = log(term_1_num / term_1_den)

        term_2_num = self.theta_3**2 * (1 - self.theta_1**4 + 4 * self.sigma * self.theta_1**2) ** 2
        term_2_den_pt_1 = (1 + beta) / (1 - beta) * (1 - self.theta_1**4) ** 2
        term_2_den_pt_2 = self.theta_3**2 * (1 - self.theta_1**4) ** 2
        term_2_den_pt_3 = 8 * self.sigma * self.theta_1**2 * self.theta_3**2 * (1 + self.theta_1**4)
        term_2_den = term_2_den_pt_1 - term_2_den_pt_2 + term_2_den_pt_3
        final_term_2 = term_2_num / term_2_den

        # Javed and Spitler (2017), Eq. 26: first-order total internal resistance R_a.
        resist_bh_total_internal = 1 / (pi * self.grout_conductivity) * (beta + final_term_1 - final_term_2)

        return resist_bh_total_internal

    def calc_grout_resistance(self, m_dot: float, temp: float) -> float:
        """
        Calculates grout resistance from local and pipe resistances.

        Javed and Spitler (2017), Equation 3, gives ``R_b = R_g + R_p/N``.
        A single U-tube has two pipe legs, so ``N = 2``.

        :param m_dot: total mass flow rate through the single U-tube, kg/s
        :param temp: mean fluid temperature, Celsius
        :return: grout resistance ``R_g``, K/(W/m)
        """

        if self.pipe_resist is None:
            raise ValueError("Pipe resistance has not been calculated yet.")

        # Javed and Spitler (2017), Eq. 3 with N = 2: R_g = R_b - R_p / 2.
        resist_bh_grout = self.calc_local_bh_resistance(m_dot, temp) - self.pipe_resist / 2.0
        return resist_bh_grout

    def calc_effective_bh_resistance_uhf(self, m_dot: float, temp: float) -> float:
        """
        Calculates effective resistance for a uniform heat-flux boundary condition.

        Javed and Spitler (2016), in Rees (ed.), Equation 3.67, accounts for
        axial fluid-temperature variation using local resistance ``R_b`` and total
        internal resistance ``R_a``.

        :param m_dot: total mass flow rate through the single U-tube, kg/s
        :param temp: mean fluid temperature, Celsius
        :return: effective borehole resistance ``R_b*`` under UHF, K/(W/m)
        """

        r_a = self.calc_total_internal_bh_resistance(m_dot, temp)
        r_b = self.calc_local_bh_resistance(m_dot, temp)

        pt_1 = 1 / (3 * r_a)
        pt_2 = (self.bh_length / (self.fluid.cp(temp) * m_dot)) ** 2
        resist_short_circuiting = pt_1 * pt_2

        # Javed and Spitler (2016), Eq. 3.67: R_b* = R_b + R_v^2 / (3 R_a).
        resist_bh_effective_uhf = r_b + resist_short_circuiting
        return resist_bh_effective_uhf

    def calc_effective_bh_resistance_ubwt(self, m_dot: float, temp: float) -> float:
        """
        Calculates effective resistance for a uniform borehole-wall temperature.

        Javed and Spitler (2016), in Rees (ed.), Equations 3.68-3.69, define
        ``R_b* = R_b eta coth(eta)``. The expression for ``eta`` below is the
        algebraically equivalent form obtained using Equations 3.12-3.14.

        :param m_dot: total mass flow rate through the single U-tube, kg/s
        :param temp: mean fluid temperature, Celsius
        :return: effective borehole resistance ``R_b*`` under UBWT, K/(W/m)
        """

        r_a = self.calc_total_internal_bh_resistance(m_dot, temp)
        r_b = self.calc_local_bh_resistance(m_dot, temp)
        # Javed and Spitler (2016), Eq. 3.69: R_v = H / (m_dot c_p).
        r_v = self.bh_length / (m_dot * self.fluid.cp(temp))

        # Javed and Spitler (2016), Eqs. 3.12-3.14 and 3.69: eta = R_v / sqrt(R_b R_a).
        eta = r_v / (r_b * r_a) ** 0.5

        # Javed and Spitler (2016), Eq. 3.68: R_b* = R_b eta coth(eta).
        resist_bh_effective_ubwt = r_b * eta * coth(eta)

        return resist_bh_effective_ubwt
