from typing import cast

from bhr.coaxial_borehole import Coaxial
from bhr.double_u_borehole import DoubleUTube
from bhr.enums import BoreholeType, BoundaryCondition
from bhr.single_u_borehole import SingleUBorehole
from bhr.utilities import set_boundary_condition_enum

AnyBHType = Coaxial | DoubleUTube | SingleUBorehole | None


class Borehole:
    def __init__(self):
        self._bh_type = None
        self._boundary_condition = None
        self._bh: AnyBHType = None
        self.length = None

    def init_single_u_borehole(
        self,
        borehole_diameter: float,
        pipe_outer_diameter: float,
        pipe_dimension_ratio: float,
        length: float,
        shank_space: float,
        pipe_conductivity: float,
        grout_conductivity: float,
        soil_conductivity: float,
        fluid_type: str,
        fluid_concentration: float = 0,
        boundary_condition: str = "UNIFORM_HEAT_FLUX",
    ) -> None:
        """
        Constructs a grouted single u-tube borehole.

        :param borehole_diameter: borehole diameter, in m.
        :param pipe_outer_diameter: outer diameter of the pipe, in m.
        :param pipe_dimension_ratio: non-dimensional ratio of pipe diameter to pipe thickness.
        :param length: length of borehole from top to bottom, in m.
        :param shank_space: radial distance from the borehole center to the pipe center, in m.
        :param pipe_conductivity: pipe thermal conductivity, in W/m-K.
        :param grout_conductivity: grout thermal conductivity, in W/m-K.
        :param soil_conductivity: soil thermal conductivity, in W/m-K.
        :param fluid_type: fluid type. "ETHYLALCOHOL", "ETHYLENEGLYCOL", "METHYLALCOHOL", "PROPYLENEGLYCOL", or "WATER"
        :param fluid_concentration: fractional concentration of antifreeze mixture, from 0-0.6.
        :param boundary_condition: borehole wall boundary condition. "UNIFORM_HEAT_FLUX" or "UNIFORM_BOREHOLE_WALL_TEMP"
        """

        self._bh_type = BoreholeType.SINGLE_U_TUBE
        self._boundary_condition = set_boundary_condition_enum(boundary_condition)
        self._bh = SingleUBorehole(
            borehole_diameter,
            pipe_outer_diameter,
            pipe_dimension_ratio,
            length,
            shank_space,
            pipe_conductivity,
            grout_conductivity,
            soil_conductivity,
            fluid_type,
            fluid_concentration,
        )

    def init_double_u_borehole(
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
        boundary_condition: str = "UNIFORM_HEAT_FLUX",
    ) -> None:
        """
        Constructs a grouted double u-tube borehole with u-tubes in parallel.

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
        :param boundary_condition: borehole wall boundary condition. "UNIFORM_HEAT_FLUX" or "UNIFORM_BOREHOLE_WALL_TEMP"
        """

        self._bh_type = BoreholeType.DOUBLE_U_TUBE
        self._boundary_condition = set_boundary_condition_enum(boundary_condition)
        self._bh = DoubleUTube(
            borehole_diameter,
            pipe_outer_diameter,
            pipe_dimension_ratio,
            length,
            shank_space,
            pipe_conductivity,
            pipe_inlet_arrangement,
            grout_conductivity,
            soil_conductivity,
            fluid_type,
            fluid_concentration,
        )

    def init_coaxial_borehole(
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
        boundary_condition: str = "UNIFORM_HEAT_FLUX",
    ) -> None:
        """
        Constructs a grouted coaxial borehole.

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
        :param boundary_condition: borehole wall boundary condition. "UNIFORM_HEAT_FLUX" or "UNIFORM_BOREHOLE_WALL_TEMP"
        """

        self._bh_type = BoreholeType.COAXIAL
        self._boundary_condition = set_boundary_condition_enum(boundary_condition)
        self._bh = Coaxial(
            borehole_diameter,
            outer_pipe_outer_diameter,
            outer_pipe_dimension_ratio,
            outer_pipe_conductivity,
            inner_pipe_outer_diameter,
            inner_pipe_dimension_ratio,
            inner_pipe_conductivity,
            length,
            grout_conductivity,
            soil_conductivity,
            fluid_type,
            fluid_concentration,
        )

    def init_from_dict(self, inputs: dict):
        """
        Constructs a borehole from a set of dictionary inputs.

        :param inputs: dict of input data.
        """

        bh_type_str = inputs["borehole_type"].upper()
        if bh_type_str == BoreholeType.SINGLE_U_TUBE.name:
            self._bh_type = BoreholeType.SINGLE_U_TUBE
        elif bh_type_str == BoreholeType.DOUBLE_U_TUBE.name:
            self._bh_type = BoreholeType.DOUBLE_U_TUBE
        elif bh_type_str == BoreholeType.COAXIAL.name:
            self._bh_type = BoreholeType.COAXIAL
        else:
            raise LookupError(f'borehole_type "{bh_type_str}" not supported')

        bc_str = inputs.get("boundary_condition", BoundaryCondition.UNIFORM_HEAT_FLUX.name).upper()
        if bc_str == BoundaryCondition.UNIFORM_HEAT_FLUX.name:
            self._boundary_condition = BoundaryCondition.UNIFORM_HEAT_FLUX
        elif bc_str == BoundaryCondition.UNIFORM_BOREHOLE_WALL_TEMP.name:
            self._boundary_condition = BoundaryCondition.UNIFORM_BOREHOLE_WALL_TEMP
        else:
            raise LookupError(f'boundary_condition "{bc_str}" not supported')

        bh_diameter = inputs["borehole_diameter"]
        length = inputs["length"]
        grout_conductivity = inputs["grout_conductivity"]
        soil_conductivity = inputs["soil_conductivity"]
        fluid_type = inputs["fluid_type"]
        fluid_concentration = inputs["fluid_concentration"]

        if self._bh_type == BoreholeType.SINGLE_U_TUBE:
            pipe_outer_dia_single = inputs["single_u_tube"]["pipe_outer_diameter"]
            dimension_ratio_single = inputs["single_u_tube"]["pipe_dimension_ratio"]
            shank_space_single = inputs["single_u_tube"]["shank_space"]
            pipe_conductivity_single = inputs["single_u_tube"]["pipe_conductivity"]

            self.init_single_u_borehole(
                bh_diameter,
                pipe_outer_dia_single,
                dimension_ratio_single,
                length,
                shank_space_single,
                pipe_conductivity_single,
                grout_conductivity,
                soil_conductivity,
                fluid_type,
                fluid_concentration,
                bc_str,
            )

        elif self._bh_type == BoreholeType.DOUBLE_U_TUBE:
            pipe_outer_dia_double = inputs["double_u_tube"]["pipe_outer_diameter"]
            dimension_ratio_double = inputs["double_u_tube"]["pipe_dimension_ratio"]
            shank_space_double = inputs["double_u_tube"]["shank_space"]
            pipe_conductivity_double = inputs["double_u_tube"]["pipe_conductivity"]
            pipe_inlet_arrangement = inputs["double_u_tube"]["pipe_inlet_arrangement"]

            self.init_double_u_borehole(
                bh_diameter,
                pipe_outer_dia_double,
                dimension_ratio_double,
                length,
                shank_space_double,
                pipe_conductivity_double,
                pipe_inlet_arrangement,
                grout_conductivity,
                soil_conductivity,
                fluid_type,
                fluid_concentration,
                bc_str,
            )

        elif self._bh_type == BoreholeType.COAXIAL:
            pipe_outer_dia_coax = inputs["coaxial"]["outer_pipe_outer_diameter"]
            outer_pipe_dimension_ratio = inputs["coaxial"]["outer_pipe_dimension_ratio"]
            pipe_conductivity_coax = inputs["coaxial"]["outer_pipe_conductivity"]
            inner_pipe_outer_diameter = inputs["coaxial"]["inner_pipe_outer_diameter"]
            inner_pipe_dimension_ratio = inputs["coaxial"]["inner_pipe_dimension_ratio"]
            inner_pipe_conductivity = inputs["coaxial"]["inner_pipe_conductivity"]

            self.init_coaxial_borehole(
                bh_diameter,
                pipe_outer_dia_coax,
                outer_pipe_dimension_ratio,
                pipe_conductivity_coax,
                inner_pipe_outer_diameter,
                inner_pipe_dimension_ratio,
                inner_pipe_conductivity,
                length,
                grout_conductivity,
                soil_conductivity,
                fluid_type,
                fluid_concentration,
                bc_str,
            )

        else:
            raise NotImplementedError(f'bh_type "{self._bh_type.name}" not implemented')

    def calc_bh_resist(self, mass_flow_rate: float, temperature: float) -> float:
        """
        Computes the effective borehole thermal resistance.

        :param mass_flow_rate: total borehole mass flow rate, in kg/s
        :param temperature: average fluid temperature, in Celsius
        :return: effective borehole resistance, in K/W-m
        """

        if self._bh is None:
            raise TypeError("Borehole not initialized")

        if self._boundary_condition == BoundaryCondition.UNIFORM_HEAT_FLUX:
            return self._bh.calc_effective_bh_resistance_uhf(mass_flow_rate, temperature)

        if self._boundary_condition == BoundaryCondition.UNIFORM_BOREHOLE_WALL_TEMP:
            return self._bh.calc_effective_bh_resistance_ubwt(mass_flow_rate, temperature)

        raise NotImplementedError(f'Boundary Condition: "{self._boundary_condition}" implemented.')

    def calc_pipe_cond_resist(self) -> float:
        """
        Computes the pipe conduction resistance.

        In the case of coaxial boreholes, the function only returns the conduction resistances for the outer pipe.

        :return: pipe conduction resistance, K/(W/m)
        """

        match self._bh_type:
            case BoreholeType.SINGLE_U_TUBE:
                return cast(SingleUBorehole, self._bh).calc_cond_resist()
            case BoreholeType.DOUBLE_U_TUBE:
                return cast(DoubleUTube, self._bh).calc_cond_resist()
            case BoreholeType.COAXIAL:
                return cast(Coaxial, self._bh).calc_cond_resist()[1]
            case _:
                raise NotImplementedError(f"{self._bh_type} not implemented.")

    def calc_fluid_resist(self, mass_flow_rate: float, temperature: float) -> float:
        """
        Computes the fluid convection resistance.

        ``mass_flow_rate`` is the total borehole flow. For a parallel double U-tube,
        it is divided equally between the two identical U-tubes, consistent with the
        flow definition used by Claesson and Javed (2019), Equations 42-46.

        In the case of coaxial boreholes, the function returns the sum of the convection resistances
        at the two surfaces of the annular region.

        :param mass_flow_rate: total borehole mass flow rate, kg/s
        :param temperature: average fluid temperature, Celsius
        :return: fluid convection, K/(W/m)
        """

        match self._bh_type:
            case BoreholeType.SINGLE_U_TUBE:
                return cast(SingleUBorehole, self._bh).calc_conv_resist(mass_flow_rate, temperature)
            case BoreholeType.DOUBLE_U_TUBE:
                # Claesson and Javed (2019), Eqs. 42-46: two identical U-tubes in parallel.
                return cast(DoubleUTube, self._bh).calc_conv_resist(mass_flow_rate / 2, temperature)
            case BoreholeType.COAXIAL:
                return sum(cast(Coaxial, self._bh).calc_conv_resist_annulus(mass_flow_rate, temperature))
            case _:
                raise NotImplementedError(f"{self._bh_type} not implemented.")

    def calc_fluid_pipe_resist(self, mass_flow_rate: float, temperature: float) -> float:
        """
        Computes the fluid convection + pipe conduction resistance.

        ``mass_flow_rate`` is the total borehole flow. For a parallel double U-tube,
        it is divided equally between the two identical U-tubes, consistent with the
        flow definition used by Claesson and Javed (2019), Equations 42-46.

        For coaxial boreholes, this returns ``R_c,ao + R_p,a`` from Grundmann
        (2016), Equation 4.5: convection at the inner surface of the outer pipe plus
        outer-pipe conduction. It excludes the grout term from Equation 4.5.

        :param mass_flow_rate: total borehole mass flow rate, kg/s
        :param temperature: average fluid temperature, Celsius
        :return: fluid convection + pipe conduction resistance, K/(W/m)
        """

        if self._bh is None:
            raise TypeError("Borehole not initialized")

        if self._bh_type is None:
            raise NotImplementedError(f"{self._bh_type} not implemented.")

        if self._bh_type == BoreholeType.DOUBLE_U_TUBE:
            # Claesson and Javed (2019), Eqs. 42-46: two identical U-tubes in parallel.
            return cast(DoubleUTube, self._bh).calc_fluid_pipe_resist(mass_flow_rate / 2, temperature)

        return self._bh.calc_fluid_pipe_resist(mass_flow_rate, temperature)
