Preferred Programmatic Usage
============================

Use :class:`bhr.borehole.Borehole` as the public interface. Initialize one
supported borehole configuration, then call
:meth:`bhr.borehole.Borehole.calc_bh_resist` with the total borehole mass
flow rate and mean fluid temperature. The result is the effective borehole
resistance, :math:`R_b^*`, in K/(W/m).

The default boundary condition is ``"UNIFORM_HEAT_FLUX"``. Set
``boundary_condition="UNIFORM_BOREHOLE_WALL_TEMP"`` during initialization
to use the uniform borehole-wall-temperature formulation.

For both U-tube configurations, ``shank_space`` is the radial distance from
the borehole center to a pipe center. For a symmetric single U-tube, it is half
the center-to-center distance between the two pipe legs. For a parallel double
U-tube, the public mass-flow argument is the combined flow through both
U-tubes; BHResist assigns half to each U-tube.

Fluid selection
---------------

BHResist uses the public API from SecondaryCoolantProps 1.5 or newer. Select a
built-in fluid with the ``fluid_type`` and ``fluid_concentration`` initializer
arguments. The canonical keys are ``"water"``, ``"ethyl_alcohol"``,
``"ethylene_glycol"``, ``"methyl_alcohol"``, and ``"propylene_glycol"``.
Compact uppercase names accepted by earlier BHResist versions remain supported.

SecondaryCoolantProps can also construct a user-defined fluid. Viscosity,
specific heat, density, and conductivity are required. Each may be a constant
or a callable that accepts temperature in Celsius::

    from scp import get_fluid

    custom_fluid = get_fluid(
        "user_defined",
        name="BoreholeFluid",
        viscosity=lambda temp: 0.003 - 1.0e-5 * temp,  # Pa-s
        specific_heat=3200.0,  # J/(kg-K)
        density=lambda temp: 1050.0 - 0.4 * temp,  # kg/m3
        conductivity=0.42,  # W/(m-K)
        freeze_point=-12.0,
        t_min=-20.0,
        t_max=80.0,
    )

Pass the resulting object through the ``fluid`` argument and omit
``fluid_type`` and ``fluid_concentration``::

    custom_bhr = Borehole()
    custom_bhr.init_single_u_borehole(
        borehole_diameter=0.127,
        pipe_outer_diameter=0.032,
        pipe_dimension_ratio=11,
        length=200,
        shank_space=0.032,
        pipe_conductivity=0.4,
        grout_conductivity=1.6,
        soil_conductivity=2.0,
        fluid=custom_fluid,
    )

An initialized model can adopt another SecondaryCoolantProps fluid without
rebuilding its geometry::

    replacement_fluid = get_fluid("water")
    custom_bhr.set_fluid(replacement_fluid)

For a coaxial model, ``set_fluid`` updates the model and both internal pipe
objects to use the same fluid instance.

Single U-tube
-------------

::

    from bhr.borehole import Borehole

    single_bhr = Borehole()
    single_bhr.init_single_u_borehole(
        borehole_diameter=0.127,
        pipe_outer_diameter=0.032,
        pipe_dimension_ratio=11,
        length=200,
        shank_space=0.032,
        pipe_conductivity=0.4,
        grout_conductivity=1.6,
        soil_conductivity=2.0,
        fluid_type="propylene_glycol",
        fluid_concentration=0.2,
    )

    mass_flow_rate = 0.5  # kg/s
    temperature = 20.0  # Celsius
    print(f"{single_bhr.calc_bh_resist(mass_flow_rate, temperature):0.5f}")

.. image:: images/single-u.webp
   :width: 600
   :alt: Cross-section of a single U-tube borehole

Parallel double U-tube
----------------------

::

    from bhr.borehole import Borehole

    double_bhr = Borehole()
    double_bhr.init_double_u_borehole(
        borehole_diameter=0.127,
        pipe_outer_diameter=0.032,
        pipe_dimension_ratio=11,
        length=200,
        shank_space=0.032,
        pipe_conductivity=0.4,
        pipe_inlet_arrangement="DIAGONAL",
        grout_conductivity=1.6,
        soil_conductivity=2.0,
        fluid_type="propylene_glycol",
        fluid_concentration=0.2,
        boundary_condition="UNIFORM_BOREHOLE_WALL_TEMP",
    )

    mass_flow_rate = 0.5  # kg/s total; 0.25 kg/s per U-tube
    temperature = 20.0  # Celsius
    print(f"{double_bhr.calc_bh_resist(mass_flow_rate, temperature):0.5f}")

.. image:: images/double-u_adjacent.webp
   :width: 600
   :alt: Cross-section of a double U-tube with adjacent inlet pipes

.. image:: images/double-u_diagonal.webp
   :width: 600
   :alt: Cross-section of a double U-tube with diagonal inlet pipes

Coaxial
-------

::

    from bhr.borehole import Borehole

    coaxial_bhr = Borehole()
    coaxial_bhr.init_coaxial_borehole(
        borehole_diameter=0.127,
        outer_pipe_outer_diameter=0.114,
        outer_pipe_dimension_ratio=17,
        outer_pipe_conductivity=0.4,
        inner_pipe_outer_diameter=0.06,
        inner_pipe_dimension_ratio=11,
        inner_pipe_conductivity=0.4,
        length=200,
        grout_conductivity=1.6,
        soil_conductivity=2.0,
        fluid_type="propylene_glycol",
        fluid_concentration=0.2,
    )

    mass_flow_rate = 0.5  # kg/s
    temperature = 20.0  # Celsius
    print(f"{coaxial_bhr.calc_bh_resist(mass_flow_rate, temperature):0.5f}")

.. image:: images/coaxial.webp
   :width: 600
   :alt: Cross-section of a concentric coaxial borehole

.. toctree::
   :maxdepth: 2
