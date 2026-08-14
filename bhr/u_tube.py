from scp.base_fluid import BaseFluid

from bhr.pipe import Pipe


class UTube(Pipe):
    def __init__(
        self,
        pipe_outer_diameter: float,
        pipe_dimension_ratio: float,
        length: float,
        shank_space: float,
        pipe_conductivity: float,
        fluid_type: str | None = None,
        fluid_concentration: float = 0,
        *,
        fluid: BaseFluid | None = None,
    ):
        super().__init__(
            pipe_outer_diameter,
            pipe_dimension_ratio,
            length * 2,
            pipe_conductivity,
            fluid_type,
            fluid_concentration,
            fluid=fluid,
        )
        self.length = length
        self.shank_space = shank_space
