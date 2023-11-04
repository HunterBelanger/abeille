from __future__ import annotations
from typing import Optional, List, Tuple, Iterable, Dict, TYPE_CHECKING
from abc import ABC, abstractmethod
import numpy as np

class Material:
    _id_counter = 1

    def __init__(self, name: str='', temperature: float=293.6, density_units: str = 'sum', fractions: str = 'atoms', density: float = 0.):
        self.name = name
        self.id = Material._id_counter
        Material._id_counter += 1
        self.temperature = temperature
        self.density = density
        self.density_units = density_units
        self.fractions = fractions
        self.nuclides = []

    def _update_density(self):
        if self.density_units == 'sum':
            self.density = 0.
            for nuc in self.nuclides:
                self.density += nuc[1]
            
    def to_string(self):
        out = "  - id: {:}\n".format(self.id)
        if self.name is not None:
            out += "    name: {:}\n".format(self.name)
        out += "    temperature: {:}\n".format(self.temperature)
        out += "    density-units: {:}\n".format(self.density_units)
        if self.density_units != 'g/cm3' and self.density_units != 'atoms/b-cm' and self.density_units != 'sum':
            raise RuntimeError("Density units can only be 'g/cm3', 'atomcs/b-cm', or 'sum'.")
        if self.density_units != 'sum':
            out += "    density: {:}\n".format(self.density)
        out += "    fractions: {:}\n".format(self.fractions)
        out += "    composition: ["
        for i, nuc in enumerate(self.nuclides):
            out += '{'
            out += "nuclide: {:}, fraction: {:}".format(nuc[0], nuc[1])
            out += '}'
            if i < len(self.nuclides)-1:
                out += ",\n                  "
        out += "]"
        return out


    def add_nuclide(self, id: str, fraction: float):
        if fraction < 0.:
            raise RuntimeError("Nuclide fraction must be >= 0.")

        self.nuclides.append((id, fraction))
        self._update_density()

    def set_density(self, units: str, density: Optional[float] = None):
        if units != 'g/cm3' and units != 'atoms/b-cm' and units != 'sum':
            raise RuntimeError("Density units can only be 'g/cm3', 'atomcs/b-cm', or 'sum'.")

        self.density_units = units

        if self.density_units == 'g/cm3' or self.density_units == 'atoms/b-cm':
            if density is None:
                raise RuntimeError("Must define a density when using units of 'g/cm3' or 'atomcs/b-cm'.")

            if density < 0.:
                raise RuntimeError("Density must be >= 0.")

            self.density = density

        self._update_density()

    def set_temperature(self, temperature: float):
        if temperature <= 0.:
            raise RuntimeError("Material temperature must be > 0.")
        self.temperature = temperature


class Region(ABC):
    def __init__(self):
        pass

    def __and__(self, other) -> Intersection:
        return Intersection(self, other)
    
    def __or__(self, other) -> Union:
        return Union(self, other)

    def __invert__(self) -> Complement:
        return Complement(self)

    def to_string(self) -> str:
        return ""

    def add_surfaces(self, surfs: Dict[int, Surface]) -> None:
        pass


class Halfspace(Region):
    def __init__(self, side: bool, surface: Surface):
        self.side = side
        self.surface = surface

    def to_string(self) -> str:
        if self.side:
            s = '+'
        else:
            s = '-'
        return "{:}{:}".format(s, self.surface.id)

    def add_surfaces(self, surfs: Dict[int, Surface]) -> None:
        if not self.surface.id in surfs:
            surfs[self.surface.id] = self.surface


class Intersection(Region):
    def __init__(self, left: Region, right: Region):
        self.left = left
        self.right = right

    def to_string(self) -> str:
        return "{:} & {:}".format(self.left.to_string(), self.right.to_string())

    def add_surfaces(self, surfs: Dict[int, Surface]) -> None:
        self.left.add_surfaces(surfs)
        self.right.add_surfaces(surfs)


class Union(Region):
    def __init__(self, left: Region, right: Region):
        self.left = left
        self.right = right

    def to_string(self) -> str:
        return "{:} U {:}".format(self.left.to_string(), self.right.to_string())

    def add_surfaces(self, surfs: Dict[int, Surface]) -> None:
        self.left.add_surfaces(surfs)
        self.right.add_surfaces(surfs)


class Complement(Region):
    def __init__(self, region: Region):
        self.region = region

    def to_string(self) -> str:
        return "~({:})".format(self.region.to_string())

    def add_surfaces(self, surfs: Dict[int, Surface]) -> None:
        self.region.add_surfaces(surfs)


class Surface(ABC):
    _id_counter = 1

    def __init__(self, name: Optional[str] = None, boundary_type: str = 'normal'):
        self.id = Surface._id_counter
        Surface._id_counter += 1
        self.name = name

        if boundary_type not in ['normal', 'vacuum', 'reflective']:
            raise RuntimeError("Boundary type must be 'normal', 'vacuum', or 'reflective'.")
        self.boundary_type = boundary_type

    def __pos__(self) -> Halfspace:
        return Halfspace(True, self)

    def __neg__(self) -> Halfspace:
        return Halfspace(False, self)

    @abstractmethod
    def to_string(self) -> str:
        pass


class XPlane(Surface):
    def __init__(self, x0: float, name: Optional[str] = None, boundary_type: str = 'normal'):
        self.x0 = x0
        super(XPlane, self).__init__(name, boundary_type)

    def to_string(self) -> str:
        out = "id: {:}, type: xplane, x0: {:}".format(self.id, self.x0)
        if self.boundary_type != 'normal':
            out += ", boundary: {:}".format(self.boundary_type)
        if self.name is not None:
            out += ", name: {:}".format(self.name)
        return "  - {"+out+"}"


class YPlane(Surface):
    def __init__(self, y0: float, name: Optional[str] = None, boundary_type: str = 'normal'):
        self.y0 = y0
        super(YPlane, self).__init__(name, boundary_type)

    def to_string(self) -> str:
        out = "id: {:}, type: yplane, y0: {:}".format(self.id, self.y0)
        if self.boundary_type != 'normal':
            out += ", boundary: {:}".format(self.boundary_type)
        if self.name is not None:
            out += ", name: {:}".format(self.name)
        return "  - {"+out+"}"


class ZPlane(Surface):
    def __init__(self, z0: float, name: Optional[str] = None, boundary_type: str = 'normal'):
        self.z0 = z0
        super(ZPlane, self).__init__(name, boundary_type)

    def to_string(self) -> str:
        out = "id: {:}, type: zplane, z0: {:}".format(self.id, self.z0)
        if self.boundary_type != 'normal':
            out += ", boundary: {:}".format(self.boundary_type)
        if self.name is not None:
            out += ", name: {:}".format(self.name)
        return "  - {"+out+"}"


class Plane(Surface):
    def __init__(self, A: float, B: float, C: float, D: float, name: Optional[str] = None, boundary_type: str = 'normal'):
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        super(Plane, self).__init__(name, boundary_type)

    def to_string(self) -> str:
        out = "id: {:}, type: plane, A: {:}, B: {:}, C: {:}, D: {:}".format(self.id, self.A, self.B, self.C, self.D)
        if self.boundary_type != 'normal':
            out += ", boundary: {:}".format(self.boundary_type)
        if self.name is not None:
            out += ", name: {:}".format(self.name)
        return "  - {"+out+"}"


class Sphere(Surface):
    def __init__(self, r: float, x0: float = 0., y0: float = 0., z0: float = 0., name: Optional[str] = None, boundary_type: str = 'normal'):
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.r = r
        super(Sphere, self).__init__(name, boundary_type)

    def to_string(self) -> str:
        out = "id: {:}, type: sphere, x0: {:}, y0: {:}, z0: {:}, r: {:}".format(self.id, self.x0, self.y0, self.z0, self.r)
        if self.boundary_type != 'normal':
            out += ", boundary: {:}".format(self.boundary_type)
        if self.name is not None:
            out += ", name: {:}".format(self.name)
        return "  - {"+out+"}"


class XCylinder(Surface):
    def __init__(self, r: float, y0: float = 0., z0: float = 0., name: Optional[str] = None, boundary_type: str = 'normal'):
        self.y0 = y0
        self.z0 = z0
        self.r = r
        super(XCylinder, self).__init__(name, boundary_type)

    def to_string(self) -> str:
        out = "id: {:}, type: xcylinder, y0: {:}, z0: {:}, r: {:}".format(self.id, self.y0, self.z0, self.r)
        if self.boundary_type != 'normal':
            out += ", boundary: {:}".format(self.boundary_type)
        if self.name is not None:
            out += ", name: {:}".format(self.name)
        return "  - {"+out+"}"


class YCylinder(Surface):
    def __init__(self, r: float, x0: float = 0., z0: float = 0., name: Optional[str] = None, boundary_type: str = 'normal'):
        self.x0 = x0
        self.z0 = z0
        self.r = r
        super(YCylinder, self).__init__(name, boundary_type)

    def to_string(self) -> str:
        out = "id: {:}, type: ycylinder, x0: {:}, z0: {:}, r: {:}".format(self.id, self.x0, self.z0, self.r)
        if self.boundary_type != 'normal':
            out += ", boundary: {:}".format(self.boundary_type)
        if self.name is not None:
            out += ", name: {:}".format(self.name)
        return "  - {"+out+"}"


class ZCylinder(Surface):
    def __init__(self, r: float, x0: float = 0., y0: float = 0., name: Optional[str] = None, boundary_type: str = 'normal'):
        self.x0 = x0
        self.y0 = y0
        self.r = r
        super(ZCylinder, self).__init__(name, boundary_type)

    def to_string(self) -> str:
        out = "id: {:}, type: zcylinder, x0: {:}, y0: {:}, r: {:}".format(self.id, self.x0, self.y0, self.r)
        if self.boundary_type != 'normal':
            out += ", boundary: {:}".format(self.boundary_type)
        if self.name is not None:
            out += ", name: {:}".format(self.name)
        return "  - {"+out+"}"


class Cell:
    _id_counter = 1
    def __init__(self, region: Region = Region(), fill = None, name: Optional[str] = None):
        self.id = Cell._id_counter
        Cell._id_counter += 1
        self.region = region
        self.name = name
        self.fill = fill

    def to_string(self) -> str:
        out = "id: {:}, region: \"{:}\"".format(self.id, self.region.to_string())
        if isinstance(self.fill, Material):
            out += ", material: {:}".format(self.fill.id)
        elif isinstance(self.fill, Universe):
            out += ", universe: {:}".format(self.fill.id)
        else:
            raise RuntimeError("Invalid cell fill type " + str(type(self.fill)) + ".")

        if self.name is not None:
            out += ", name: {:}".format(self.name)
        return "  - {"+out+"}"

    def add_surfaces(self, surfs: Dict[int, Surface]) -> None:
        self.region.add_surfaces(surfs)

    def add_cells(self, cells: Dict[int, Cell]) -> None:
        if isinstance(self.fill, Universe):
            self.fill.add_cells(cells)

    def add_universes(self, unis: Dict[int, Universe]) -> None:
        if isinstance(self.fill, Universe):
            if self.fill.id not in unis:
                unis[self.fill.id] = self.fill
                self.fill.add_universes(unis)

    def add_lattices(self, lats: Dict[int, Lattice]) -> None:
        if isinstance(self.fill, Universe):
            self.fill.add_lattices(lats)


class Universe(ABC):
    _id_counter = 1

    def __init__(self, name: Optional[str] = None):
        self.id = Universe._id_counter
        Universe._id_counter += 1
        self.name = name

    @abstractmethod
    def to_string(self) -> str:
        pass
    
    @abstractmethod
    def add_surfaces(self, surfs: Dict[int, Surface]) -> None:
        pass

    @abstractmethod
    def add_cells(self, cells: Dict[int, Cell]) -> None:
        pass

    @abstractmethod 
    def add_lattices(self, lats: Dict[int, Lattice]) -> None:
        pass
    
    @abstractmethod
    def add_universes(self, unis: Dict[int, Universe]) -> None:
        pass


class Lattice(ABC):
    _id_counter = 1

    def __init__(self, name: Optional[str] = None, outer_universe: Optional[Universe] = None, origin: Tuple[float, float, float] = (0., 0., 0.), universes: Iterable[Universe] = []):
        self.id = Lattice._id_counter
        Lattice._id_counter += 1
        self.name = name
        self.outer_universe = outer_universe
        self.origin = origin
        self.universes = universes

    @abstractmethod
    def to_string(self) -> str:
        pass

    @abstractmethod
    def add_surfaces(self, surfs: Dict[int, Surface]) -> None:
        pass

    @abstractmethod
    def add_cells(self, cells: Dict[int, Cell]) -> None:
        pass

    @abstractmethod 
    def add_lattices(self, lats: Dict[int, Lattice]) -> None:
        pass
    
    @abstractmethod
    def add_universes(self, unis: Dict[int, Universe]) -> None:
        pass


class CellUniverse(Universe):
    def __init__(self, cells: List[Cell], name: Optional[str] = None):
        self.cells = cells
        super(CellUniverse, self).__init__(name)

    def to_string(self):
        out = "id: {:}, cells: [".format(self.id)
        for i in range(len(self.cells)):
            out += str(self.cells[i].id)
            if i != len(self.cells)-1:
                out += ", "
        out += "]"
        if self.name is not None:
            out += ", name: {:}".format(self.name)
        return "  - {"+out+"}"
    
    def add_surfaces(self, surfs: Dict[int, Surface]) -> None:
        for cell in self.cells:
            cell.add_surfaces(surfs)

    def add_cells(self, cells: Dict[int, Cell]) -> None:
        for cell in self.cells:
            if cell.id not in cells:
                cells[cell.id] = cell
                cell.add_cells(cells)

    def add_lattices(self, lats: Dict[int, Lattice]) -> None:
        for cell in self.cells:
            cell.add_lattices(lats)
    
    def add_universes(self, unis: Dict[int, Universe]) -> None:
        for cell in self.cells:
            cell.add_universes(unis)


class LatticeUniverse(Universe):
    def __init__(self, lattice: Lattice, name: Optional[str] = None):
        self.lattice = lattice
        super(LatticeUniverse, self).__init__(name)

    def to_string(self):
        out = "id: {:}, lattice: {:}".format(self.id, self.lattice.id)
        if self.name is not None:
            out += ", name: {:}".format(self.name)
        return "  - {"+out+"}"

    def add_surfaces(self, surfs: Dict[int, Surface]) -> None:
        self.lattice.add_surfaces(surfs)

    def add_cells(self, cells: Dict[int, Cell]) -> None:
        self.lattice.add_cells(cells)

    def add_lattices(self, lats: Dict[int, Lattice]) -> None:
        if self.lattice.id not in lats:
            lats[self.lattice.id] = self.lattice
            self.lattice.add_lattices(lats)
    
    def add_universes(self, unis: Dict[int, Universe]) -> None:
        self.lattice.add_universes(unis)


class RectLattice(Lattice):
    def __init__(self, shape: Tuple[int, int, int], pitch: Tuple[float, float, float], origin: Tuple[float, float, float] = (0., 0., 0.), name: Optional[str] = None, outer_universe: Optional[Universe] = None, universes: Iterable[Universe] = []):
        if len(shape) != 3:
            raise RuntimeError("HexLattice shape must have 2 dimensions.")

        if shape[0] < 1 or shape[1] < 1 or shape[2] < 1:
            raise RuntimeError("All elements of HexLattice shape must be >= 1.")
        self.shape = shape

        if len(pitch) != 3:
            raise RuntimeError("RectLattice pitch must have 3 dimensions.")

        if pitch[0] <= 0. or pitch[1] <= 0. or pitch[2] <= 0.:
            raise RuntimeError("All elements of HexLattice pitch must > 0.")

        self.pitch = pitch
        super(RectLattice, self).__init__(name, outer_universe, origin, universes)

    def to_string(self) -> str:
        out = "  - id: {:}\n".format(self.id)
        if self.name is not None:
            out += "    name: {:}\n".format(self.name)
        out += "    type: rectlinear\n"
        out += "    shape: [{:}, {:}, {:}]\n".format(self.shape[0], self.shape[1], self.shape[2])
        out += "    pitch: [{:}, {:}, {:}]\n".format(self.pitch[0], self.pitch[1], self.pitch[2])
        out += "    origin: [{:}, {:}, {:}]\n".format(self.origin[0], self.origin[1], self.origin[2])
        if self.outer_universe is not None:
            out += "    outer: {:}\n".format(self.outer_universe.id)
        out += "    universes: ["
        indx = 0
        for z in range(self.shape[2]):
            for y in range(self.shape[1]):
                for x in range(self.shape[0]):
                    uni_id = -1
                    if isinstance(self.universes[indx], Universe):
                        uni_id = self.universes[indx].id

                    if indx < len(self.universes)-1:
                        out += "{:},".format(uni_id)
                    else:
                        out += "{:}]\n".format(uni_id)

                    indx += 1

                if indx < len(self.universes)-1:
                    out += "\n                "
                #else:
                #    out += "\n"
            if z < self.shape[1]-1:
                out += "\n                "

        return out

    def add_surfaces(self, surfs: Dict[int, Surface]) -> None:
        if self.outer_universe is not None:
            self.outer_universe.add_surfaces(surfs)
        
        for uni in self.universes:
            if uni is not None:
                uni.add_surfaces(surfs)

    def add_cells(self, cells: Dict[int, Cell]) -> None:
        if self.outer_universe is not None:
            self.outer_universe.add_cells(cells)
        
        for uni in self.universes:
            if uni is not None:
                uni.add_cells(cells)

    def add_lattices(self, lats: Dict[int, Lattice]) -> None:
        if self.outer_universe is not None:
            self.outer_universe.add_lattices(lats)
        
        for uni in self.universes:
            if uni is not None:
                uni.add_lattices(lats)
    
    def add_universes(self, unis: Dict[int, Universe]) -> None:
        if self.outer_universe is not None:
            if self.outer_universe.id not in unis:
                unis[self.outer_universe.id] = self.outer_universe
                self.outer_universe.add_universes(unis)
        
        for uni in self.universes:
            if uni is not None:
                if uni.id not in unis:
                    unis[uni.id] = uni
                    uni.add_universes(unis)


"""
class HexLattice(Lattice):
    def __init__(self, shape: Tuple[int, int], pitch: Tuple[float, float], origin: Tuple[float, float, float] = (0., 0., 0.), top: str = 'pointy', name: Optional[str] = None, outer_universe: Optional[Universe] = None, universes: Iterable[Universe] = []):
        if len(shape) != 2:
            raise RuntimeError("HexLattice shape must have 2 dimensions.")

        if shape[0] < 1 or shape[1] < 1:
            raise RuntimeError("All elements of HexLattice shape must be >= 1.")
        self.shape = shape

        if len(pitch) != 2:
            raise RuntimeError("HexLattice pitch must have 2 dimensions.")

        if pitch[0] <= 0. or pitch[1] <= 0.:
            raise RuntimeError("All elements of HexLattice pitch must > 0.")
        self.pitch = pitch

        if top not in ['pointy', 'flat']:
            raise RuntimeError("HexLattice top can either be 'point' or 'flat'.")
        self.top = top

        super(HexLattice, self).__init__(name, outer_universe, origin, universes)
"""

class SpatialDistribution(ABC):
    def __init__(self, fissile_only: bool):
        self.fissile_only = fissile_only
    
    @abstractmethod
    def to_string(self) -> str:
        pass


class Point(SpatialDistribution):
    def __init__(self, x: float, y: float, z: float, fissile_only: bool = False):
        self.x = x
        self.y = y
        self.z = z

        super(Point, self).__init__(fissile_only)

    def to_string(self) -> str:
        out = "type: point, position: [{}, {}, {}], fissile-only: {}".format(self.x, self.y, self.z, self.fissile_only)
        return '{'+out+'}'


class Box(SpatialDistribution):
    def __init__(self, low: Point, hi: Point, fissile_only: bool = False):
        self.low = low
        self.hi = hi

        super(Box, self).__init__(fissile_only)

    def to_string(self) -> str:
        out = "type: box, low: [{}, {}, {}], hi: [{}, {}, {}], fissile-only: {}".format(self.low.x, self.low.y, self.low.z, self.hi.x, self.hi.y, self.hi.z, self.fissile_only)
        return '{'+out+'}'


class DirectionDistribution(ABC):
    def __init__(self):
        pass
    
    @abstractmethod
    def to_string(self) -> str:
        pass


class Isotropic(DirectionDistribution):
    def __init__(self):
        pass

    def to_string(self):
        return "{type: isotropic}"


class MonoDirectional(DirectionDistribution):
    def __init__(self, x: float, y: float, z: float):
        self.x = x
        self.y = y
        self.z = z
    
    def to_string(self) -> str:
        out = "type: mono-directional, direction: [{}, {}, {}]".format(self.x, self.y, self.z)
        return '{'+out+'}'


# aperture is in radians, and is HALF of the beam
class Cone(DirectionDistribution):
    def __init__(self, x: float, y: float, z: float, aperture: float):
        self.x = x
        self.y = y
        self.z = z
        self.aperture = aperture
    
    def to_string(self) -> str:
        out = "type: cone, direction: [{}, {}, {}], aperture: {}".format(self.x, self.y, self.z, self.aperture)
        return '{'+out+'}'


class EnergyDistribution(ABC):
    def __init__(self):
        pass
    
    @abstractmethod
    def to_string(self) -> str:
        pass


class MonoEnergetic(EnergyDistribution):
    def __init__(self, energy: float):
        self.energy = energy

    def to_string(self) -> str:
        out = "type: mono-energetic, energy: {}".format(self.energy)
        return '{'+out+'}'


class Maxwellian(EnergyDistribution):
    def __init__(self, a: float):
        self.a = a

    def to_string(self) -> str:
        out = "type: maxwellian, a: {}".format(self.a)
        return '{'+out+'}'


class Watt(EnergyDistribution):
    def __init__(self, a: float, b: float):
        self.a = a
        self.b = b

    def to_string(self) -> str:
        out = "type: watt, a: {}, b: {}".format(self.a, self.b)
        return '{'+out+'}'


class Tabulated(EnergyDistribution):
    def __init__(self, values: List[float], pdf: List[float]):
        self.values = values
        self.pdf = pdf

    def to_string(self) -> str:
        out = "\n"
        out += "      type: tabulated\n"
        out += "      values: {}\n".format(self.values)
        out += "      pdf: {}".format(self.pdf)
        return out


class Source:
    def __init__(self, spatial: SpatialDistribution, direction: DirectionDistribution, energy: EnergyDistribution, weight: float):
        self.spatial = spatial
        self.direction = direction
        self.energy = energy
        self.weight = weight

    def to_string(self):
        out  = "  - spatial: {}\n".format(self.spatial.to_string())
        out += "    energy: {}\n".format(self.energy.to_string())
        out += "    direction: {}\n".format(self.direction.to_string())
        out += "    weight: {}\n".format(self.weight)
        return out

class Entropy:
    def __init__(self, low: Point, hi: Point, shape: Tuple[int, int, int]):
        self.low = low
        self.hi = hi
        self.shape = shape

    def to_string(self) -> str:
        if self.low.x >= self.hi.x or self.low.y >= self.hi.y or self.low.z >= self.hi.z:
            raise RuntimeError('Entropy low point must be <= hi point.')

        if self.shape[0] <= 0 or self.shape[1] <= 0 or self.shape[2] <= 0:
            raise RuntimeError('Entropy shape values must all be >= 1.')

        out = 'entropy:\n'
        out += '  low: [{}, {}, {}]\n'.format(self.low.x, self.low.y, self.low.z)
        out += '  hi: [{}, {}, {}]\n'.format(self.hi.x, self.hi.y, self.hi.z)
        out += '  shape: [{}, {}, {}]\n'.format(self.shape[0], self.shape[1], self.shape[2])

        return out


class Tally:
    def __init__(self, low: Point, hi: Point, shape: Tuple[int, int, int], energy_bounds: List[float], quantity: str, name: str, estimator: str='collision'):
        self.low = low
        self.hi = hi
        self.shape = shape
        self.energy_bounds = energy_bounds
        self.name = name
        self.quantity = quantity
        self.estimator = estimator

    def to_string(self):
        if self.shape[0] <= 0 or self.shape[1] <= 0 or self.shape[2] <= 0:
            raise RuntimeError('Entropy shape values must all be >= 1.')

        out = '  - name: {}\n'.format(self.name)
        out += '    low: [{}, {}, {}]\n'.format(self.low.x, self.low.y, self.low.z)
        out += '    hi: [{}, {}, {}]\n'.format(self.hi.x, self.hi.y, self.hi.z)
        out += '    shape: [{}, {}, {}]\n'.format(self.shape[0], self.shape[1], self.shape[2])
        out += '    energy-bounds: {}\n'.format(self.energy_bounds)
        out += '    quantity: {}\n'.format(self.quantity)
        out += '    estimator: {}\n'.format(self.estimator)

        return out


class Settings:
    __simulation_modes = ['k-eigenvalue', 'branchless-k-eigenvalue', 'fixed-source', 'noise']
    __transport_modes = ['surface-tracking', 'delta-tracking', 'carter-tracking', 'implicit-leakage-delta-tracking']
    __energy_modes = ['continuous-energy', 'multi-group']
    __temperature_interp = ['exact', 'nearest', 'linear']

    def __init__(self):
        # General run settings
        self.simulation = 'k-eigenvalue'
        self.transport = 'surface-tracking'
        self.energy_mode = 'continuous-energy' # Default
        self.nparticles = 10000
        self.ngenerations = 2100
        self.nignored = 100
        self.max_run_time = None # Default is inf
        self.rng_stride_warnings = None # Bool, default is True
        self.seed = None # Default is 19073486328125
        self.stride = None # Default is 152917
        self.cancellation = None # Bool

        self.wgt_cutoff = None # Default is 0.25
        self.wgt_survival = None # Default is 1.0
        self.wgt_split = None # Default is 2.0

        # Noise settings
        self.noise_angular_frequency = 2.*np.pi
        self.nskip = 10 # Number of generations skipped between noise sampling
        self.keff = None
        self.inner_generations = True
        self.normalize_noise_source = True
        self.noise_cancellation = False
        self.cancel_noise_gens = None # Int

        # Power Iteration settings
        self.pair_distance_sqrd = False
        self.nfamilies = False
        self.empty_entropy_bins = False

        # MG Settings
        self.ngroups = 0
        self.energy_bounds = []

        # Carter tracking settings
        self.sampling_xs_ratio = []
        
        # General physics settings
        self.use_dbrc = True
        self.dbrc_nuclides = []
        self.use_urr_ptables = True
        self.temperature_interpolation = 'linear' # Default
        self.nuclear_data = None # Default then to environment variable

        # Branchless settings
        self.branchless_combing = True
        self.branchless_material = True
        self.branchless_splitting = False

    def to_string(self) -> str:
        out = 'settings:\n'

        if self.simulation not in Settings.__simulation_modes:
            raise RuntimeError('Unkown simualtion mode {}'.format(self.simulation))
        out += '  simulation: {}\n'.format(self.simulation)

        if self.transport not in Settings.__transport_modes:
            raise RuntimeError('Unkown transport mode {}'.format(self.transport))
        out += '  transport: {}\n'.format(self.transport)

        if self.energy_mode not in Settings.__energy_modes:
            raise RuntimeError('Unknown energy mode {}'.format(self.energy_mode))
        if self.energy_mode == 'multi-group':
            out += '  energy-mode: {}\n'.format(self.energy_mode)

        out += '  nparticles: {}\n'.format(self.nparticles)
        out += '  ngenerations: {}\n'.format(self.ngenerations)
        
        if self.simulation in ['k-eigenvalue', 'branchless-k-eigenvalue', 'noise']:
            out += '  nignored: {}\n'.format(self.nignored)

        if self.max_run_time is not None:
            out += '  max-run-time: {}\n'.format(self.max_run_time)

        if self.rng_stride_warnings is not None:
            out += '  rng-stride-warnings: {}\n'.format(self.rng_stride_warnings)

        if self.seed is not None:
            out += '  seed: {}\n'.format(self.seed)

        if self.stride is not None:
            out += '  stride: {}\n'.format(self.stride)

        if self.cancellation is not None:
            out += '  cancellation: {}\n'.format(self.cancellation)

        if self.wgt_cutoff is not None:
            out += '  wgt-cutoff: {}\n'.format(self.wgt_cutoff)
        
        if self.wgt_survival is not None:
            out += '  wgt-survival: {}\n'.format(self.wgt_survival)

        if self.wgt_split is not None:
            out += '  wgt-split: {}\n'.format(self.wgt_split)

        # General physics
        if self.energy_mode == 'continuous-energy':
            if self.use_dbrc == False:
                out += '  use-dbrc: {}\n'.format(self.use_dbrc)
            elif len(self.dbrc_nuclides) > 0:
                out += '  dbrc-nuclides: {}\n'.format(self.dbrc_nuclides)
            
            if self.use_urr_ptables == False:
                out += '  use-urr-ptables: {}\n'.format(self.use_urr_ptables)
            
            if self.temperature_interpolation not in Settings.__temperature_interp:
                raise RuntimeError('Unkown temperature-interpolation {}'.format(self.temperature_interpolation))
            if self.temperature_interpolation != 'linear':
                out += '  temperature-interpolation: {}\n'.format(self.temperature_interpolation)
            
            if self.nuclear_data is not None:
                out += '  nuclear-data: {}\n'.format(self.nuclear_data)

        if self.simulation == 'noise':
            out += '  noise-angular-frequency: {}\n'.format(self.noise_angular_frequency)
            out += '  nskip: {}\n'.format(self.nskip)

            if self.keff is None:
                raise RuntimeError('Must assign keff in Settings for noise simulations.')
            out += '  keff: {}\n'.format(self.keff)

            out += '  inner-generations: {}\n'.format(self.inner_generations)
            out += '  normalize-noise-source: {}\n'.format(self.normalize_noise_source)
            out += '  noise-cancellation: {}\n'.format(self.noise_cancellation)

            if self.cancel_noise_gens is not None:
                out += '  cancel-noise-gens: {}\n'.format(self.cancel_noise_gens)
        
        elif self.simulation in ['k-eigenvalue', 'branchless-k-eigenvalue']:
            out += '  pair-distance-sqrd: {}\n'.format(self.pair_distance_sqrd)
            out += '  nfamilies: {}\n'.format(self.nfamilies)
            out += '  empty-entropy-bins: {}\n'.format(self.empty_entropy_bins)

            if self.simulation == 'branchless-k-eigenvalue':
                out += '  branchless-combing: {}\n'.format(self.branchless_combing)
                out += '  branchless-splitting: {}\n'.format(self.branchless_splitting)
                out += '  branchless-material: {}\n'.format(self.branchless_material)

        if self.transport == 'carter-tracking':
            out += '  sampling-xs-ratio: {}\n'.format(self.sampling_xs_ratio)
        
        if self.energy_mode == 'multi-group':
            out += '  ngroups: {}\n'.format(self.ngroups)
            out += '  energy-bounds: {}\n'.format(self.energy_bounds)

        return out


class Input:
    def __init__(self, root_universe: Universe, sources: List[Source], settings: Settings, tallies: List[Tally]=[], entropy: Optional[Entropy] = None):
        self.root_universe = root_universe
        self.sources = sources
        self.settings = settings
        self.tallies = tallies
        self.entropy = entropy

    def to_file(self, fname: str):
        # Get everything we need for the input
        surfaces = {}
        self.root_universe.add_surfaces(surfaces)

        cells = {}
        self.root_universe.add_cells(cells)

        lattices = {}
        self.root_universe.add_lattices(lattices)

        universes = {}
        self.root_universe.add_universes(universes)
        universes[self.root_universe.id] = self.root_universe

        materials = {}
        for cell_id in cells:
            if isinstance(cells[cell_id].fill, Material):
                if cells[cell_id].fill.id not in materials:
                    materials[cells[cell_id].fill.id] = cells[cell_id].fill
        
        input_str = ""
        # We first write materials
        if len(materials) == 0:
            raise RuntimeError("No materials present in root universe.")
        input_str += "materials:\n"
        for mat_id in sorted(materials.keys()):
            input_str += "{:}\n\n".format(materials[mat_id].to_string())

        # Now we do surfaces
        if len(surfaces) == 0:
            raise RuntimeError("No surfaces present in root universe.")
        input_str += "surfaces:\n"
        for surf_id in sorted(surfaces.keys()):
            input_str += "{:}\n".format(surfaces[surf_id].to_string())
        input_str += "\n"

        # Now cells
        if len(cells) == 0:
            raise RuntimeError("No cells present in root universe.")
        input_str += "cells:\n"
        for cell_id in sorted(cells.keys()):
            input_str += "{:}\n".format(cells[cell_id].to_string())
        input_str += "\n"

        # Now universes
        if len(universes) == 0:
            raise RuntimeError("No universes present in root universe.")
        input_str += "universes:\n"
        for uni_id in sorted(universes.keys()):
            input_str += "{:}\n".format(universes[uni_id].to_string())
        input_str += "\n"

        # Finally lattices
        if len(lattices) > 0:
            input_str += "lattices:\n"
            for lat_id in sorted(lattices.keys()):
                input_str += "{:}\n".format(lattices[lat_id].to_string())

        # Write root universe id 
        input_str += "root-universe: {:}\n".format(self.root_universe.id)
        input_str += "\n"

        # Write tallies
        if len(self.tallies) > 0:
            input_str += 'tallies:\n'
            for taly in self.tallies:
                input_str += "{:}\n".format(taly.to_string())

        if self.entropy is not None:
            input_str += self.entropy.to_string()

        # Write sources
        input_str += "sources:\n"
        for src in self.sources:
            input_str += src.to_string()
            input_str += "\n"

        # Write settings
        input_str += self.settings.to_string()

        fl = open(fname, 'w')
        fl.write(input_str)
        fl.close()

if __name__ == "__main__":
    mod = Material("moderator", density=1., density_units='g/cm3')
    mod.add_nuclide("H1", 2.)
    mod.add_nuclide("O16", 1.)

    fuel = Material(name="fuel", density=92., density_units='g/cm3')
    fuel.add_nuclide("U235", 1.)
    fuel.add_nuclide("O16", 2.)

    fuel_rad = ZCylinder(0.4058)
    fuel_cell = Cell(region=-fuel_rad, fill=fuel)
    mod_cell = Cell(region=+fuel_rad, fill=mod)

    p = CellUniverse(cells=[fuel_cell, mod_cell])

    mod_uni_cell = Cell(fill=mod)
    m = CellUniverse(cells=[mod_uni_cell])


    assembly = RectLattice(shape=[17, 17, 2], pitch=[1.26, 1.26, 192.], name="assembly")
    assembly.universes = [p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,p,p,p,m,p,p,m,p,p,m,p,p,p,p,p,
                          p,p,p,m,p,p,p,p,p,p,p,p,p,m,p,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,m,p,p,m,p,p,m,p,p,m,p,p,m,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,m,p,p,m,p,p,m,p,p,m,p,p,m,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,m,p,p,m,p,p,m,p,p,m,p,p,m,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,p,m,p,p,p,p,p,p,p,p,p,m,p,p,p,
                          p,p,p,p,p,m,p,p,m,p,p,m,p,p,p,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,p,p,p,m,p,p,m,p,p,m,p,p,p,p,p,
                          p,p,p,m,p,p,p,p,p,p,p,p,p,m,p,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,m,p,p,m,p,p,m,p,p,m,p,p,m,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,m,p,p,m,p,p,m,p,p,m,p,p,m,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,m,p,p,m,p,p,m,p,p,m,p,p,m,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,p,m,p,p,p,p,p,p,p,p,p,m,p,p,p,
                          p,p,p,p,p,m,p,p,m,p,p,m,p,p,p,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,
                          p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p,p]

    assembly_uni = LatticeUniverse(assembly)


    spatial = Box(Point(0., 0., 0.), Point(1., 2., 3.))
    energy = Watt(2., 3.)
    direc = MonoDirectional(0., 1., 0.)

    sources = [Source(spatial, direc, energy, 1.)]
    
    spatial = Point(0., 0., 0.)
    energy = Tabulated([1., 2., 3.], [0.3, 0.6, 0.1])
    direc = Isotropic()
    sources.append(Source(spatial, direc, energy, 2.))

    settings = Settings()

    inpt = Input(assembly_uni, sources, settings)
    inpt.to_file("fl.yaml")
    