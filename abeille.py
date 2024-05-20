from __future__ import annotations
from typing import Optional, List, Tuple, Iterable, Dict, TYPE_CHECKING
from abc import ABC, abstractmethod
import numpy as np
from enum import Enum

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


class Lattice(Universe):
    def __init__(self, name: Optional[str] = None, outer_universe: Optional[Universe] = None, origin: Tuple[float, float, float] = (0., 0., 0.), universes: Iterable[Universe] = []):
        super(Lattice, self).__init__(name)
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
                        out += "{:}]".format(uni_id)

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
        out += "        type: tabulated\n"
        out += "        values: {}\n".format(self.values)
        out += "        pdf: {}".format(self.pdf)
        return out


class Source:
    def __init__(self, spatial: SpatialDistribution, direction: DirectionDistribution, energy: EnergyDistribution, weight: float):
        self.spatial = spatial
        self.direction = direction
        self.energy = energy
        self.weight = weight

    def to_string(self):
        out  = "    - spatial: {}\n".format(self.spatial.to_string())
        out += "      energy: {}\n".format(self.energy.to_string())
        out += "      direction: {}\n".format(self.direction.to_string())
        out += "      weight: {}\n".format(self.weight)
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

        out = '  entropy:\n'
        out += '    low: [{}, {}, {}]\n'.format(self.low.x, self.low.y, self.low.z)
        out += '    hi: [{}, {}, {}]\n'.format(self.hi.x, self.hi.y, self.hi.z)
        out += '    shape: [{}, {}, {}]\n\n'.format(self.shape[0], self.shape[1], self.shape[2])

        return out


class TransportOperator(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def to_string(self) -> str:
        pass


class SurfaceTracker(TransportOperator):
    def __init__(self):
        pass

    def to_string(self) -> str:
        out = '  transport-operator:\n'
        out += '    type: surface-tracking\n'
        return out


class DeltaTracker(TransportOperator):
    def __init__(self):
        pass

    def to_string(self) -> str:
        out = '  transport-operator:\n'
        out += '    type: delta-tracking\n'
        return out


class ImplicitLeakageDeltaTracker(TransportOperator):
    def __init__(self):
        pass

    def to_string(self) -> str:
        out = '  transport-operator:\n'
        out += '    type: implicit-leakage-delta-tracking\n'
        return out


class CarterTracker(TransportOperator):
    def __init__(self, xs : Iterable[float], energy : Optional[Iterable[float]]):
        self.xs = xs
        self.energy = energy

    def to_string(self) -> str:
        out = '  transport-operator:\n'
        out += '    type: carter-tracking\n'
        out += '    xs: {}\n'.format(self.xs)
        if self.energy is not None:
            out += '    energy: {}\n'.format(self.energy)
        return out


class CollisionOperator(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def to_string(self) -> str:
        pass


class BranchingCollisions(CollisionOperator):
    def __init__(self):
        pass

    def to_string(self) -> str:
        out = '  collision-operator:\n'
        out += '    type: branching\n'
        return out


class BranchlessIosotopeCollisions(CollisionOperator):
    def __init__(self, splitting: bool = False):
        self.splitting = splitting

    def to_string(self) -> str:
        out = '  collision-operator:\n'
        out += '    type: branchless-isotope\n'
        if self.splitting:
            out += '    splitting: {}\n'.format(self.splitting)
        return out


class BranchlessMaterialCollisions(CollisionOperator):
    def __init__(self, splitting: bool = False):
        self.splitting = splitting

    def to_string(self) -> str:
        out = '  collision-operator:\n'
        out += '    type: branchless-material\n'
        if self.splitting:
            out += '    splitting: {}\n'.format(self.splitting)
        return out


class Simulation(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def to_string(self) -> str:
        pass


class FixedSource(Simulation):
    def __init__(self, nparticles : int, nbatches : int, sources : List[Source], transport_op : TransportOperator = SurfaceTracker, collision_op : CollisionOperator = SurfaceTracker):
        self.nparticles = nparticles
        self.nbatches = nbatches
        self.transport_operator = transport_op
        self.collision_operator = collision_op
        self.sources = sources

    def to_string(self) -> str:
        out = 'simulation:\n'
        out += '  mode: fixed-source\n'
        out += '  nparticles: {}\n'.format(self.nparticles)
        out += '  nbatches: {}\n\n'.format(self.nbatches)
        out += self.transport_operator.to_string() + '\n'
        out += self.collision_operator.to_string() + '\n'
        out += "  sources:\n"
        for src in self.sources:
            out += src.to_string()
            out += "\n"
        return out


class PowerIterator(Simulation):
    def __init__(self, nparticles : int, ngenerations : int, nignored : int, sources : List[Source], transport_op : TransportOperator = SurfaceTracker(), collision_op : CollisionOperator = BranchingCollisions()):
        self.nparticles = nparticles
        self.ngenerations = ngenerations
        self.nignored = nignored
        self.combing = False
        self.families = False
        self.pair_distance_sqrd = False
        self.empty_entropy_bins = False
        self.transport_operator = transport_op
        self.collision_operator = collision_op
        self.entropy = None
        self.cancelator = None
        self.sources = sources

    def to_string(self) -> str:
        out = 'simulation:\n'
        out += '  mode: k-eigenvalue\n'
        out += '  nparticles: {}\n'.format(self.nparticles)
        out += '  ngenerations: {}\n'.format(self.ngenerations)
        out += '  nignored: {}\n'.format(self.nignored)

        if self.combing:
            out += '  combing: {}\n'.format(self.combing)

        if self.families:
            out += '  families: {}\n'.format(self.families)

        if self.pair_distance_sqrd:
            out += '  pair-distance-sqrd: {}\n'.format(self.pair_distance_sqrd)

        if self.empty_entropy_bins:
            out += '  empty-entropy-bins: {}\n'.format(self.empty_entropy_bins)
        
        out += '\n'
        out += self.transport_operator.to_string() + '\n'
        out += self.collision_operator.to_string() + '\n'
        
        if self.entropy is not None and not isinstance(self.entropy, Entropy):
            raise RuntimeError("Provided entropy entry is an unknown type")
        elif self.entropy is not None:
            out += self.entropy.to_string()

        out += "  sources:\n"
        for src in self.sources:
            out += src.to_string()
            out += "\n"
        return out


class Settings:
    __energy_modes = ['continuous-energy', 'multi-group']
    __temperature_interp = ['exact', 'nearest', 'linear']

    def __init__(self):
        # Energy mode settings
        self.energy_mode = 'continuous-energy'
        self.ngroups = 0
        self.energy_bounds = []

        self.max_run_time = None # Default is inf

        # RNG settings
        self.rng_stride_warnings = None # Bool, default is True
        self.seed = None # Default is 19073486328125
        self.stride = None # Default is 152917
        
        # Roulette settings
        self.wgt_cutoff = None # Default is 0.25
        self.wgt_survival = None # Default is 1.0
        self.wgt_split = None # Default is 2.0

        # General physics settings
        self.dbrc_nuclides = []
        self.use_urr_ptables = True
        self.temperature_interpolation = 'linear' # Default
        self.nuclear_data = None # Default then to environment variable

    def to_string(self) -> str:
        out = 'settings:\n'

        if self.energy_mode not in Settings.__energy_modes:
            raise RuntimeError('Unknown energy mode {}'.format(self.energy_mode))
        if self.energy_mode == 'multi-group':
            out += '  energy-mode: {}\n'.format(self.energy_mode)

        if self.max_run_time is not None:
            out += '  max-run-time: {}\n'.format(self.max_run_time)

        if self.rng_stride_warnings is not None:
            out += '  rng-stride-warnings: {}\n'.format(self.rng_stride_warnings)

        if self.seed is not None:
            out += '  seed: {}\n'.format(self.seed)

        if self.stride is not None:
            out += '  stride: {}\n'.format(self.stride)

        if self.wgt_cutoff is not None:
            out += '  wgt-cutoff: {}\n'.format(self.wgt_cutoff)
        
        if self.wgt_survival is not None:
            out += '  wgt-survival: {}\n'.format(self.wgt_survival)

        if self.wgt_split is not None:
            out += '  wgt-split: {}\n'.format(self.wgt_split)

        # General physics
        if self.energy_mode == 'continuous-energy':
            if len(self.dbrc_nuclides) > 0:
                out += '  use-dbrc: True\n'
                out += '  dbrc-nuclides: {}\n'.format(self.dbrc_nuclides)
            
            if self.use_urr_ptables == False:
                out += '  use-urr-ptables: {}\n'.format(self.use_urr_ptables)
            
            if self.temperature_interpolation not in Settings.__temperature_interp:
                raise RuntimeError('Unkown temperature-interpolation {}'.format(self.temperature_interpolation))
            if self.temperature_interpolation != 'linear':
                out += '  temperature-interpolation: {}\n'.format(self.temperature_interpolation)
            
            if self.nuclear_data is not None:
                out += '  nuclear-data: {}\n'.format(self.nuclear_data)
        
        elif self.energy_mode == 'multi-group':
            if self.ngroups == 0:
                raise RuntimeError('Cannot have 0 energy groups in multi-group mode')
            
            if len(self.energy_bounds) != self.ngroups + 1:
                raise RuntimeError('Number of energy bounds does not agree with the number of energy groups')

            out += '  ngroups: {}\n'.format(self.ngroups)
            out += '  energy-bounds: {}\n'.format(self.energy_bounds)

        if out == 'settings:\n':
            return ''

        out += '\n'
        return out


class Input:
    def __init__(self, root_universe: Universe, settings: Settings, simulation : Simulation):
        self.root_universe = root_universe
        self.settings = settings
        self.simulation = simulation

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

        # Write settings
        input_str += self.settings.to_string()

        # Write simulation
        input_str += self.simulation.to_string()

        fl = open(fname, 'w')
        fl.write(input_str)
        fl.close()
