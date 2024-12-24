# Immersed Boundary Method Solver for OpenFOAM

## Overview

This project implements Peskin's continuous forcing immersed boundary method in OpenFOAM. The solver facilitates fluid-structure interaction by managing the coupling between fluid flow and solid deformation through interpolation using discrete Dirac delta functions.

## Method Description

The implementation follows Peskin's immersed boundary method, which uses a Eulerian description for the fluid and a Lagrangian description for the structure. The key features include:

1. **Velocity Interpolation**: Fluid velocities are interpolated from the Eulerian grid to the Lagrangian points representing the solid structure.
2. **Force Spreading**: Forces computed at the Lagrangian points are spread back to the Eulerian grid using the same interpolation kernel.
3. **Structure Movement**: The Lagrangian points are moved according to the interpolated fluid velocities.

The method uses a regularized discrete Dirac delta function for both interpolation operations, ensuring smooth coupling between the fluid and solid domains.

## Code Structure

### Main Components

- `pointSolver.C`: Main solver implementing the IB method
- `createFields.H`: Initializes the velocity (U) and force (F) fields
- `createIbPoints.H`: Sets up the Lagrangian points cloud
- `interpolateForces.H`: Implements force spreading from points to mesh
- `interpolateVelocity.H`: Implements velocity interpolation from mesh to points
- `diracdelta.H`: Implements the discrete Dirac delta function
- `UEqn.H`: Defines the prescribed velocity field (currently a rotation)

### Key Features

#### Particle Cloud Management
```cpp
Cloud<passiveParticle> ibpoints(mesh, "pointsCloud", false);
```
The solid structure is represented by Lagrangian points stored in an OpenFOAM particle cloud.

#### Spatial Search
```cpp
const indexedOctree<treeDataCell>& cellTree = mesh.cellTree();
```
An octree data structure is used for efficient spatial searching when identifying cells near Lagrangian points.

#### Discrete Dirac Delta Function
The implementation uses a three-dimensional regularized delta function with compact support:
```cpp
scalar diracdelta(vector x, scalar h)
```
The function implements a smoothed approximation with support over 2 mesh cells in each direction.

## Current Implementation Status

The current version includes:
- Single point demonstration
- Prescribed rotational velocity field
- Force spreading and velocity interpolation using regularized delta functions
- Point tracking and movement

## Planned Extensions

1. **Fortran FEM Integration**: Integration with an external Fortran library for computing deformation forces of a solid comprising of multiple points using the finite element method.
2. **Fluid Solver Coupling**: Full coupling with an OpenFOAM fluid solver for two-way interaction.

## Usage

### Prerequisites
- OpenFOAM installation
- Compatible C++ compiler

### Running the Solver
1. Compile the solver using wmake
2. Set up the case directory with appropriate boundary and initial conditions
3. Execute the solver using `pointSolver`

## Implementation Details

### Velocity Interpolation
The velocity at each Lagrangian point is computed as:
```cpp
pu[i] = pu[i] + U[icell][i]*diracdelta(ddpoint,hh)*(hh*hh*hh);
```
where `hh` is the mesh spacing and `ddpoint` is the distance vector between the cell center and the Lagrangian point.

### Force Spreading
Forces are spread to the Eulerian grid using:
```cpp
F[icell][i] = F[icell][i] + pf[i]*diracdelta(ddpoint,hh);
```

### Point Movement
Points are moved using the interpolated velocity:
```cpp
dispvec.x() = pu[0]*runTime.deltaTValue();
iter().track(dispvec,1.0);
```

## Notes on Numerical Implementation

1. **Search Region**: The code uses a 3-cell radius for force spreading and velocity interpolation:
```cpp
scalar minX = pos.x() - 3*meshWidth;
```

2. **Time Integration**: Currently implements first-order Euler time integration for point movement.

## Future Work and Improvements

1. Implementation of higher-order time integration schemes
2. Addition of multiple Dirac delta function options
3. Parallel computation support
4. Integration with dynamic mesh handling
5. Advanced FEM coupling through the planned Fortran library
