/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pointSolver

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cloud.H"
#include "passiveParticle.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "diracdelta.H"
#include "cellSet.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Read the mesh
    #include "createMesh.H"
   
    #include "createFields.H"
    #include "createIbPoints.H"
    
    Info<< nl;
    runTime.printExecutionTime(Info);

    // Get neighbouring cells around a point
    // Define the bounding box
    scalar minX = 0.025;
    scalar minY = 0.025;
    scalar minZ = 0.0;
    scalar maxX = 0.075;
    scalar maxY = 0.075;
    scalar maxZ = 0.01;

    treeBoundBox searchBox(point(minX, minY, minZ), point(maxX, maxY, maxZ));

    // Access the octree from the mesh
    const indexedOctree<treeDataCell>& cellTree = mesh.cellTree();

    // Query the cells in the bounding box
    labelHashSet foundCells;
    label count = cellTree.findBox(searchBox, foundCells);


    // Spread the force at that point
    vector pf(0.025, 0.0, 0.00);


   while (runTime.loop())
   {
        #include "UEqn.H"
        F = F*0;
        forAllConstIter(labelHashSet, foundCells, iter){
            label icell = iter.key();
            vector rr(0.0, 0.0, 0.0);
            rr[0] = pos1.x();
            rr[1] = pos1.y();
            rr[2] = pos1.z();
            #include "interpolateForces.H"
        }
        
        runTime.write();
   }

   Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
