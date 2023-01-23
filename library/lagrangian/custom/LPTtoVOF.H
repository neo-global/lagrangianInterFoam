/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "mathematicalConstants.H"
//#include "volScalarField.H"
using namespace Foam::constant::mathematical;

bool closeVOF= false;

vector posP = p.position();
label cellP = cloud.pMesh().findCell(posP);
scalar alpha1P = td.alpha1Interp().interpolate(posP, cellP);

scalar distanceP = d_/2.0+(pow(cloud.pMesh().cellVolumes()[cellP],1.0/3.0))/2.0;
//-RM: This assumes that closeness is at a distance of radius + half of cell and then deletes, instead of which we could use d_/2.0+p.U()*dt to calculate a point where, if there interface, we delete the particle and convert it to equivalent VOF. This might fudge the interface and needs to be checked if implemented.
//-TODO.


//-RM: collate deletion data cell-wise: no. of particles in cell.
volScalarField& nPInCell = (cloud.pMesh().objectRegistry::template lookupObjectRef<volScalarField>("nPInCell"));	//check no of particles in cell
volScalarField& surfaceCell = (cloud.pMesh().objectRegistry::template lookupObjectRef<volScalarField>("surfaceCell"));
volScalarField& cellsWithConv = (cloud.pMesh().objectRegistry::template lookupObjectRef<volScalarField>("cellsWithConv"));	//number of particles for deletion in the cell

//-RM: calculate void fraction to find alpha Source.
cloud.alpha1Su()[cellP] += p.volume()*nParticle()/cloud.pMesh().cellVolumes()[cellP]; //iterate over particles and add to cell; should be same as theta. will eventually change this to save only for closeVOF particle containing cells.

//-RM: calculate slip Velocity in a center-of-mass sense
volVectorField& slipVel = (cloud.pMesh().objectRegistry::template lookupObjectRef<volVectorField>("slipVel"));
slipVel[cellP] += U() - cloud.U()[cellP]; 		//vector addition over multiple (possible) particles in the cell.

/*
//-RM: save bubble dia as public field using d_ which is private.
volScalarField& Db = (
			cloud.pMesh().objectRegistry::template
			lookupObjectRef<volScalarField>("Db")
		     );
Db[cellP] = d_;
*/

nPInCell[cellP]++;

//-RM: Loop over all 6 directions to detect interface crossing.
for (int j=0; j<6;j++)
{
  int sign=pow(-1,j);
  int coord=j/2;

	/*
	Coords:
	coord	sign
	+x	0		 1
	-x	0		-1
	+y	1		 1
	-y	1		-1
	+z	2		 1
	-z	2		-1
	*/

  vector possurfP=posP;
  possurfP.component(coord) += sign * distanceP ;	//pos of surf = posP + (d+a)/2
  label cellsurfP=cloud.pMesh().findCell(possurfP);

  if (cellsurfP != (-1))
  {
    scalar alpha1surfP = td.alpha1Interp().interpolate(possurfP, cellsurfP);
    if (alpha1surfP > 0.9)
    {
	closeVOF= true;
	if(cellsurfP != cellP)
	{
	    surfaceCell[cellP] = cellsurfP;	//increment if cell different.
	}
      //Does not trigger for parallel, non-conforming meshes.
    }
  }
}

if (closeVOF)
{
    td.keepParticle = false;
    cellsWithConv[cellP] += closeVOF;

    scalar VParticle = p.volume()*p.nParticle();
    scalar VCell = cloud.pMesh().cellVolumes()[cellP];	//Volume of cell

    if (p.volume()*p.nParticle() <= cloud.pMesh().cellVolumes()[cellP])
    // fill only the eulerian cell at the particle position
    //Should be non-problematic as on coarse grid, dParticle < cellSize by design
    {
        cloud.alpha1Su()[cellP] += VParticle/VCell;
        cloud.alpha1Su()[cellP] += p.volume()*p.nParticle()/cloud.pMesh().cellVolumes()[cellP];

        if(cloud.alpha1Su()[cellP]+cloud.alpha1()[cellP] > 1.0)
        {
          Info<<"total lpt and cnt volume fracton > 1.0"<<endl;
        }
        else
        {

        }
   // td.cloud().correctU()[cellP] += U_  ;
   // Info<< "CellP has alpha= " <<alpha1P << " substract(add) liquid(gas) to cellP: "
   // << cellP << "  "<< td.cloud().correctalpha1()[cellP] << "  "
   // << td.cloud().correctU()[cellP] <<endl;
    }
    else

  // start to fill the neighbour cells
    {
      Info<<"Vparticle > Vcell"<<endl;
      //Start with basics.
    }
}
