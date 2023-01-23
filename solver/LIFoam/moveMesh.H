// Get the patch ID where impingement occurs: patch name must be "wall"
    label patchWallID = mesh.boundaryMesh().findPatchID("wall");
    const fvPatch& patchWallFaces = mesh.boundary()[patchWallID];

    Info<< "#--------------- Displacement ---------------# " << endl;

//Find the reference to the location of pointDisplacement field
pointVectorField& PointDisplacement = const_cast<pointVectorField&>
(
    mesh.objectRegistry::lookupObject<pointVectorField>
    (
    "pointDisplacement"
    )
);

//Get the vector field of the patch
vectorField &pDisp=refCast<vectorField>(PointDisplacement.boundaryField()[patchWallID]);

//Find the relevant size of the vector and declare a vectorField.
int Psize= pDisp.size();
vectorField dispVals(Psize);

// ########################################################################## //

//- set-up interpolator
primitivePatchInterpolation patchInterpolator (mesh.boundaryMesh()[patchWallID] );

scalarField MimpPatch = Mimp.boundaryField()[patchWallID];

//- perform interpolation
scalarField faceValues = patchInterpolator.faceToPointInterpolate(MimpPatch);

vectorField &PointPointer = refCast<vectorField>(PointDisplacement.boundaryField()[patchWallID]);
vectorField PointNormalVector = mesh.boundaryMesh()[patchWallID].pointNormals();

// loop over points to move the nodes
forAll(dispVals, index)
{
    if(faceValues[index] > 0.001)
    {
    dispVals[index].x() = PointPointer[index].x() 
     - faceValues[index] * PointNormalVector[index].x() * runTime.deltaT().value() / rho_mat.value();
    dispVals[index].y() = PointPointer[index].y() 
     - faceValues[index] * PointNormalVector[index].y() * runTime.deltaT().value() / rho_mat.value();

        if(Problem2D.value() == 1)
        {
        dispVals[index].z() = 0.0;
        }
        else
        {
        dispVals[index].z() = PointPointer[index].z() 
         - faceValues[index] * PointNormalVector[index].z() * runTime.deltaT().value() / rho_mat.value();
        }
    }
    else
    {
    dispVals[index].x() = PointPointer[index].x();
    dispVals[index].y() = PointPointer[index].y();

        if(Problem2D.value() == 1)
        {
        dispVals[index].z() = 0.0;
        }
        else
        {
        dispVals[index].z() = PointPointer[index].z();
        }
    }
}

// ########################################################################## //

//Once the values have been assigned to dispVals, assign them to cellDisplacement boundaryField
PointDisplacement.boundaryField()[patchWallID] == dispVals;

mesh.update();
mesh.moving(false); //set moving = false to allow localEuler ddt - only for negligible displacement compared to the cell size
        
// ########################################################################## //