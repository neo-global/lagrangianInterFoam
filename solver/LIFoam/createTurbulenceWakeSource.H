dimensionedScalar vsmall("vsmall", dimLength, VSMALL);
volScalarField SkWake = 3.0/4.0*Cd/(Db+vsmall)*rho*voidage*pow(pos(mag(slipVel)),3);

volScalarField SepsilonWake = CepsWakeDelta*SkWake*pos(mag(slipVel))/(kappa*(Db+vsmall));
//SkWake.primitiveField() = 0.0;
//SepsilonWake.primitiveField() = 0.0;

