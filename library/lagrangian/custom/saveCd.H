volScalarField& Cd = (
			this->mesh().objectRegistry::template
			lookupObjectRef<volScalarField>("Cd")
		     );

vector posP = p.position();
label cellP = this->mesh().findCell(posP);
Cd[cellP] = CdRe(Re)/Re;

