//Nada2 = 1.25;
//Sl2 = 40;
//Kw1 = 0.000000000000012;
//Ks = 0.000000012;

forAll(mesh.cells(),celli)
{	
	k[celli] = Ks/(1+(Ks/Kw1-1.0)*Foam::pow(S[celli]/Sl2,Nada2));

	if(k[celli] < 0.000000000011) k[celli] = 0.000000000011;
}
