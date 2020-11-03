// SWCC


//ZMAX = 5.25;
//mark = 0;


forAll(mesh.cells(),celli)
{
     //S[celli] = Theta[celli]*1.0;             
     
     if (Theta[celli] < Vws)
        {
            if (Theta[celli] < Vws+0.001) 
		{
		Theta[celli] = Vwr + 0.001;  
		}
	mark=0;
	temp1 = (Vws - Vwr)/(Vw1 - Vwr)-1.0;
	temp2 = (Vws - Vwr)/(Theta[celli]-Vwr)-1.0;
	S[celli] = Foam::pow(2.718282,Foam::log(temp2/temp1)/Nada1+Foam::log(Sl1));
	Uw[celli]=-1.0*S[celli];
	//Hw[celli]=Uw[celli]/9.8+ZMAX-Z[celli];                                                                                      
        }
     else
        { 
	   if (mark==0)
              {
	      isaturate = celli;
	      mark=1;
	      }
	   Theta[celli] = Vws;
	   S[celli] = 0;	
	   Uw[celli] = 0;
	   for( j = isaturate; j < celli; j++){
	      Uw[celli] = Uw[celli] + 9.8 * (mesh.C()[celli][j+1]-mesh.C()[celli][j]);
	      //Hw[celli] = ZMAX-Z[isaturate]+Uw[celli]/9.8;
	   }
	}
}


