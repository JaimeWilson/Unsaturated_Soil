// SWCC

Vws = 0.49;
Vwr = 0.08;
Vw1 = 0.24;
Nada1 = 1.10;
Sl1 = 700;
//ZMAX = 5.25;
//mark = 0;



forAll(mesh.cells(),celli)
{
                   
     
     if (Vw[celli] < Vws)
        {
            if (Vw < Vws+0.001) 
		{
		Vw = Vwr + 0.001;  
		}
	//mark=0;
	temp1=(Vws - Vwr)/(Vw1 - Vwr)-1.0;
	temp2=(Vws - Vwr)/(Vw[celli]-Vwr)-1.0;
	S[celli] = Foam::pow(2.718282,log(temp2/temp1)/Nada1+log(Sl1));
	//Uw[celli]=-1.0*S[celli];
	//Hw[celli]=Uw[celli]/9.8+ZMAX-Z[celli];                                                                                      
        //}
     else
        { 
	   S[celli] = 0;
	   //if (mark==0)
              //{
	      //isaturate=celli;
	      //mark=1;
	      //}
	   //Vw[celli] = Vws;
	   //S[celli] = 0;	
	   //Uw[celli] = 0;
	   //for(int j = isaturate; j < celli; j++)
	      //Uw[celli] = Uw[celli] +9.8*(Z[j+1]-Z[j]);
	      //Hw[celli] = ZMAX-Z[isaturate]+Uw[celli]/9.8;
	   
	}
}


