float CalcCos(float dir1[], float dir2[]){
  float cos = ( dir1[0]*dir2[0] + dir1[1]*dir2[1] + dir1[2]*dir2[2] )
    / ( sqrt( pow(dir1[0], 2) + pow(dir1[1], 2) + pow(dir1[2], 2) )
	* sqrt( pow(dir2[0], 2) + pow(dir2[1], 2) + pow(dir2[2], 2) ) );
  return cos;
}

float Transmom(float dir[][3], int nring, float tot_energy, float amome[], int mer )
{ 
  float cos_angle[10], sin_angle[10];
  float ltransmom=0;

  for(int j=0; j< nring; j++)
    {
      cos_angle[j] = CalcCos(dir[mer], dir[j]);
      if(cos_angle[j]>=1)
	{
	  cos_angle[j]=1;
	}

      sin_angle[j]= sqrt(1-cos_angle[j]*cos_angle[j]);
      
      ltransmom = ltransmom + amome[j]*sin_angle[j];
    }
  
  return ltransmom/tot_energy;
}

float Fracmom(float tot_energy, float amome[], int mer)
{
  // fractional momentum carried by mering
  float fmom = amome[mer] / tot_energy;

  return fmom;
}


#include <iostream>
float Dpose(Int_t & nmue, float evis, float etime[], UInt_t etype[], float epos[][3], float ehit[], float pos[], float egood[], float maxamome, int skgen )
{

  //  return longest distance between primary vertex position
  //  and decay electron position, divided by mering's energy

  float dpos;
  float dx;
  float x, y , z;
  int   nmuemax;
 
  int ehit_cut;
  if( skgen == 0 ) ehit_cut = 60;
  if( skgen == 1 ) ehit_cut = 30;
  if( skgen == 2 ) ehit_cut = 60;
  if( skgen == 3 ) ehit_cut = 60;
  if( skgen == 4 ) ehit_cut = 60;

  int ehit_cut_2[5] = { 40, 20, 40, 40, 40 };

  nmuemax = ( nmue > 10 ? 10 :  nmue );
  if ( nmue <= 0 ) return 0.;

  dpos=0.;

  //  TBW 2023: nmue in h1 is standard for SK IV+
  //  Dpose calculates nmue based on events used in dpose calculation, which
  //  sometimes differs from nmue
  //  here, we have dpose compute nmue, but only use it for sk I, II, III
  if (skgen < 3) nmue = 0 ;

  for( int i = 0 ; i < nmuemax ; i++ )
    {
      if( skgen >= 3) //for SK4 only
          // note: use value of nmue passed in, do not modify
	{
	  if( etime[i] < 0.1) continue;
   
	  if( etime[i] < 0.6) continue;

	  if( etype[i] == 1 )
	    {
	      x = pos[0]-epos[i][0];
	      y = pos[1]-epos[i][1];
	      z = pos[2]-epos[i][2];
	      dx = sqrt(x*x + y*y + z*z);

              //std::cout << i << "   dx: " << dx << std::endl;
	      if( dx > dpos)
		{
		  dpos = dx;
		}
	    }
	}
      else
	{
        // sk I, II, III
        // we calculate nmue ourselves

	  if( evis < 1330.0 && etime[i] <  0.1 ) continue;
	  if( evis >= 1330.0 && etime[i] <  1.2 ) continue;
	  if( etime[i] >  0.8 && etime[i] < 1.2 ) continue;

	  if( etype[i] == 1   && ehit[i] >= ehit_cut && egood[i] > 0.5 )
	    {
              nmue++;
	      x  = pos[0] - epos[i][0];
	      y  = pos[1] - epos[i][1];
	      z  = pos[2] - epos[i][2];
	      dx = sqrt( x*x + y*y + z*z );

	      if( dx > dpos ) dpos = dx;
	    }

	  if( etype[i] >= 2  &&  etype[i] <=4 &&  ehit[i] >= ehit_cut_2[skgen] )
	    {
              // bye bye and in-gate decay electrons are counted but the vertex is
              // not reliable, exclude for dpos!
              nmue++;
            }
     
	}
    }// end of loop on decay-e's


  //return dpos / llh->GetEmax() ;
  //return dpos / maxamome ;
  return dpos;
}

void SearchMER(int nring, float amome[], float amomm[], float prmslg[][6], int *mer, float *merp){//return values are mer and merp
  *mer = 0;
  *merp = 0;
  for(int searchid = 0; searchid < nring; searchid++){
    float p = ( prmslg[searchid][1] < prmslg[searchid][2] ? amome[searchid] : amomm[searchid] );
    if(p > *merp){
      *merp = p;
      *mer = searchid;
    }
  }
}

