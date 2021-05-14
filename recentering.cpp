#include "recentering.h"

void recenter(const float xc_fof, const float yc_fof, const float zc_fof, vector <float> x, vector <float> y, vector <float> z, float xc_rc, float yc_rc, float zc_rc, float r_max){

    int np = x.size();
    
    int ncentermin = 10; //min np for recentering
    int ncentertmp = np; //+1 for passing the first while loop below
    int nbin_rc = 10; //log(np);
    
    double drc_crit = 0.1; //r_max/float(nbin_rc);
    
    float xc, yc, zc;
    
    //double drc[nbin_rc];
    //for(j = 0; j < nbin_rc; j++){drc[j] = 0.;}
    
    int j = 0;
    float ri;
    float r_samp;
    int ncenter;
    
    ncenter = np;

    xc = 0; yc = 0; zc = 0;

    // Compute the max distance from the fof centre

    for (int k = 0; k < np; k++) {//loop over particles in halo
    
        ri = sqrt(x[k]*x[k] + y[k]*y[k]+ z[k]*z[k]); // distance to halo center
    
        if(ri > r_max){ 
        r_max = ri;
        }
    }


    xc_rc = xc_fof;
    yc_rc = yc_fof;
    zc_rc = zc_fof;
    
    while(ncentermin < ncentertmp && j < nbin_rc ){  // Here iterates until the max distance includes more than necentertmp particles or up to nbin_rc
            
            // Maximum radius rescaled up to which particles are consider
            // to compute the centre of mass
            r_samp =r_max*(1 - float(j)/float(nbin_rc)); 
            
            ncenter = 0; 
            xc = 0; yc = 0; zc = 0; 
        
            for (int k = 0; k < np; k++) {//loop over particles in halo
            
                ri = sqrt(x[k]*x[k] + y[k]*y[k]+ z[k]*z[k]); // distance to halo center
            
                if(ri < r_samp){
                
                    xc = xc + x[k];
                    yc = yc + y[k];
                    zc = zc + z[k];
                
                    ncenter = ncenter + 1;
                
                }
            }
    
            ncentertmp = ncenter;
    
            xc = xc / double(ncenter);//center of all particles that lie within rbin
            yc = yc / double(ncenter);
            zc = zc / double(ncenter);
        
            if(ncentermin  < ncentertmp){
                
                for(int k = 0; k < np; k++){
    
                    x[k] = x[k] - xc;
                    y[k] = y[k] - yc;
                    z[k] = z[k] - zc;
                
                }

                
                xc_rc = xc_rc + xc;
                yc_rc = yc_rc + yc;
                zc_rc = zc_rc + zc;

            }

        
            j++;
        }
    
    //------------------------------------------------------
}
