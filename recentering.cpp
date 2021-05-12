 //--------------- recentering coordinates---------------


void recenter(const vector <float> xc_fof, const vector <float> yc_fof, const vector <float> zc_fof, const vector <float> x, const vector <float> y, const vector <float> z, vector <float> xc_rc, vector <float> yc_rc, vector <float> zc_rc, float r_max){

    int np = xc_fof.size();
    
    int ncentermin = 10; //min np for recentering
    int ncentertmp = np; //+1 for passing the first while loop below
    int nbin_rc = 10; //log(np);
    double drcbin = 0;
    double drcmax = 0;
    double r_samp_drcmax=0;
    
    double drc_crit = 0.1; //r_max/float(nbin_rc);
    
    //double drc[nbin_rc];
    //for(j = 0; j < nbin_rc; j++){drc[j] = 0.;}
    
    j = 0;
    
    ncenter = np;


    // Compute the max distance from the fof centre

    for (k = 0; k < np; k++) {//loop over particles in halo
    
        ri = sqrt(x[k]*x[k] + y[k]*y[k]+ z[k]*z[k]); // distance to halo center
    
        if(ri < r_max){
        
        r_max = ri;
        
        }
    }


    
    if(ncenter<=ncentermin){
        xc_rc = xc_fof;
        yc_rc = yc_fof;
        zc_rc = zc_fof;
    }else{
        while(ncentermin < ncentertmp && j < nbin_rc ){  // Here iterates until the max distance includes more than necentertmp particles or up to nbin_rc
            
            // Maximum radius rescaled up to which particles are consider
            // to compute the centre of mass
            r_samp =r_max*(1 - float(j)/float(nbin_rc)); // EJG: Como defino r_max
        
            xc_rc = xc_rc + xc;
            yc_rc = yc_rc + yc;
            zc_rc = zc_rc + zc;
        
            xc = 0; yc = 0; zc = 0; // EJG: Esto no deberia estar definido mas arriba
    
            ncentertmp = ncenter;
            ncenter = 0; 
        
            for (k = 0; k < np; k++) {//loop over particles in halo
            
                ri = sqrt(x[k]*x[k] + y[k]*y[k]+ z[k]*z[k]); // distance to halo center
            
                if(ri < r_samp){
                
                    xc = xc + x[k];
                    yc = yc + y[k];
                    zc = zc + z[k];
                
                    ncenter = ncenter + 1;
                
                }
            }
    
            xc = xc / double(ncenter);//center of all particles that lie within rbin
            yc = yc / double(ncenter);
            zc = zc / double(ncenter);
        
        
            if(ncentermin <=  ncenter){
                for(k = 0; k < np; k++){
    
                    x[k] = x[k] - xc;
                    y[k] = y[k] - yc;
                    z[k] = z[k] - zc;
                
                }
                //drc[j] = sqrt((xc)*(xc) + (yc)*(yc) + (zc)*(zc));
                drcbin = sqrt((xc)*(xc) + (yc)*(yc) + (zc)*(zc))/r_max;
    
                if(drcmax <= drcbin){
                    drcmax = drcbin;
                    r_samp_drcmax = r_samp;
                }
            }
        
            j++;
        }
    }
    
    //------------------------------------------------------
}
