/**
 * @file main.cpp
 * @author DGC@lide
 * @brief  test of the Kalman uNavIns with Arduino data
 * @version
 * @date 2019-08
 * 
 * @copyright Copyright © LIDE 2019
 * 
 */

#include <cmath>
#include <iostream>

#include "uNavINS.h"

int main(int argc, char* argv[])
{
    // Input file
    char filename[512] = "";
    sprintf(filename,argv[1]);//"../Trip_612525837_.txt");
    FILE* readfile = fopen(filename,"r");
    if(readfile == NULL){
        printf("ERROR: unable to open read file %s\n",filename);
        exit(1);
    }

    // Output file
    char filename2[512] = "";
    sprintf(filename2,argv[2]); //"../Trip_612525837_filtered.txt");
    FILE* exportfile = fopen(filename2,"w+");
    if(exportfile == NULL){
        printf("ERROR: unable to open %s\n",filename2);
        exit(1);
    }

    //Constants
    float R_e = 6378137.0;

    //Data to be read
    
    //Arduino's
    double t;
    int nsat;
    long lon_i, lat_i, alt_i, v_i;
    float ax, ay, az, p,q,r; 
    float vn,ve,vd,vnorm;
 
    //Data to be exported
    float phi,psi;

    uNavINS filter;

    //Working variables
    char buff[512];
    double t_prev=0.;
    double lon,lat,alt,v;
    double lon_prev=0., lat_prev=0., alt_prev=0.,theta, bgx,bgy,bgz,bax,bay,baz;
    bool init = true;

    printf("Start reading %s and exporting to %s...\n",filename, filename2);

    //Reading and working
    while (fscanf(readfile,"%s,", buff) != EOF) {

        sscanf(buff,"%lf,%li,%li,%li,%li,%i,%f,%f,%f,%f,%f,%f,", &t, &lat_i, &lon_i, &alt_i, &v_i, &nsat, &ax, &ay, &az, &p,&q,&r  );
        
        // printf("%s\n",buff);
        //printf("%f,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,\n\n", t, lon_i, lat_i, alt_i, v_i, nsat, ax, ay, az, p,q,r  );

        lat = lat_i*1e-7 *M_PI/180.; 
        lon = lon_i*1e-7 *M_PI/180.; 
        alt = alt_i/100.0;
        v   = v_i*1.852 /1000.; //mkn 2 m/s

    

        vn = (lat-lat_prev)/(t-t_prev)*R_e;
        ve = (lon-lon_prev)/(t-t_prev)*R_e *sin(lat);
        vd = (alt-alt_prev)/(t-t_prev);
        vnorm = sqrt(ve*ve+vn*vn+vd*vd);
        vn /= vnorm;
        ve /= vnorm;
        vd /= vnorm;

        if (init){
            vn   = 0.;
            ve   = 0. ;
            vd   = 0. ;
        }
        init = false;

        t_prev   = t;
        lat_prev = lat;
        lon_prev = lon;
        alt_prev = alt;


        //Perform filter operation:

#ifdef DEBUG        
        printf("ARGS: %12.10lg,%12.4e,%12.4e,%12.4e,%12.10e,%12.10e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e\n", \
                            t, vn,ve,vd, lon,lat,alt, ax, ay, az, p,q,r);
#endif
        filter.update(t, vn, ve, vd, lat, lon, alt, p, q, r, ax, ay, az, 0.0,1.0,0.0);

        theta = filter.getPitch_rad();
        phi   = filter.getRoll_rad();
        //psi = filter.getYaw_rad();
        psi   = filter.getHeading_rad();
        lat   = filter.getLatitude_rad();
        lon   = filter.getLongitude_rad();
        alt   = filter.getAltitude_m();
        vn    = filter.getVelNorth_ms();
        ve    = filter.getVelEast_ms();
        vd    = filter.getVelDown_ms();
        // trck  = filter.getGroundTrack_rad();
        bgx   = filter.getGyroBiasX_rads();
        bgy   = filter.getGyroBiasY_rads();
        bgz   = filter.getGyroBiasZ_rads();
        bax   = filter.getAccelBiasX_mss();
        bay   = filter.getAccelBiasY_mss();
        baz   = filter.getAccelBiasZ_mss();

        // fprintf(exportfile,"%12.12e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e\n", \
        //                     t, lon, lat, alt, v,vnorm, vn,ve,vd, ax, ay, az, p,q,r);

        fprintf(exportfile,"%12.12e,%12.6e,%12.6e,%12.6e,%12.12e,%12.12e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e,%12.6e\n", \
                            t, theta, phi, psi, lat, lon, alt, vn, ve, vd, bgx  ,bgy  , bgz  ,bax  ,bay  , baz  );
    }   
        
    fclose(readfile);
    fclose(exportfile);

    return 0;
}