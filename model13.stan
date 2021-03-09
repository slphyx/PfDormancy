//
// Stan code by Sompob Saralamba
// 
functions{

//concentration function
//t-time is in days Cj-concentration
real Conc(real t, real Cj){
    return(( 0<=t&&t<=0.25)?Cj:0.0);
}

//drug efficacy
//C-concentration Em-% max eff EC50-concentration at 50% eff gamma-slope of the EC curve
real Eff(real C, real Em, real EC50, real gamma){
    return(Em*pow(C,gamma)/(pow(C,gamma)+pow(EC50,gamma)));
}

//system function
vector system1(real t, vector y, real[] parms){

    vector[4] dydt;

    // parameters (MUST BE IN THIS ORDER)
    // Mw, EmW, EmR, EC50W, EC50R, gamma, rhoW, rhoR, xW,xR,f, Cj
    real Mw = parms[1];
    real EmW = parms[2];
    real EmR = parms[3];
    real EC50W = parms[4];
    real EC50R = parms[5];
    real gamma = parms[6];
    real rhoW = parms[7];
    real rhoR = parms[8];
    real xW = parms[9];
    real xR = parms[10];
    real f = parms[11];
    real Cj = parms[12];

    // system
    real W = y[1];
    real Wd = y[2];
    real R = y[3];
    real Rd = y[4];

    dydt[1] = (log(Mw)/2.0)*W - Eff(Conc(t, Cj), EmW, EC50W, gamma)*W + rhoW*Wd;
    dydt[2] = xW*Eff(Conc(t,Cj), EmW, EC50W, gamma)*W - rhoW*Wd;
    dydt[3] = (log(f*Mw)/2.0)*R - Eff(Conc(t, Cj),EmR,EC50R,gamma)*R + rhoR*Rd;
    dydt[4] = xR*Eff(Conc(t,Cj),EmR,EC50R,gamma)*R - rhoR*Rd;

    return dydt;
}

// //for converting ode output to a matrix
// matrix Vec2Mat(vector[] odeout){
//     int dimensions[2] = dims(odeout);
//     int nrow = dimensions[1];
//     int ncol = dimensions[2];

//     matrix[nrow,ncol] nmat;

//     for(i in 1:nrow){
//         nmat[i] = to_row_vector(odeout[i]);        
//     }
//     return(nmat);
// }

//to return Pij, yij 
// vector[] Py(vector[] odeout){
//     int dimensions[2] = dims(odeout);
//     int nrow=dimensions[1];
//     real Pij[nrow];  //total parasites
//     real yij[nrow];  //percentages
//     vector[2] output[nrow];

//     for(i in 1:dimensions[1]){
//         // Pij = Wij +Wdij + Rij +Rdij
//         Pij[i] = odeout[i,1]+odeout[i,2]+odeout[i,3]+odeout[i,4];
        
//         // yij = 100(Rij +Rdij)/Pij
//         yij[i] = 100.0*(odeout[i,3]+odeout[i,4])/Pij[i];
        
//     }

//     output[,1] = Pij;
//     output[,2] = yij;
//     return output;
// }

//to return Pij, yij (matrix version)
matrix Py(vector[] odeout){
    int dimensions[2] = dims(odeout);
    int nrow=dimensions[1];
    vector[nrow] Pij;  //total parasites
    vector[nrow] yij;  //percentages
    matrix[nrow,2] output;

    for(i in 1:dimensions[1]){
        // Pij = Wij +Wdij + Rij +Rdij
        Pij[i] = odeout[i,1]+odeout[i,2]+odeout[i,3]+odeout[i,4];
        
        // yij = 100(Rij +Rdij)/Pij
        yij[i] = 100.0*(odeout[i,3]+odeout[i,4])/Pij[i];
        
    }

    output[1:nrow,1] = Pij;
    output[1:nrow,2] = yij;
    return output;
}

matrix RunModel(int nts, real[] ts, real[] parms, vector inits){
    vector[4] outode[nts];
    matrix[nts,2] outPy;
    
    outode = ode_rk45(system1, inits,0,ts, parms);
    // outode = ode_bdf(system1, inits,0,ts, parms);

    outPy = Py(outode);

    return(outPy);
}


} //FUNCTION BLOCK

data {
    int nsets;  //number of data sets
    int resTsteps[nsets]; //number of res. time steps
    int totTsteps[nsets]; //number of tot. time steps
    
    matrix[15,nsets] y;     //% resistant parasites
    matrix[18,nsets] log10tps;   // log10 total parasites
 
    real t0;            // initial time
    matrix[15,nsets] resTime;    // res. output times 
    matrix[18,nsets] totTime;    // res. output times 
    real P0;            // intial parasites
    real ppi[nsets];           // proportion of the parasites
    real MCj[nsets];           // concentration
}

transformed data {
    // real x_r[0];   //for sharing real variables in ode
    // int x_i[0];    //for sharing integer variables in ode

    // vector[4] inits = [(1-ppi)*P0, 0, ppi*P0, 0]';

    //### total parasites in log scale !!!!

    matrix[nsets,4] inits;
    for(i in 1:nsets)
        inits[i,1:4] = [(1-ppi[i])*P0, 0, ppi[i]*P0, 0];
}   

parameters {
    real<lower = 0> Mw; 
    real<lower = 0> EmW;
    real<lower = 0> EmR;
    real<lower = 0> EC50W;
    real<lower = 0> EC50R;
    real<lower = 0> gamma;
    real<lower = 0> rhoW;
    real<lower = 0> rhoR;
    real<lower = 0> xW;
    real<lower = 0> xR;
    real<lower = 0> f;  
    
    real<lower = 0> sig[nsets];
    real<lower = 0> sig2[nsets];
    
}

model {
    real parms[nsets,12]; 

    matrix[15,4] res_outode[nsets];
    matrix[18,4] tot_outode[nsets];
    matrix[15,2] res_outPy[nsets]; 
    matrix[18,2] tot_outPy[nsets]; 
    
    Mw ~ uniform(1,5);  
    EmW ~ uniform(50, 500);
    EmR ~ uniform(50,500);
    EC50W ~ uniform(1,200);
    EC50R ~ uniform(1,200);
    gamma ~ uniform(0.5,5);
    rhoW ~ uniform(0,1);
    rhoR ~ uniform(0,1);
    xW ~ uniform(0,1);
    xR ~ uniform(0,1);
    f ~ uniform(0,1);
        
    sig ~ uniform(0,20);
    sig2 ~ uniform(0,4);

    for(i in 1:nsets){
        parms[i] = {Mw,EmW,EmR,EC50W,EC50R,gamma,rhoW,rhoR,xW,xR,f,MCj[i]};

        res_outPy[i,1:resTsteps[i]] = RunModel(resTsteps[i],to_array_1d(resTime[1:resTsteps[i],i]),parms[i],inits[i]');
        tot_outPy[i,1:totTsteps[i]] = RunModel(totTsteps[i],to_array_1d(totTime[1:totTsteps[i],i]),parms[i],inits[i]');

        for(j in 1:resTsteps[i]){
            y[j,i] ~ normal(res_outPy[i,j,2],sig[i]);
        }
        for(k in 1:totTsteps[i]){
            log10tps[k,i] ~ normal(log10(tot_outPy[i,k,1]),sig2[i]);
        }
        
    }
}

