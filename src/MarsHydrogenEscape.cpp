#include "MarsHydrogenEscape.h"

double steam_from_impact(double& N_H2O_ocean, double& N_CO2, 
                      double& N_N2, double& m_i, double& area)
{
  const double v_i = 1700000.0;
  const double impactor_energy_frac = 0.5;
  const double T_prime = 2000;
  const double T_0 = 300; // inital temperature
  const double C_H2O = 18649674.13020014;
  const double C_w = C_H2O;
  const double C_CO2 = 8456820.456951989;
  const double C_N2 = 10290056.559022075;
  // const double C_CO = 10405216.637589749;
  const double Q_w = 2.5e10; // latent heat of water vaporization ergs/g;

  // convert to grams from mol/cm2
  double M_H2O = N_H2O_ocean*18.015*area;
  double M_CO2 = N_CO2*44.009*area;
  double M_N2 = N_N2*28.014*area;

  // heat capacity of dry atmosphere (erg/K)
  double M_AC_A = M_CO2*C_CO2 + M_N2*C_N2;

  // energy of impactor (ergs)
  double E_i = 0.5*m_i*v_i*v_i;

  // assume that the heated atm/steam is heated to T_prime
  double M_steam = (E_i*impactor_energy_frac-M_AC_A*(T_prime-T_0))/(Q_w + (T_prime-T_0)*C_w);

  // Can't be more steam than water
  if (M_steam > M_H2O){ 
    M_steam = M_H2O;
  }

  // convert to moles/cm2 of steam
  double N_H2O_steam = M_steam/18.015/area;
  return N_H2O_steam;
}


void react_iron(double& N_H2O_steam,double& N_CO2,double& M_i, double& area, double& Fe_react_frac,
                double& N_H2, double& N_H2O, double& N_CO, double& N_CO2_out)
{
  const double impactor_Fe_frac = 0.33;
  double Moles_Fe = Fe_react_frac*impactor_Fe_frac*M_i/56.;
  double N_Fe = Moles_Fe/area;
  double xxx = N_Fe/(N_H2O_steam +N_CO2); //moles Fe/moles O
  if (xxx > 1.0){
    xxx = 1.0;
  }
  N_H2 = xxx*N_H2O_steam;
  N_H2O = (1-xxx)*N_H2O_steam;
  N_CO = xxx*N_CO2;
  N_CO2_out = (1-xxx)*N_CO2;
  // }
  // else{
  //   throw "More Fe than O!";
  // }   
}

double gel2N(double& gel){
  double mu = 18; // g/mol H2O
  double rho_water = 1; // g/cm3 density of water
  return gel * rho_water / mu; // moles/cm2
}

double bars2N(double& P_bar, double mu, double& g){
  double P_cgi = P_bar * 1.0e6; // convert bars to dynes/cm2 (the cgi unit system)
  double N = P_cgi/(mu * g);
  return N;
}

// km
double mass(double& D){
  const double rho = 3.0e12; // kg/km3 (Morbidelli et al. 2012)
  return rho*(4./3.)*3.14159*std::pow(D/2.0,3.0)*1.0e3; // mass in g
}

void impact(double& ocean_GEL_cm, double& P_CO2_bars, double& P_N2_bars, 
            double& Fe_react_frac, double& M_i, double& grav, double& area,
            double& N_H2, double& N_t, double& N_x, double& mu_x)
{
  double N_H2O_ocean = gel2N(ocean_GEL_cm);
  double N_CO2 = bars2N(P_CO2_bars, 44.0, grav);
  double N_N2 = bars2N(P_N2_bars, 28.0, grav);
  double N_H2O_steam = steam_from_impact(N_H2O_ocean, N_CO2, N_N2, M_i, area);
  // React the atmosphere with impactor Fe
  double N_H2O, N_CO, N_CO2_out;
  react_iron(N_H2O_steam, N_CO2, M_i, area, Fe_react_frac,
             N_H2, N_H2O, N_CO, N_CO2_out);

  // total is
  N_t = N_H2 + N_CO + N_CO2_out + N_N2;
  N_x = N_CO + N_CO2_out + N_N2;
  mu_x = N_CO/N_x*28. + N_CO2_out/N_x*44. + N_N2/N_x*28.;

}

MarsClimate::MarsClimate(){
  
  std::vector<double> Psurf;
  std::vector<double> h2;
  std::vector<double> Tsurf;

  // double data;
  // std::ifstream infile; 
  // std::vector<std::string> filenames(4);
  // filenames[0] = "data/0.5.csv";
  // filenames[1] = "data/1.0.csv";
  // filenames[2] = "data/1.5.csv";
  // filenames[3] = "data/2.0.csv";
  // std::vector<double> pressures = {0.5, 1.0, 1.5, 2.0};
  
  // for (int i = 0; i < filenames.size(); i++){
  //   infile.open(filenames[i]);
  //   while ( infile.peek() != EOF ){
  //     Psurf.push_back(pressures[i]);
  //     infile >> data; 
  //     h2.push_back(data);
  //     infile >> data; 
  //     Tsurf.push_back(data);
  //   }
  //   infile.close();
  // }
  
  h2 = {-0.0008847827149630894,0.016025810236105     ,0.04033729285819035   ,
        0.06888302292613775   ,0.09390917534340323   ,0.09955027618907954   ,
        0.1                   ,-0.0009277997751380891,0.012754558341887861  ,
        0.027515275944664425  ,0.046162193870068924  ,0.07291782763846116   ,
        0.09581854621889817   ,0.09933910153003862   ,0.1                   ,
        0.00009678838539375553,0.007077284059246226  ,0.016873441853644238  ,
        0.030918512000782137  ,0.04287236642713986   ,0.05659774160434082   ,
        0.06926919880725424   ,0.08441511463068876   ,0.09921493865180624   ,
        0.1                   ,-0.0002776555702204593,0.004931319352788777  ,
        0.009445177689788347  ,0.01609620178911864   ,0.028367795864496264  ,
        0.03750696583076697   ,0.04980984504081733   ,0.07023121669844064   ,
        0.09313975656254583   ,0.0987798797477636    ,0.1                   };
  
  Psurf = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,
           1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
           1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,
           2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0};
           
  Tsurf = {218.89117009662544,223.2238353619788 ,229.09713056655426,
           235.11870427400564,239.74711182806215,240.59115869058678,
           240.59115869058678,224.9852536214173 ,236.65118704274005,
           245.5495266493947 ,253.90282055042286,263.5213700281892 ,
           269.2529044662789 ,270.5075687213831 ,270.5075687213831 ,
           229.83526421273893,240.93171041697218,253.1426895439214 ,
           263.42441869938574,269.9617082987079 ,275.5335581952388 ,
           280.41045444916983,284.7390461292793 ,288.097309804305  ,
           288.097309804305  ,232.88149125808607,244.94337716511058,
           255.48011275683956,263.2516986850467 ,274.7758713398837 ,
           280.06012611819915,287.152238027733  ,294.12458653109775,
           298.7481057828616 ,299.730654543677  ,300.0             ,};
    
  nd = h2.size();
  triangle = new int[3*2*nd];
  element_neighbor = new int[3*2*nd];
  xyd = new double[2*nd];
  zd = new double[nd];
  
  for (int i = 0; i < nd; i++ )
  {
    xyd[0+i*2] = h2[i];
    xyd[1+i*2] = Psurf[i];
    zd[i] = Tsurf[i];
  }
  
  r8tris2(nd, xyd, element_num, triangle, element_neighbor);
  for (int j = 0; j < element_num; j++ )
  {
    for (int i = 0; i < 3; i++ )
    {
      if ( 0 < element_neighbor[i+j*3] )
      {
        element_neighbor[i+j*3] = element_neighbor[i+j*3] - 1;
      }
    }
  }  
}

MarsClimate::~MarsClimate(){
  delete[] triangle;
  delete[] element_neighbor;
  delete[] xyd;
  delete[] zd;
}

double MarsClimate::MarsSurfTemp(double h2, double Psurf){
    
  const int ni = 1;
  double xyi[2];
  
  if (h2 < 0.0){
    h2 = 0.0;
  }
  if (h2 > 0.1){
    h2 = 0.1;
  }
  if (Psurf < 0.5){
    Psurf = 0.5;
  }
  if (Psurf > 2.0){
    Psurf = 2.0;
  }
  
  xyi[0] = h2;
  xyi[1] = Psurf;
  
  double *zi;
  zi = pwl_interp_2d_scattered_value(nd, xyd, zd, element_num, 
                                     triangle, element_neighbor, ni, xyi );
  double SurfTemp = zi[0];
  delete[] zi;
  return SurfTemp;
}

  
void MarsHydrogenEscape::rhs(double t, double *y, double *ydot, void *data){
  const double A_escape = 2e12;
  const double B_escape = 0.006;
  
  double *data_{static_cast<double*>(data)}; 
  
  double N_avo_ = data_[0];
  double S_1_ = data_[1];
  double N_x_ = data_[2];
  
  double N_H2 = y[0];
  double N_t = N_x_ + N_H2;
  double fH2 = N_H2/N_t;
  
  double dNH2_dt = -(A_escape*S_1_)/(std::sqrt(1.0 + B_escape*S_1_*S_1_)) * fH2;
  ydot[0] = dNH2_dt/N_avo_;
}
  
double MarsHydrogenEscape::WarmTimeAfterImpact(double& ocean_GEL_cm, double& P_CO2_bars, double& P_N2_bars, \
                                              double& Fe_react_frac, double& M_i){
    
  double N_H2_init, N_t_init, N_x, mu_x;
  impact(ocean_GEL_cm, P_CO2_bars, P_N2_bars, 
         Fe_react_frac, M_i, grav, area,
         N_H2_init, N_t_init, N_x, mu_x);
  
  int neq = 1;
  double t, tout;
  t = 0e0;
  tout = 1.0;
  int istate = 1;
  
  std::vector<double> y = {N_H2_init};
  std::vector<double> yout;
  double data[] = {N_avo, S_1, N_x};
  
  lsoda.lsoda_update(rhs, neq, y, yout, &t, tout, &istate, &data, rtol, atol);
  double N_t = N_x + yout[1];
  double fH2 = yout[1]/N_t;
  double mubar = fH2*2.0 + N_x/N_t*mu_x;
  double Psurf = N_t*mubar*grav/1.0e6;
  double Temp = mars.MarsSurfTemp(fH2,Psurf);
  
  if (Temp < 273.0){
    return 0.0;
  }
  
  while (true){
    lsoda.stoda(neq, yout, rhs, data);
    N_t = N_x + yout[1];
    fH2 = yout[1]/N_t;
    mubar = fH2*2.0 + N_x/N_t*mu_x;
    Psurf = N_t*mubar*grav/1.0e6;
    Temp = mars.MarsSurfTemp(fH2,Psurf);
    // std::cout << fH2 << " " << Psurf << " " << Temp << " " << lsoda.tn_ << std::endl;
    if (Temp < 273.0){
      break;
    }
  }
  return lsoda.tn_;
}
  
double MarsHydrogenEscape::WarmTimeAfterSeveralImpacts(std::vector<double>& x){
  double ocean_GEL_cm = x[0];
  double P_CO2_bars = x[1];
  double P_N2_bars = x[2];
  double Fe_react_frac = x[3];
  
  double tot_warm_time = 0.0;
  for (int i = 0; i < impact_diameters.size(); i++){
    double M_i = mass(impact_diameters[i]);
    tot_warm_time += WarmTimeAfterImpact(ocean_GEL_cm, P_CO2_bars, P_N2_bars, \
                                         Fe_react_frac, M_i);
  }
  return tot_warm_time;
}


// pybind11
double MarsHydrogenEscape::get_S_1(){
  return S_1;
}
void MarsHydrogenEscape::set_S_1(double value){
  S_1 = value;
}

double MarsHydrogenEscape::get_grav(){
  return grav;
}
void MarsHydrogenEscape::set_grav(double value){
  grav = value;
}

double MarsHydrogenEscape::get_area(){
  return area;
}
void MarsHydrogenEscape::set_area(double value){
  area = value;
}

double MarsHydrogenEscape::get_rtol(){
  return rtol;
}
void MarsHydrogenEscape::set_rtol(double value){
  rtol = value;
}

double MarsHydrogenEscape::get_atol(){
  return atol;
}
void MarsHydrogenEscape::set_atol(double value){
  atol = value;
}

std::vector<double> MarsHydrogenEscape::get_impact_diameters(){
  return impact_diameters;
}
void MarsHydrogenEscape::set_impact_diameters(std::vector<double> arr){
  impact_diameters = arr;
}



