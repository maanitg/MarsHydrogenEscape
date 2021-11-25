#pragma once
#include <vector>

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "libsoda/LSODA.h"
#include "interpolation/pwl_interp_2d_scattered.hpp"
#include "interpolation/r8lib.hpp"

double steam_from_impact(double& N_H2O_ocean, double& N_CO2, 
                      double& N_N2, double& m_i, double& area);

double react_iron(double& N_H2O_steam,double& N_CO2,double& M_i, double& area, double& Fe_react_frac,
                double& N_H2, double& N_H2O, double& N_CO, double& N_CO2_out);

double gel2N(double& gel);

double bars2N(double& P_bar, double mu, double& g);

double mass(double& D);

double impact(double& ocean_GEL_cm, double& P_CO2_bars, double& P_N2_bars, 
            double& Fe_react_frac, double& M_i, double& grav, double& area,
            double& N_H2, double& N_t, double& N_x, double& mu_x);

class MarsClimate
{
private:
  int nd;
  int *triangle;
  int *element_neighbor;
  int element_num;
  
  double *xyd;
  double *zd;
  
public:  
  MarsClimate();
  ~MarsClimate();
  double MarsSurfTemp(double h2, double Psurf);
};

class MarsHydrogenEscape
{
public:
  const double N_avo = 6.022e23;
  
  double S_1 = 20.0*0.4305564166683889;
  double grav = 372.1; 
  double area = 1.448e18;
  
  double rtol = 1.0e-6;
  double atol = 1.0e-100;
  
  std::vector<double> impact_diameters = {108.8, 134.2, 134.9,
                                          151.7, 163.3, 168.4,
                                          172.8, 173.5, 180.6,
                                          216.4, 216.8, 225.5,
                                          274.4, 300.2, 356.4,
                                          378.4, 465.1};
  MarsClimate mars;
  LSODA lsoda;
  
  static void rhs(double t, double *y, double *ydot, void *data);
  
  std::vector<double> WarmTimeAfterImpact(double& ocean_GEL_cm, double& P_CO2_bars, double& P_N2_bars, \
                             double& Fe_react_frac, double& M_i);
  
  double WarmTimeAfterSeveralImpacts(std::vector<double>& x);
  
  // pybind11
  double get_S_1();
  void set_S_1(double value);
  
  double get_grav();
  void set_grav(double value);
  
  double get_area();
  void set_area(double value);
  
  double get_rtol();
  void set_rtol(double value);
  
  double get_atol();
  void set_atol(double value);
  
  std::vector<double> get_impact_diameters();
  void set_impact_diameters(std::vector<double> value);

};
