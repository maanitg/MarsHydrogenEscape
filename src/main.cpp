#include "MarsHydrogenEscape.h"

int main(){
  double ocean_GEL_cm = 200e2;
  double P_CO2_bars = 1;
  double P_N2_bars = 1;
  double Fe_react_frac = 1;
  double M_i = 2.0e23;
  
  std::vector<double> x(4);
  x[0] = ocean_GEL_cm;
  x[1] = P_CO2_bars;
  x[2] = P_N2_bars;
  x[3] = Fe_react_frac;

  MarsHydrogenEscape Hesc; 
  double WarmTime;
  WarmTime = Hesc.WarmTimeAfterSeveralImpacts(x);
  
  std::cout<< WarmTime << std::endl;

}