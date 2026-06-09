// Blacklight simulation reader - HDF5 general structure interface

// C++ headers
#include <cstring>  // memcpy
#include <fstream>  // ifstream
#include <ios>      // streamoff
#include <string>   // string

// Blacklight headers
#include "mc_reader.hpp"
#include "../blacklight.hpp"        // enums
#include "../utils/array.hpp"       // Array
#include "../utils/exceptions.hpp"  // BlacklightException
//--------------------------------------------------------------------------------------------------

// Function for calculating covariant metric components in simulation coordinates
// Inputs:
//   x, y, z: Cartesian Kerr-Schild coordinates
// Outputs:
//   gcov: components set
// Notes:
//   Assumes gcov is allocated to be 4*4.
void MCReader::CovariantSimulationMetric(double x, double y, double z, double gcov[4][4])
    const
{
  // Calculate Cartesian Kerr-Schild metric
  if (simulation_coord == Coordinates::cks)
  {
    // Calculate useful quantities
    double a2 = bh_a * bh_a;
    double rr2 = x * x + y * y + z * z;
    double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
    double r = std::sqrt(r2);
    double f = 2.0 * bh_m * r2 * r / (r2 * r2 + a2 * z * z);

    // Calculate null vector
    double l_0 = 1.0;
    double l_1 = (r * x + bh_a * y) / (r2 + a2);
    double l_2 = (r * y - bh_a * x) / (r2 + a2);
    double l_3 = z / r;

    // Calculate metric components
    gcov[0][0] = f * l_0 * l_0 - 1.0;
    gcov[0][1] = f * l_0 * l_1;
    gcov[0][2] = f * l_0 * l_2;
    gcov[0][3] = f * l_0 * l_3;
    gcov[1][0] = f * l_1 * l_0;
    gcov[1][1] = f * l_1 * l_1 + 1.0;
    gcov[1][2] = f * l_1 * l_2;
    gcov[1][3] = f * l_1 * l_3;
    gcov[2][0] = f * l_2 * l_0;
    gcov[2][1] = f * l_2 * l_1;
    gcov[2][2] = f * l_2 * l_2 + 1.0;
    gcov[2][3] = f * l_2 * l_3;
    gcov[3][0] = f * l_3 * l_0;
    gcov[3][1] = f * l_3 * l_1;
    gcov[3][2] = f * l_3 * l_2;
    gcov[3][3] = f * l_3 * l_3 + 1.0;
  }

  // Calculate spherical Kerr-Schild metric
  else if (simulation_coord == Coordinates::sks or simulation_coord == Coordinates::fmks)
  {
    // Calculate useful quantities
    double a2 = bh_a * bh_a;
    double rr2 = x * x + y * y + z * z;
    double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
    double r = std::sqrt(r2);
    double cth = z / r;
    double cth2 = cth * cth;
    double sth2 = 1.0 - cth2;
    double sigma = r2 + a2 * cth2;

    // Calculate metric components
    gcov[0][0] = -(1.0 - 2.0 * bh_m * r / sigma);
    gcov[0][1] = 2.0 * bh_m * r / sigma;
    gcov[0][2] = 0.0;
    gcov[0][3] = -2.0 * bh_m * bh_a * r * sth2 / sigma;
    gcov[1][0] = 2.0 * bh_m * r / sigma;
    gcov[1][1] = 1.0 + 2.0 * bh_m * r / sigma;
    gcov[1][2] = 0.0;
    gcov[1][3] = -(1.0 + 2.0 * bh_m * r / sigma) * bh_a * sth2;
    gcov[2][0] = 0.0;
    gcov[2][1] = 0.0;
    gcov[2][2] = sigma;
    gcov[2][3] = 0.0;
    gcov[3][0] = -2.0 * bh_m * bh_a * r * sth2 / sigma;
    gcov[3][1] = -(1.0 + 2.0 * bh_m * r / sigma) * bh_a * sth2;
    gcov[3][2] = 0.0;
    gcov[3][3] = (r2 + a2 + 2.0 * bh_m * a2 * r * sth2 / sigma) * sth2;
  } else if(simulation_coord == Coordinates::spm){
    //TEGAN: put here the covariant metric for SPM coordinates 
    double rr2 = x * x + y * y + z * z;
    double r = std::sqrt(rr2);
    double cth = z / r;
    double cth2 = cth * cth;
    double sth2 = 1.0 - cth2;

    gcov[0][0] = -1.0;
    gcov[0][1] = 0.0;
    gcov[0][2] = 0.0;
    gcov[0][3] = 0.0;
    gcov[1][0] = 0.0;
    gcov[1][1] = 1.0;
    gcov[1][2] = 0.0;
    gcov[1][3] = 0.0;
    gcov[2][0] = 0.0;
    gcov[2][1] = 0.0;
    gcov[2][2] = rr2;
    gcov[2][3] = 0.0;
    gcov[3][0] = 0.0;
    gcov[3][1] = 0.0;
    gcov[3][2] = 0.0;
    gcov[3][3] = rr2*sth2;
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating contravariant metric components in simulation coordinates
// Inputs:
//   x, y, z: Cartesian Kerr-Schild coordinates
// Outputs:
//   gcon: components set
// Notes:
//   Assumes gcon is allocated to be 4*4.
void MCReader::ContravariantSimulationMetric(double x, double y, double z,double gcon[4][4]) 
    const
{
  // Calculate Cartesian Kerr-Schild metric
  if (simulation_coord == Coordinates::cks)
  {
    // Calculate useful quantities
    double a2 = bh_a * bh_a;
    double rr2 = x * x + y * y + z * z;
    double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
    double r = std::sqrt(r2);
    double f = 2.0 * bh_m * r2 * r / (r2 * r2 + a2 * z * z);

    // Calculate null vector
    double l0 = -1.0;
    double l1 = (r * x + bh_a * y) / (r2 + a2);
    double l2 = (r * y - bh_a * x) / (r2 + a2);
    double l3 = z / r;

    // Calculate metric components
    gcon[0][0] = -f * l0 * l0 - 1.0;
    gcon[0][1] = -f * l0 * l1;
    gcon[0][2] = -f * l0 * l2;
    gcon[0][3] = -f * l0 * l3;
    gcon[1][0] = -f * l1 * l0;
    gcon[1][1] = -f * l1 * l1 + 1.0;
    gcon[1][2] = -f * l1 * l2;
    gcon[1][3] = -f * l1 * l3;
    gcon[2][0] = -f * l2 * l0;
    gcon[2][1] = -f * l2 * l1;
    gcon[2][2] = -f * l2 * l2 + 1.0;
    gcon[2][3] = -f * l2 * l3;
    gcon[3][0] = -f * l3 * l0;
    gcon[3][1] = -f * l3 * l1;
    gcon[3][2] = -f * l3 * l2;
    gcon[3][3] = -f * l3 * l3 + 1.0;
  }

  // Calculate spherical Kerr-Schild metric
  else if (simulation_coord == Coordinates::sks or simulation_coord == Coordinates::fmks)
  {
    // Calculate useful quantities
    double a2 = bh_a * bh_a;
    double rr2 = x * x + y * y + z * z;
    double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
    double r = std::sqrt(r2);
    double cth = z / r;
    double cth2 = cth * cth;
    double sth2 = 1.0 - cth2;
    double delta = r2 - 2.0 * bh_m * r + a2;
    double sigma = r2 + a2 * cth2;

    // Calculate metric components
    gcon[0][0] = -(1.0 + 2.0 * bh_m * r / sigma);
    gcon[0][1] = 2.0 * bh_m * r / sigma;
    gcon[0][2] = 0.0;
    gcon[0][3] = 0.0;
    gcon[1][0] = 2.0 * bh_m * r / sigma;
    gcon[1][1] = delta / sigma;
    gcon[1][2] = 0.0;
    gcon[1][3] = bh_a / sigma;
    gcon[2][0] = 0.0;
    gcon[2][1] = 0.0;
    gcon[2][2] = 1.0 / sigma;
    gcon[2][3] = 0.0;
    gcon[3][0] = 0.0;
    gcon[3][1] = bh_a / sigma;
    gcon[3][2] = 0.0;
    gcon[3][3] = 1.0 / (sigma * sth2);
  }

  else if(simulation_coord == Coordinates::spm){
    //TEGAN: put here the contravariant metric for SPM coordinates 
    double rr2 = x * x + y * y + z * z;
    double r = std::sqrt(rr2);
    double cth = z / r;
    double cth2 = cth * cth;
    double sth2 = 1.0 - cth2;

    gcon[0][0] = -1.0;
    gcon[0][1] = 0.0;
    gcon[0][2] = 0.0;
    gcon[0][3] = 0.0;
    gcon[1][0] = 0.0;
    gcon[1][1] = 1.0;
    gcon[1][2] = 0.0;
    gcon[1][3] = 0.0;
    gcon[2][0] = 0.0;
    gcon[2][1] = 0.0;
    gcon[2][2] = 1.0/(rr2);
    gcon[2][3] = 0.0;
    gcon[3][0] = 0.0;
    gcon[3][1] = 0.0;
    gcon[3][2] = 0.0;
    gcon[3][3] = 1.0/(rr2*sth2);
  }
  return;
}
