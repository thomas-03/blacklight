// Blacklight radiation integrator - formula radiative transfer coefficients

// C++ headers
#include <cmath>   // abs, acos, atan, atan2, cos, exp, pow, sin, sqrt
#include <limits>  // numeric_limits

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"         // Math
#include "../utils/array.hpp"        // Array

//--------------------------------------------------------------------------------------------------

// Function for integrating radiative transfer equation based on formula
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes sample_flags[adaptive_level], sample_num[adaptive_level], sample_pos[adaptive_level],
//       sample_dir[adaptive_level], and momentum_factors[adaptive_level] have been set.
//   Allocates and initializes j_i[adaptive_level] and alpha_i[adaptive_level].
//   References code comparison paper 2020 ApJ 897 148 (C).
void RadiationIntegrator::CalculateFormulaCoefficients()
{
  // Allocate arrays
  int num_pix = camera_num_pix;
  if (adaptive_level > 0)
    num_pix = block_counts[adaptive_level] * block_num_pix;
  if (first_time or adaptive_level > 0)
  {
    j_i[adaptive_level].Allocate(image_num_frequencies, num_pix,
        geodesic_num_steps[adaptive_level]);
    alpha_i[adaptive_level].Allocate(image_num_frequencies, num_pix,
        geodesic_num_steps[adaptive_level]);
  }
  j_i[adaptive_level].Zero();
  alpha_i[adaptive_level].Zero();

  // Go through rays in parallel
  #pragma omp parallel for schedule(static)
  for (int m = 0; m < num_pix; m++)
  {
    // Check number of steps
    int num_steps = sample_num[adaptive_level](m);
    if (num_steps <= 0)
      continue;

    // Set pixel to NaN if ray has problem
    if (fallback_nan and sample_flags[adaptive_level](m))
    {
      for (int n = 0; n < num_steps; n++)
      {
        j_i[adaptive_level](m,n) = std::numeric_limits<double>::quiet_NaN();
        alpha_i[adaptive_level](m,n) = std::numeric_limits<double>::quiet_NaN();
      }
      continue;
    }

    // Go through samples
    for (int n = 0; n < num_steps; n++)
    {
      // Extract geodesic position and momentum
      double x = sample_pos[adaptive_level](m,n,1);
      double y = sample_pos[adaptive_level](m,n,2);
      double z = sample_pos[adaptive_level](m,n,3);
      double k_0 = sample_dir[adaptive_level](m,n,0);
      double k_1 = sample_dir[adaptive_level](m,n,1);
      double k_2 = sample_dir[adaptive_level](m,n,2);
      double k_3 = sample_dir[adaptive_level](m,n,3);

      // Cut outside camera radius
      double r = RadialGeodesicCoordinate(x, y, z);
      if (r > camera_r)
        continue;

      // Cut camera plane
      if (cut_omit_near or cut_omit_far)
      {
        double dot_product = x * camera_x[1] + y * camera_x[2] + z * camera_x[3];
        if ((cut_omit_near and dot_product > 0.0) or (cut_omit_far and dot_product < 0.0))
          continue;
      }

      // Cut spheres
      if ((cut_omit_in >= 0.0 and r < cut_omit_in) or (cut_omit_out >= 0.0 and r > cut_omit_out))
        continue;

      // Cut with respect to midplane
      if (cut_midplane_theta > 0.0 or cut_midplane_theta < 0.0)
      {
        double th = std::acos(z / r);
        if ((cut_midplane_theta > 0.0 and std::abs(th - Math::pi / 2.0) > cut_midplane_theta)
            or (cut_midplane_theta < 0.0 and std::abs(th - Math::pi / 2.0) < -cut_midplane_theta))
        {
          sample_cut[adaptive_level](m,n) = true;
          continue;
        }
      }
      if ((cut_midplane_z > 0.0 and std::abs(z) > cut_midplane_z)
          or (cut_midplane_z < 0.0 and std::abs(z) < -cut_midplane_z))
      {
        sample_cut[adaptive_level](m,n) = true;
        continue;
      }

      // Cut arbitrary plane
      if (cut_plane)
      {
        double dot_product = (x - cut_plane_origin_x) * cut_plane_normal_x
            + (y - cut_plane_origin_y) * cut_plane_normal_y
            + (z - cut_plane_origin_z) * cut_plane_normal_z;
        if (dot_product < 0.0)
          continue;
      }


      // Calculate curvilinear coordinates
      double rr = std::sqrt(r * r - z * z);
      double cth = z / r;
      double sth = std::sqrt(1.0 - cth * cth);
      double ph = std::atan2(y, x) - std::atan(bh_a / r);
      double sph = std::sin(ph);
      double cph = std::cos(ph);

      // Calculate metric
      double delta = r * r - 2.0 * bh_m * r + bh_a * bh_a;
      double sigma = r * r + bh_a * bh_a * cth * cth;
      double gtt_bl = -(1.0 + 2.0 * bh_m * r * (r * r + bh_a * bh_a) / (delta * sigma));
      double gtph_bl = -2.0 * bh_m * bh_a * r / (delta * sigma);
      double grr_bl = delta / sigma;
      double gthth_bl = 1.0 / sigma;
      double gphph_bl = (sigma - 2.0 * bh_m * r) / (delta * sigma * sth * sth);

      // Calculate angular momentum (C 6)
      double ll = formula_l0 / (1.0 + rr) * std::pow(rr, 1.0 + formula_q);

      // Calculate 4-velocity (C 7-8)
      double u_norm = 1.0 / std::sqrt(-gtt_bl + 2.0 * gtph_bl * ll - gphph_bl * ll * ll);
      double u_t_bl = -u_norm;
      double u_r_bl = 0.0;
      double u_th_bl = 0.0;
      double u_ph_bl = u_norm * ll;
      double ut_bl = gtt_bl * u_t_bl + gtph_bl * u_ph_bl;
      double ur_bl = grr_bl * u_r_bl;
      double uth_bl = gthth_bl * u_th_bl;
      double uph_bl = gtph_bl * u_t_bl + gphph_bl * u_ph_bl;
      double ut = ut_bl + 2.0 * bh_m * r / delta * ur_bl;
      double ur = ur_bl;
      double uth = uth_bl;
      double uph = uph_bl + bh_a / delta * ur_bl;
      double u0 = ut;
      double u1 =
          sth * cph * ur + cth * (r * cph - bh_a * sph) * uth + sth * (-r * sph - bh_a * cph) * uph;
      double u2 =
          sth * sph * ur + cth * (r * sph + bh_a * cph) * uth + sth * (r * cph - bh_a * sph) * uph;
      double u3 = cth * ur - r * sth * uth;

      // Calculate fluid-frame number density (C 5)
      double n_n0_fluid =
          std::exp(-0.5 * (r * r / (formula_r0 * formula_r0) + formula_h * formula_h * cth * cth));
      

      // Go through frequencies
      for (int l = 0; l < image_num_frequencies; l++)
      {
        if(formula_name=="Gold+2020"){
        // Calculate frequency in CGS units
        double nu_fluid_cgs = -(u0 * k_0 + u1 * k_1 + u2 * k_2 + u3 * k_3) * image_frequencies(l)
            * momentum_factors[adaptive_level](m);

        // Calculate emission coefficient in CGS units (C 9-10)
        double j_nu_fluid_cgs =
            formula_cn0 * n_n0_fluid * std::pow(nu_fluid_cgs / formula_nup, -formula_alpha);
        j_i[adaptive_level](l,m,n) = j_nu_fluid_cgs / (nu_fluid_cgs * nu_fluid_cgs);
        
        // Calculate absorption coefficient in CGS units (C 11-12)
        double alpha_nu_fluid_cgs = formula_a * formula_cn0 * n_n0_fluid
            * std::pow(nu_fluid_cgs / formula_nup, -formula_beta - formula_alpha);
        alpha_i[adaptive_level](l,m,n) = alpha_nu_fluid_cgs * nu_fluid_cgs;
        }

        //for spherical case. we just want a spherical ball of unmoving gas with constant rho and T
        //the below is modified from the simulation coefficients code to achieve this. you simply set uu1,2,3 to 0 and use minkowski values for both simulation and geodesic metrics
        if(formula_name == "spherical" || formula_name == "disk"){
          if(formula_name == "spherical"){
            if(r > formula_r_out){
              //outside of the sphere, no emission or absorption
              j_i[adaptive_level](l,m,n) = 0.0;
              alpha_i[adaptive_level](l,m,n) = 0.0;
              continue;
            }
          }else{
            if(r > formula_r_out){
              //outside of the disk, no emission or absorption
              j_i[adaptive_level](l,m,n) = 0.0;
              alpha_i[adaptive_level](l,m,n) = 0.0;
              continue;
            }
            if(z > formula_height/2 || z < -formula_height/2){
              //outside of the disk, no emission or absorption
              j_i[adaptive_level](l,m,n) = 0.0;
              alpha_i[adaptive_level](l,m,n) = 0.0;
              continue;
            }
          }
        
        //std::printf("calculating spherical/disk formula coefficients at r=%.5f, theta=%.5f\n",r,std::acos(z/r));
        double kcov[4];
        kcov[0] = sample_dir[adaptive_level](m,n,0);
        kcov[1] = sample_dir[adaptive_level](m,n,1);
        kcov[2] = sample_dir[adaptive_level](m,n,2);
        kcov[3] = sample_dir[adaptive_level](m,n,3);

        double gcon[4][4];
        double gcov[4][4];
        //set to minkowski
        for(int mu=0; mu<4; mu++)
          for(int nu=0; nu<4; nu++)
          {
            gcov[mu][nu]=0.0;
            gcon[mu][nu]=0.0;
          }
        gcov[0][0]=-1.0;
        gcov[1][1]=1.0;
        gcov[2][2]=1.0;
        gcov[3][3]=1.0;
        gcon[0][0]=-1.0;
        gcon[1][1]=1.0;
        gcon[2][2]=1.0;
        gcon[3][3]=1.0;

        //set the velocity to 0 
        double ucon[4];
        ucon[0]=1.0;
        ucon[1]=0.0;
        ucon[2]=0.0;
        ucon[3]=0.0;
        double ucov[4];
        ucov[0]=-1.0;
        ucov[1]=0.0;
        ucov[2]=0.0;
        ucov[3]=0.0;

        //set the magnetic field to 0
        double bcon[4];
        bcon[0]=0.0;
        bcon[1]=0.0;
        bcon[2]=0.0;
        bcon[3]=0.0;
        

        // Calculate geodesic contravariant momentum
        double kcon[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            kcon[mu] += gcon[mu][nu] * kcov[nu];


            
          double nu_cgs = 0.0;
          for (int mu = 0; mu < 4; mu++)
            nu_cgs -= kcov[mu] * ucon[mu];//this gives the fluid frame frequency
          nu_cgs *= image_frequencies(l) * momentum_factors[adaptive_level](m);
          double nu_2_cgs = nu_cgs * nu_cgs;

          double kb_tt_e_cgs = formula_T*Physics::k_b; //formula temperature in K
          double rho_cgs = formula_rho; //formula density in g/cm^3
          double n_cgs = rho_cgs / (plasma_mu * Physics::m_p);

          //plasma_ne_ni is set through our input parameters as 1 so basically number density for both is equal everywhere
          double n_e_cgs = n_cgs*plasma_ne_ni;
          double n_i_cgs = n_cgs;


          //TEGAN: here is the table opacity 
          double table_opacity_value=0.0;
          
          bool default_to_free_free = false;
          if(opacity_table){
            //if we have NaN frequency, just treat it as blacklight does normally
            if(nu_cgs != nu_cgs){
              default_to_free_free = true;
            }else{
            double log_freq = std::log10(nu_cgs*Physics::h);
            double log_temp = std::log10(kb_tt_e_cgs/Physics::k_b);
            double log_rho = std::log10(rho_cgs);
            if(log_freq < p_opacity_table_reader->fmin || log_freq > p_opacity_table_reader->fmax){
              default_to_free_free = true;
            }
            if(log_temp < p_opacity_table_reader->tmin || log_temp > p_opacity_table_reader->tmax){
              default_to_free_free = true;
            }
            if(log_rho < p_opacity_table_reader->rmin || log_rho > p_opacity_table_reader->rmax){
              default_to_free_free = true;
            }
          
            if(!default_to_free_free){
              double xi = (log_rho - p_opacity_table_reader->rmin)/p_opacity_table_reader->dlr;
              double xj = (log_temp - p_opacity_table_reader->tmin)/p_opacity_table_reader->dlt;
              double xk = (log_freq - p_opacity_table_reader->fmin)/p_opacity_table_reader->dlf;
              
              int i = std::floor(xi);
              int j = std::floor(xj);
              int k = std::floor(xk);

              double xd = xi - i;
              double yd = xj - j;
              double zd = xk - k;
              double c00 = p_opacity_table_reader->plan_tab(k,j,i)*(1 - xd) + p_opacity_table_reader->plan_tab(k+1,j,i)*xd;
              double c01 = p_opacity_table_reader->plan_tab(k,j,i+1)*(1 - xd) + p_opacity_table_reader->plan_tab(k+1,j,i+1)*xd;
              double c10 = p_opacity_table_reader->plan_tab(k,j+1,i)*(1 - xd) + p_opacity_table_reader->plan_tab(k+1,j+1,i)*xd;
              double c11 = p_opacity_table_reader->plan_tab(k,j+1,i+1)*(1 - xd) + p_opacity_table_reader->plan_tab(k+1,j+1,i+1)*xd;
              double c0 = c00*(1 - yd) + c10*yd;
              double c1 = c01*(1 - yd) + c11*yd;
              table_opacity_value = c0*(1 - zd) + c1*zd;

              alpha_i[adaptive_level](l,m,n)+= table_opacity_value*nu_cgs;
              
              double partA = 4*pow(Physics::e,6.)/(3*Physics::m_e*Physics::c*Physics::h);
              double partB = std::sqrt(2.0*Math::pi/(3.0*kb_tt_e_cgs*Physics::m_e));
              double gaunt_factor = 1.0; //approximate it as this because shouldn't impact too much
              
              double coefficient = partA*partB*n_e_cgs*n_i_cgs*(1.0 - std::exp(-Physics::h*nu_cgs/kb_tt_e_cgs))*gaunt_factor/(nu_cgs*nu_cgs*nu_cgs);
              

              double planck_function = 2.0 * Physics::h * nu_cgs * nu_cgs * nu_cgs
                  / (Physics::c * Physics::c) / std::expm1(Physics::h * nu_cgs / kb_tt_e_cgs);
              j_i[adaptive_level](l,m,n) += table_opacity_value* planck_function/(nu_cgs*nu_cgs);
            }
          }
          }

          //Calculate thermal free-free emissivities (Rybicki & Lightman, eqn 5.14a)
          if (image_free_free and (!opacity_table || default_to_free_free))
          {
           double partA = 16.0*std::pow(Physics::e,6.)/(3.0*Physics::m_e*std::pow(Physics::c,3.));
           double partB = std::sqrt(2.0*Math::pi/(3.0*kb_tt_e_cgs*Physics::m_e));
           double gaunt_factor = 1.0; //approximate it as this because shouldn't impact too much

           double coefficient = partA*partB*n_e_cgs*n_i_cgs*std::exp(-Physics::h*nu_cgs/kb_tt_e_cgs)*gaunt_factor;
           
            if (image_light or image_emission or image_emission_ave)
              j_i[adaptive_level](l,m,n) += coefficient/(nu_cgs*nu_cgs);

            if (image_light and image_polarization)
            {
              //assume that there's no polarizational bremsstrahlung
              j_q[adaptive_level](l,m,n) += 0.0;
              j_v[adaptive_level](l,m,n) += 0.0;
            }
          }

          //Calculate thermal free-free absorptivities (Rybicki & Lightman, eqn 5.18a)
          if ((image_light or image_tau or image_tau_int) and image_free_free and (!opacity_table || default_to_free_free))
          {
           double partA = 4*pow(Physics::e,6.)/(3*Physics::m_e*Physics::c*Physics::h);
           double partB = std::sqrt(2.0*Math::pi/(3.0*kb_tt_e_cgs*Physics::m_e));
           double gaunt_factor = 1.0; //approximate it as this because shouldn't impact too much
           
           double coefficient = partA*partB*n_e_cgs*n_i_cgs*(1.0 - std::exp(-Physics::h*nu_cgs/kb_tt_e_cgs))*gaunt_factor/(nu_cgs*nu_cgs*nu_cgs);

            if (image_light or image_emission or image_emission_ave)
            alpha_i[adaptive_level](l,m,n) += coefficient*nu_cgs;
            if (image_light and image_polarization)
            {
              alpha_q[adaptive_level](l,m,n) += 0.0;
              alpha_v[adaptive_level](l,m,n) += 0.0;
            }

            // Account for numerical issues later arising from absorptivities being too small
            //(taken from the synchrotron treatment above mostly)
            if ((image_light or image_tau or image_tau_int)
                and 1.0 / (alpha_i[adaptive_level](l,m,n) * alpha_i[adaptive_level](l,m,n))
                == std::numeric_limits<double>::infinity())
            {
              alpha_i[adaptive_level](l,m,n) += 0.0;
              if (image_light and image_polarization)
              {
                alpha_q[adaptive_level](l,m,n) += 0.0;
                alpha_v[adaptive_level](l,m,n) += 0.0;
              }
            }
          }
          
        }
      }
    }
  }
  return;
}
