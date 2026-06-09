// Blacklight opacity table reader

// C++ headers
#include <algorithm>  // remove
#include <cctype>     // tolower
#include <cmath>      // abs, pow
#include <cstdint>    // int32_t
#include <cstdio>     // snprintf
#include <cstring>    // strncmp, strtok
#include <fstream>    // ifstream
#include <ios>        // ios_base, streamoff
#include <iosfwd>     // streampos
#include <optional>   // optional
#include <sstream>    // ostringstream
#include <string>     // getline, stod, stoi, string

// Library headers
#include <omp.h>  // pragmas, omp_get_wtime

// Blacklight headers
#include "mc_reader.hpp"
#include "../blacklight.hpp"                 // Math, enums
#include "../input_reader/input_reader.hpp"  // InputReader
#include "../simulation_reader/simulation_reader.hpp"
#include "../utils/array.hpp"                // Array
#include "../utils/exceptions.hpp"           // BlacklightException, BlacklightWarning
#include "../utils/file_io.hpp"              // ReadBinary

//--------------------------------------------------------------------------------------------------

// Simulation reader constructor
// Inputs:
//   p_input_reader_: pointer to object containing input parameters
// Notes:
//   File is not opened for writing until Read() function is called because the file name might be
//       reformatted after this constructor is called.
MCReader::MCReader(const InputReader *p_input_reader_, const SimulationReader *p_simulation_reader_)
  : p_input_reader(p_input_reader_), p_simulation_reader(p_simulation_reader_)
{
  // Copy general input data
  model_type = p_input_reader->model_type.value();
  mc_input = p_input_reader->mc_input.value();
  
  if(mc_input){
    mc_file_name = p_input_reader->mc_file.value();
    mc_freq_file_name = p_input_reader->mc_freq_file.value();
    simulation_coord = p_input_reader->simulation_coord.value();
    simulation_m_msun = p_input_reader->simulation_m_msun.value();
    simulation_all_cgs = p_input_reader->simulation_all_cgs.value();
    compton = p_input_reader->compton.value();
    stimulated_compton = p_input_reader->stimulated_compton.value();
  
    if(simulation_all_cgs){
      simulation_rho_cgs = 1.0;
      simulation_r_rg = Physics::c*Physics::c/(simulation_m_msun*Physics::gg_msun);
      simulation_v_c = 1/Physics::c;
    }else{
      simulation_rho_cgs = p_input_reader->simulation_rho_cgs.value();
      simulation_r_rg = p_input_reader->simulation_r_rg.value();
      simulation_v_c = p_input_reader->simulation_v_c.value();
    }
    simulation_hd_only = p_input_reader->simulation_hd_only.value();
    simulation_mc_temp = p_input_reader->simulation_mc_temp.value();
    plasma_model = p_input_reader->plasma_model.value();
    if (plasma_model == PlasmaModel::ti_te_beta)
      {
        plasma_use_p = p_input_reader->plasma_use_p.value();
        plasma_rat_low = p_input_reader->plasma_rat_low.value();
        plasma_rat_high = p_input_reader->plasma_rat_high.value();
      }

    plasma_mu = p_input_reader->plasma_mu.value();
    plasma_ne_ni = p_input_reader->plasma_ne_ni.value();
    plasma_power_frac = p_input_reader->plasma_power_frac.value();
    plasma_kappa_frac = p_input_reader->plasma_kappa_frac.value();

    plasma_thermal_frac = 1.0 - (plasma_power_frac + plasma_kappa_frac);
    if (plasma_thermal_frac < 0.0 or plasma_thermal_frac > 1.0)
      BlacklightWarning("Fraction of thermal electrons outside [0, 1].");

    simulation_format = p_input_reader->simulation_format.value();
    scattering_source_terms = new Array<float>[1];
    grid_prim = p_simulation_reader->prim;
    ind_rho = p_simulation_reader->ind_rho;
    ind_pgas = p_simulation_reader->ind_pgas;
    ind_kappa = p_simulation_reader->ind_kappa;
    ind_uu1 = p_simulation_reader->ind_uu1;
    ind_uu2 = p_simulation_reader->ind_uu2;
    ind_uu3 = p_simulation_reader->ind_uu3;
    if( !simulation_hd_only){
      ind_bb1 = p_simulation_reader->ind_bb1;
      ind_bb2 = p_simulation_reader->ind_bb2;
      ind_bb3 = p_simulation_reader->ind_bb3;
    }
    bh_m = 1.0;
    bh_a = p_input_reader->simulation_a.value();
  }
}

//--------------------------------------------------------------------------------------------------

// Simulation reader destructor
MCReader::~MCReader()
{
  if (num_freqs > 0){
    freq_grid.Deallocate();
    ln_freq_grid.Deallocate();
  }

  if (num_dataset_names > 0)
    delete[] dataset_names;
  if (num_variable_names > 0)
    delete[] variable_names;
  scattering_source_terms[0].Deallocate();
  delete[] scattering_source_terms;
}

//--------------------------------------------------------------------------------------------------

// Opacity table reader read and initialize function
// Inputs:
//   snapshot: index (starting at 0) of which snapshot is about to be prepared
// Outputs:
//   returned value: execution time in seconds
// Notes:
//   Reads opacity table from file specified in input file.
//   This function is largely adapted from code written by Shane Davis
//   specifically the function MonteCarlo::InitUserMonteCarloData() in athena/src/pgen/mc_tde.cpp
double MCReader::Read(int snapshot)
{
  // Only proceed if needed
  if (model_type != ModelType::simulation)
    return 0.0;
  if(!mc_input)
    return 0.0;

  std::cout << "Reading MC results from file..." << mc_file_name << std::endl;
  double time_start = omp_get_wtime();

  int num_read = 1;
  ReadFreqFile();
  dlf = std::log10(freq_grid(1))-std::log10(freq_grid(0));


  // Open input file

  data_stream =
        std::ifstream(mc_file_name, std::ios_base::in | std::ios_base::binary);
  if (not data_stream.is_open())
      throw BlacklightException("Could not open MC file for reading.");

  if (simulation_format == SimulationFormat::athenak or simulation_format == SimulationFormat::iharm3d)
    {
      throw BlacklightException("Reading MC files is only implemented for Athena++ formats.");
    }
  

    // Read coordinates
    if (first_time)
    {
        //the following way of getting the hdf5 file fails because it doesn't match the expected signature
        data_stream = std::ifstream(mc_file_name,std::ios_base::in|std::ios_base::binary);
        if (not data_stream.is_open())
              throw BlacklightException("Could not open file for reading.");
        
        ReadHDF5Superblock();
        root_data_segment_address = ReadHDF5Heap(root_name_heap_address);
        ReadHDF5RootObjectHeader();


        ReadHDF5IntArray("Levels", levels);

        ReadHDF5FloatArray("x1f", x1f);
        ReadHDF5FloatArray("x2f", x2f);
        ReadHDF5FloatArray("x3f", x3f);
        ReadHDF5FloatArray("x1v", x1v);
        ReadHDF5FloatArray("x2v", x2v);
        ReadHDF5FloatArray("x3v", x3v);

        //Update radial boundaries and check all of the boundaries match
        for(int j=0;j<x1f.n1;j++){
          for(int i=0;i<x1f.n2;i++){
            x1f(i,j) = simulation_r_rg*x1f(i,j);
            x1v(i,j) = simulation_r_rg*x1v(i,j);
            if(std::fabs(x1f(i,j)-p_simulation_reader->x1f(i,j))>1e-8 ){
              std::printf("i: %d j: %d , MC x1f and x1v values : %.15f %.15f Regular x1f and x1v values : %.15f %.15f \n",i,j,x1f(i,j),x1v(i,j),p_simulation_reader->x1f(i,j),p_simulation_reader->x1v(i,j));
              throw BlacklightException("MC X1 grid bounds don't match file bounds within 1e-8.");
            }
          }
        }
        for(int j=0;j<x2f.n1;j++){
          for(int i=0;i<x2f.n2;i++){
            if(simulation_coord==Coordinates::cks){
              x2f(i,j) = simulation_r_rg*x2f(i,j);
              x2v(i,j) = simulation_r_rg*x2v(i,j);
            }
            if(std::fabs(x2f(i,j)-p_simulation_reader->x2f(i,j))>1e-8 ){
              std::printf("i: %d j: %d , MC x2f and x2v values : %.15f %.15f Regular x2f and x2v values : %.15f %.15f \n",i,j,x2f(i,j),x2v(i,j),p_simulation_reader->x2f(i,j),p_simulation_reader->x2v(i,j));
              throw BlacklightException("MC X2 grid bounds don't match file bounds within 1e-8");
            }
          }
        }
        for(int j=0;j<x3f.n1;j++){
          for(int i=0;i<x3f.n2;i++){
            if(simulation_coord==Coordinates::cks){
              x3f(i,j) = simulation_r_rg*x3f(i,j);
              x3v(i,j) = simulation_r_rg*x3v(i,j);
            }
            if(std::fabs(x3f(i,j)-p_simulation_reader->x3f(i,j))>1e-8 ){
              throw BlacklightException("MC X3 grid bounds don't match file bounds within 1e-8");
            }
          }
        }
    } 

    Array<float> *scattering;
    scattering = new Array<float>[1];
    Array<float> *scattering_first_derivs;
    scattering_first_derivs = new Array<float>[1];
    Array<float> *scattering_second_derivs;
    scattering_second_derivs = new Array<float>[1];
    
    if (first_time)
      {
        int n5 = num_freqs;
        int n3 = x3v.n1; 
        int n2 = x2v.n1; 
        int n1 = x1v.n1;
        int n4 = levels.n1;
        scattering[0].Allocate(n5,n4, n3, n2, n1);
        if(compton){
          scattering_first_derivs[0].Allocate(n5,n4, n3, n2, n1);
          scattering_second_derivs[0].Allocate(n5,n4, n3, n2, n1);
        }else{
          scattering_first_derivs[0].Allocate(n1);
          scattering_second_derivs[0].Allocate(n1);
        }
        scattering_source_terms[0].Allocate(n5,n4, n3, n2, n1);
      }
      //Array<float> source_terms(scattering_source_terms[0]);
      Array<float> shallow_scatter(scattering[0]);
      Array<float> shallow_first_deriv(scattering_first_derivs[0]);
      Array<float> shallow_second_deriv(scattering_second_derivs[0]);
      Array<float> shallow_source_terms(scattering_source_terms[0]);
      
      ReadHDF5FloatArray("mcscat",shallow_scatter);
      if(compton){
        Gradient(shallow_first_deriv, shallow_scatter, ln_freq_grid);
        Gradient(shallow_second_deriv, shallow_first_deriv, ln_freq_grid);
      }
      CalculateSourceTerm(shallow_source_terms, shallow_scatter, shallow_first_deriv, shallow_second_deriv);
      

    // Close input file
    data_stream.close();

    scattering[0].Deallocate();
    delete[] scattering;
    scattering_first_derivs[0].Deallocate();
    delete[] scattering_first_derivs;
    scattering_second_derivs[0].Deallocate();
    delete[] scattering_second_derivs;

    // Update first time flag
    first_time = false;

  // Calculate elapsed time
  return omp_get_wtime() - time_start;

}



// Opacity table reader read and initialize function
// Inputs:
//   snapshot: index (starting at 0) of which snapshot is about to be prepared
// Outputs:
//   returned value: execution time in seconds
// Notes:
//   Reads opacity table from file specified in input file.
//   This function is largely adapted from code written by Shane Davis
//   specifically the function MonteCarlo::InitUserMonteCarloData() in athena/src/pgen/mc_tde.cpp
void MCReader::ReadFreqFile()
{
  // Open input file
  if ( (mc_freq_file=fopen(mc_freq_file_name.c_str(),"r"))==NULL) {
    std::stringstream msg;
    msg << "FATAL ERROR: Could not open " << mc_freq_file_name << "." << std::endl;
    throw BlacklightException(msg.str().c_str());
  }


  
  fscanf(mc_freq_file,"%d",&(num_freqs));


  // allocate arrays
  freq_grid.Allocate(num_freqs);
  ln_freq_grid.Allocate(num_freqs);

  double trash;
  for(int i=0; i<num_freqs; ++i){
    fscanf(mc_freq_file,"%d",&(trash));
    fscanf(mc_freq_file,"%lf",&(freq_grid(i)));
    //std::cout<<freq_grid(i)<<std::endl;
    fscanf(mc_freq_file,"%lf",&(trash));
    fscanf(mc_freq_file,"%lf",&(trash));
    ln_freq_grid(i) = std::log(freq_grid(i));
  }

  fclose(mc_freq_file);

  return;

}


//--------------------------------------------------------------------------------------------------

// Function for evaluating the gradient of f over x 
// Used specifically to calculate the gradient of J over frequency for the Compton source term
void MCReader::Gradient(Array<float> &grad,Array<float> &f, Array<double> &x){
  int nx = x.n1;
  if(nx!=num_freqs){
    std::printf("nx: %d mc_num_freqs: %d",nx,num_freqs);
  }
  #pragma omp parallel for schedule(static) collapse(4)
  for(int i=0;i<f.n1;i++){
    for(int j=0;j<f.n2;j++){
      for(int k=0;k<f.n3;k++){
        for(int b=0;b<f.n4;b++){
          grad(0,b,k,j,i) = (f(1,b,k,j,i) - f(0,b,k,j,i))/(x(1) - x(0));
          grad(nx-1,b,k,j,i) = (f(nx-1,b,k,j,i) - f(nx-2,b,k,j,i))/(x(nx-1) - x(nx-2));
          for(int l=1;l<(f.n5-1);l++){
            grad(l,b,k,j,i) = (f(l+1,b,k,j,i) - f(l-1,b,k,j,i))/(x(l+1) - x(l-1));
          }
        }
      }
    }
  }
} 

void MCReader::CalculateSourceTerm(Array<float> &source_term,Array<float> &scattering,Array<float> &scattering_prime,Array<float> &scattering_prime_prime){
  // Calculate the source term using the scattering and its derivatives

  double e_unit = simulation_rho_cgs * Physics::c * Physics::c * simulation_v_c*simulation_v_c;
  double b_unit = std::sqrt(4.0 * Math::pi * e_unit);

  #pragma omp parallel for schedule(static) collapse(4)
  for(int i=0;i<source_term.n1;i++){
    for(int j=0;j<source_term.n2;j++){
      for(int k=0;k<source_term.n3;k++){
        for(int b=0;b<source_term.n4;b++){
          double rho_cgs = simulation_rho_cgs*grid_prim[0](ind_rho,b,k,j,i);
          double pgas_cgs = e_unit*grid_prim[0](p_simulation_reader->ind_pgas,b,k,j,i);
          double gcov_sim[4][4];
          double gcon_sim[4][4];
          

          // Calculate electron temperature for model with T_i/T_e a function of beta (E1 1)
          double kb_tt_e_cgs = std::numeric_limits<double>::quiet_NaN();
          double theta_e = std::numeric_limits<double>::quiet_NaN();
          if (plasma_thermal_frac != 0.0 and plasma_model == PlasmaModel::ti_te_beta)
          {
            double uu1_sim = grid_prim[0](p_simulation_reader->ind_uu1,b,k,j,i);
            double uu2_sim = grid_prim[0](p_simulation_reader->ind_uu2,b,k,j,i);
            double uu3_sim = grid_prim[0](p_simulation_reader->ind_uu3,b,k,j,i);
            double bb1_sim = grid_prim[0](p_simulation_reader->ind_bb1,b,k,j,i);
            double bb2_sim = grid_prim[0](p_simulation_reader->ind_bb2,b,k,j,i);
            double bb3_sim = grid_prim[0](p_simulation_reader->ind_bb3,b,k,j,i);
            uu1_sim *=simulation_v_c;
            uu2_sim *=simulation_v_c;
            uu3_sim *=simulation_v_c;

            // Calculate simulation metric
            if(simulation_coord == Coordinates::cks){
              CovariantSimulationMetric(simulation_r_rg*x1v(i,j), simulation_r_rg*x2v(i,j), x3v(i,j), gcov_sim);
              ContravariantSimulationMetric(simulation_r_rg*x1v(i,j), simulation_r_rg*x2v(i,j), x3v(i,j), gcon_sim);
            }else{
              CovariantSimulationMetric(simulation_r_rg*x1v(i,j), x2v(i,j), x3v(i,j), gcov_sim);
              ContravariantSimulationMetric(simulation_r_rg*x1v(i,j), x2v(i,j), x3v(i,j), gcon_sim);
            }

            // Calculate simulation velocity
            double uu0_sim = std::sqrt(1.0 + gcov_sim[1][1] * uu1_sim * uu1_sim
                + 2.0 * gcov_sim[1][2] * uu1_sim * uu2_sim + 2.0 * gcov_sim[1][3] * uu1_sim * uu3_sim
                + gcov_sim[2][2] * uu2_sim * uu2_sim + 2.0 * gcov_sim[2][3] * uu2_sim * uu3_sim
                + gcov_sim[3][3] * uu3_sim * uu3_sim);
            double lapse_sim = 1.0 / std::sqrt(-gcon_sim[0][0]);
            double shift1_sim = -gcon_sim[0][1] / gcon_sim[0][0];
            double shift2_sim = -gcon_sim[0][2] / gcon_sim[0][0];
            double shift3_sim = -gcon_sim[0][3] / gcon_sim[0][0];
            double ucon_sim[4];
            ucon_sim[0] = uu0_sim / lapse_sim;
            ucon_sim[1] = uu1_sim - shift1_sim * uu0_sim / lapse_sim;
            ucon_sim[2] = uu2_sim - shift2_sim * uu0_sim / lapse_sim;
            ucon_sim[3] = uu3_sim - shift3_sim * uu0_sim / lapse_sim;
            double ucov_sim[4] = {};
            for (int mu = 0; mu < 4; mu++)
              for (int nu = 0; nu < 4; nu++)
                ucov_sim[mu] += gcov_sim[mu][nu] * ucon_sim[nu];

            // Calculate simulation magnetic field
            double bcon_sim[4];
            bcon_sim[0] = ucov_sim[1] * bb1_sim + ucov_sim[2] * bb2_sim + ucov_sim[3] * bb3_sim;
            bcon_sim[1] = (bb1_sim + bcon_sim[0] * ucon_sim[1]) / ucon_sim[0];
            bcon_sim[2] = (bb2_sim + bcon_sim[0] * ucon_sim[2]) / ucon_sim[0];
            bcon_sim[3] = (bb3_sim + bcon_sim[0] * ucon_sim[3]) / ucon_sim[0];
            double bcov_sim[4] = {};
            for (int mu = 0; mu < 4; mu++)
              for (int nu = 0; nu < 4; nu++)
                bcov_sim[mu] += gcov_sim[mu][nu] * bcon_sim[nu];
            double b_sq = 0.0;
            for (int mu = 0; mu < 4; mu++)
              b_sq += bcov_sim[mu] * bcon_sim[mu];
            double bb_cgs = std::sqrt(b_sq) * b_unit;
            double sigma = b_sq / grid_prim[0](ind_rho,b,k,j,i);
            double beta_inv = b_sq / (2.0 * grid_prim[0](p_simulation_reader->ind_pgas,b,k,j,i));

            double tti_tte = (plasma_rat_high + plasma_rat_low * beta_inv * beta_inv)
                / (1.0 + beta_inv * beta_inv);
            double kb_tt_tot_cgs = plasma_mu * Physics::m_p * pgas_cgs / rho_cgs;
            
            if (plasma_use_p)
              kb_tt_e_cgs = (1.0 + plasma_ne_ni) / (tti_tte + plasma_ne_ni) * kb_tt_tot_cgs;
            else
            {
              kb_tt_e_cgs = (1.0 + plasma_ne_ni) * kb_tt_tot_cgs / (plasma_gamma - 1.0);
              kb_tt_e_cgs /= tti_tte / (plasma_gamma_i - 1.0) + plasma_ne_ni / (plasma_gamma_e - 1.0);
            }
            
            theta_e = kb_tt_e_cgs / (Physics::m_e * Physics::c * Physics::c);
          }
          if(plasma_thermal_frac!=0.0 and plasma_model == PlasmaModel::one_temp)
          {
            double kb_tt_tot_cgs = 0.0;
            if(simulation_mc_temp){
              kb_tt_tot_cgs =  Physics::m_p *pgas_cgs / rho_cgs;
            }else{
              kb_tt_tot_cgs = plasma_mu * Physics::m_p *pgas_cgs / rho_cgs;
            }
            
            kb_tt_e_cgs = kb_tt_tot_cgs;
            
            theta_e = kb_tt_e_cgs / (Physics::m_e * Physics::c * Physics::c);

          }
          // Calculate electron temperature for given electron entropy (E2 13)
          if (plasma_thermal_frac != 0.0 and plasma_model == PlasmaModel::code_kappa)
          {
            double mu_e = plasma_mu * (1.0 + 1.0 / plasma_ne_ni);
            double rho_e = grid_prim[0](ind_rho,b,k,j,i) * Physics::m_e / (mu_e * Physics::m_p);
            double rho_kappa_e_cbrt = std::cbrt(rho_e * grid_prim[0](ind_kappa,b,k,j,i));
            theta_e = 1.0 / 5.0 * (std::sqrt(1.0 + 25.0 * rho_kappa_e_cbrt * rho_kappa_e_cbrt) - 1.0);
            kb_tt_e_cgs = theta_e * Physics::m_e * Physics::c * Physics::c;
          }

          for(int l=0;l<source_term.n5;l++){
            //calculate frequency dependent compton source term
            double x = Physics::h*freq_grid(l)/(Physics::m_e*Physics::c*Physics::c);
            if(compton){
              source_term(l,b,k,j,i) = (1-x)*scattering(l,b,k,j,i)+ (x-3*theta_e)*scattering_prime(l,b,k,j,i)+theta_e*scattering_prime_prime(l,b,k,j,i);
              if(stimulated_compton){
                source_term(l,b,k,j,i) += Physics::c*Physics::c*x*scattering(l,b,k,j,i)*(scattering_prime(l,b,k,j,i)-scattering(l,b,k,j,i))/(Physics::h*freq_grid(l)*freq_grid(l)*freq_grid(l));
              }
            }else{
              source_term(l,b,k,j,i) = scattering(l,b,k,j,i);
            }
          }
        }
      }
    }
  }
}

