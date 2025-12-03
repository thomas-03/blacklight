// Blacklight simulation reader

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
#include "opacity_table_reader.hpp"
#include "../blacklight.hpp"                 // Math, enums
#include "../input_reader/input_reader.hpp"  // InputReader
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
OpacityTableReader::OpacityTableReader(const InputReader *p_input_reader_)
  : p_input_reader(p_input_reader_)
{
  // Copy general input data
  model_type = p_input_reader->model_type.value();
  opacity_table = p_input_reader->opacity_table.value();
  opacity_file = p_input_reader->opacity_file.value();

  // Copy simulation parameters
  /*if (model_type == ModelType::simulation)
  {
    simulation_file = p_input_reader->simulation_file.value();
    simulation_multiple = p_input_reader->simulation_multiple.value();
    if (simulation_multiple)
    {
      simulation_start = p_input_reader->simulation_start.value();
      if (simulation_start < 0)
        throw BlacklightException("Must have nonnegative index simulation_start.");
      simulation_end = p_input_reader->simulation_end.value();
      if (simulation_end < simulation_start)
        throw
            BlacklightException("Must have simulation_end at least as large as simulation_start.");
    }
    simulation_coord = p_input_reader->simulation_coord.value();
    simulation_a = p_input_reader->simulation_a.value();
    simulation_m_msun = p_input_reader->simulation_m_msun.value();
    simulation_rho_cgs = p_input_reader->simulation_rho_cgs.value();
    simulation_v_cgs = p_input_reader->simulation_v_cgs.value();
    simulation_r_rg = p_input_reader->simulation_r_rg.value();
  }*/

  // Copy slow-light parameters
  /*if (model_type == ModelType::simulation)
  {
    slow_light_on = p_input_reader->slow_light_on.value();
    if (slow_light_on)
    {
      if (not simulation_multiple)
        throw BlacklightException("Must enable simulation_multiple to use slow light.");
      slow_chunk_size = p_input_reader->slow_chunk_size.value();
      if (slow_chunk_size < 2)
        throw BlacklightException("Must have slow_chunk_size be at least 2.");
      if (slow_chunk_size > simulation_end - simulation_start + 1)
        throw BlacklightException("Not enough simulation files for given slow_chunk_size.");
      slow_t_start = p_input_reader->slow_t_start.value();
      slow_dt = p_input_reader->slow_dt.value();
      if (slow_dt <= 0.0)
        throw BlacklightException("Must have positive time interval slow_dt.");
    }
  }*/

  // Copy plasma parameters
  /*if (model_type == ModelType::simulation)
  {
    plasma_mu = p_input_reader->plasma_mu.value();
    plasma_model = p_input_reader->plasma_model.value();
    if (plasma_model == PlasmaModel::ti_te_beta || plasma_model == PlasmaModel::one_temp)
    {
      plasma_use_p = p_input_reader->plasma_use_p.value();
      if (plasma_use_p)
      {
        if (p_input_reader->plasma_gamma.has_value())
        {
          plasma_gamma = p_input_reader->plasma_gamma.value();
          gamma_set = true;
        }
        if (p_input_reader->plasma_gamma_i.has_value())
          BlacklightWarning("Ignoring plasma_gamma_i selection.");
        if (p_input_reader->plasma_gamma_e.has_value())
          BlacklightWarning("Ignoring plasma_gamma_e selection.");
      }
      else
      {
        if (simulation_format == SimulationFormat::athena
            or p_input_reader->plasma_gamma.has_value())
        {
          plasma_gamma = p_input_reader->plasma_gamma.value();
          gamma_set = true;
        }
        if (simulation_format == SimulationFormat::iharm3d)
        {
          if (p_input_reader->plasma_gamma_i.has_value())
          {
            plasma_gamma_i = p_input_reader->plasma_gamma_i.value();
            gamma_i_set = true;
          }
          if (p_input_reader->plasma_gamma_e.has_value())
          {
            plasma_gamma_e = p_input_reader->plasma_gamma_e.value();
            gamma_e_set = true;
          }
        }
        else
        {
          plasma_gamma_i = p_input_reader->plasma_gamma_i.value();
          plasma_gamma_e = p_input_reader->plasma_gamma_e.value();
        }
      }
    }
  }*/

  // Determine how many files will be held in memory simultaneously
  
}

//--------------------------------------------------------------------------------------------------

// Simulation reader destructor
OpacityTableReader::~OpacityTableReader()
{
  if (num_freqs > 0){
    freq_grid.Deallocate();
  }
  if (num_temps > 0){
    temp_grid.Deallocate();
  }
  if (num_rho > 0){
    rho_grid.Deallocate();
  }
  if (num_freqs > 0 && num_temps > 0 && num_rho > 0){
    ross_tab.Deallocate();
    plan_tab.Deallocate();
  }
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
double OpacityTableReader::Read(int snapshot)
{

  std::cout << "Reading opacity table from file..." << opacity_file << std::endl;
  // Only proceed if needed
  if (model_type != ModelType::simulation)
    return 0.0;
  if(!opacity_table)
    return 0.0;
  double time_start = omp_get_wtime();

  int num_read = 1;

  // Read new files

  // Determine file name
  std::string opacity_file_formatted = opacity_file;
  std::cout << "Formatted opacity file name: " << opacity_file_formatted << std::endl;

  // Open input file
  if ( (opac_file=fopen(opacity_file_formatted.c_str(),"r"))==NULL) {
    std::stringstream msg;
    msg << "FATAL ERROR: Could not open " << opacity_file_formatted << "." << std::endl;
    throw BlacklightException(msg.str().c_str());
  }
  // Read basic data about file
  fscanf(opac_file,"%d",&(num_freqs));
  fscanf(opac_file,"%d",&(num_temps));
  fscanf(opac_file,"%d",&(num_rho));
  
  // allocate arrays
  freq_grid.Allocate(num_freqs+1);
  temp_grid.Allocate(num_temps);
  rho_grid.Allocate(num_rho);
  ross_tab.Allocate(num_freqs,num_temps,num_rho);
  plan_tab.Allocate(num_freqs,num_temps,num_rho);

  for(int i=1; i<=num_freqs; ++i){
    fscanf(opac_file,"%lf",&(freq_grid(i)));
  }

  freq_grid(0) = freq_grid(1)*freq_grid(1)/freq_grid(2);
  // convert to erg
  double keverg = 1.602176634e-9;
  for(int i=0; i<=num_freqs; ++i)
    freq_grid(i) *= keverg;
  fmin = std::log10(freq_grid(1));
  fmax = std::log10(freq_grid(num_freqs-1));
  dlf = (fmax-fmin)/(num_freqs-2);

  // temperature grid (keV)
  for(int i=0; i<num_temps; ++i){
    fscanf(opac_file,"%lf",&(temp_grid(i)));
  }

  // convert to kelvin
  for(int i=0; i<num_temps; ++i)
    temp_grid(i) *= keverg/Physics::k_b;
  tmin = std::log10(temp_grid(0));
  tmax = std::log10(temp_grid(num_temps-1));
  dlt = (tmax-tmin)/(num_temps-1);
  // note temperature grid not uniform in log

  // density grid (g/cm^3)
  for(int i=0; i<num_rho; ++i) {
    fscanf(opac_file,"%lf",&(rho_grid(i)));
  }
  rmin = std::log10(rho_grid(0));
  rmax = std::log10(rho_grid(num_rho-1));
  dlr = (rmax-rmin)/(num_rho-1);


  std::cout<<"Opacity table frequency range from "<<freq_grid(0)/Physics::h<<" to "<<freq_grid(num_freqs)/Physics::h<< " Hz with "<<num_freqs<<" frequencies in total."<<std::endl;
  std::cout<<"Opacity table temperature range from "<<temp_grid(0)<<" to "<<temp_grid(num_temps-1)<<" K with "<<num_temps<<" temperatures in total."<<std::endl;
  std::cout<<"Opacity table density range from "<<rho_grid(0)<<" to "<<rho_grid(num_rho-1)<<" g/cm^3 with "<<num_rho<<" densities in total."<<std::endl;


  // frequency integrated rosseland mean
  // Read in but not used
  double buf;
  for(int j=0; j<num_temps; ++j) {
    for(int i=0; i<num_rho; ++i) {
      fscanf(opac_file,"%lf",&buf);
    }
  }

  // frequency integrated planck mean
  // Read in but not used
  for(int j=0; j<num_temps; ++j) {
    for(int i=0; i<num_rho; ++i) {
      fscanf(opac_file,"%lf",&buf);
    }
  }

  // ross mean for each frequency group
  for(int k=0; k<num_freqs; ++k) {
    for(int j=0; j<num_temps; ++j) {
      for(int i=0; i<num_rho; ++i) {
        fscanf(opac_file,"%lf",&(ross_tab(k,j,i)));
        ross_tab(k,j,i) *= rho_grid(i);
        //std::cout << "Rosseland opacity " << k << " " << j << " " << i << " : " << ross_tab(k,j,i) << std::endl;
      }
    }
  }

  // planck mean for each frequency group
  for(int k=0; k<num_freqs; ++k) {
    for(int j=0; j<num_temps; ++j) {
      for(int i=0; i<num_rho; ++i) {
        fscanf(opac_file,"%lf",&(plan_tab(k,j,i)));
        plan_tab(k,j,i) *= rho_grid(i);
        //std::cout << "Planck opacity " << k << " " << j << " " << i << " : " << plan_tab(k,j,i) << std::endl;
      }
    }
  }

  // planck mean for each frequency group
  for(int j=0; j<num_temps; ++j) {
    for(int i=0; i<num_rho; ++i) {
      double min = 1.e40;
      double max = 1.e-40;
      for(int k=0; k<num_freqs; ++k) {
        min = (min > plan_tab(k,j,i)) ? plan_tab(k,j,i) : min;
        max = (max < plan_tab(k,j,i)) ? plan_tab(k,j,i) : max;
      }
      if (max/min < 1.2) {
        double ffnrm = 3.692146e8;
        double heabund = 0.09; //hardcode for now (should be parameter)
        double tgas = temp_grid(j);
        double nh = rho_grid(i) / (Physics::m_p*(1.+4.*heabund));
        double nhe = nh*heabund;
        double ne = nh + 2.*nhe;
        for(int k=0; k<num_freqs; ++k) {
          double nu = freq_grid(k) / Physics::h;
          double ehnu = exp(-Physics::h*nu / (Physics::k_b * tgas) );
          double aff = ffnrm/sqrt(tgas)/pow(nu,3);
          double opac = ne * (nh + 4. * nhe) * aff * (1. - ehnu);
          //printf("%d %g %g %g %g\n",k,temp_grid(i),rho_grid(j),plan_tab(k,j,i),opac);
          plan_tab(k,j,i) = opac;
        }
      }
    }
  }
  
  
    // Read time
    

      // Convert internal energy to pressure
      /*#pragma omp parallel for schedule(static) collapse(3)
      for (int block = 0; block < athenak_num_blocks; block++)
        for (int k = 0; k < athenak_block_nz; k++)
          for (int j = 0; j < athenak_block_ny; j++)
            for (int i = 0; i < athenak_block_nx; i++)
              prim[n](ind_pgas,block,k,j,i) *= static_cast<float>(plasma_gamma - 1.0);
    */


    

  // Close input file
  fclose(opac_file);

  // Calculate elapsed time
  return omp_get_wtime() - time_start;

}

void OpacityTableReader::InterpolateToUniformTLog(){
  //have a holder for the original, non-uniform temp grid so that we can refer to it
  Array<double> temp_grid_holder = temp_grid;
  int temp_indices[num_temps];
  int num_found = 0;

  temp_indices[0] = 0;
  //we start at the second index because we know the first and the last must be the same btwn both grids
  for(int i=1; i<num_temps; ++i){
    //update the true temp_grid to be uniform in log
    temp_grid(i) = pow(10.0, tmin + i * dlt);

    //identify which two points in the non-uniform grid you are between by performing binary search
    temp_indices[i] = temp_indices[i-1]+1;
    int high = num_temps - 1;
    while (temp_indices[i] < high) {
        int mid = temp_indices[i] + (high - temp_indices[i]) / 2;
        if (temp_grid_holder(mid)< temp_grid(i) && temp_grid_holder(mid+1) >= temp_grid(i)) {
            temp_indices[i] = mid;
            break;
        }
        else if (temp_grid_holder(mid) < temp_grid(i)) {
            temp_indices[i] = mid + 1;
        }
        else {
            high = mid - 1;
        }
    }

    std::cout<<temp_grid(i)<<" vs "<<temp_grid_holder(i)<<std::endl;

    //now that you have your upper and lower t bounds, interpolate the opacities to the uniform log T grid
    for(int f=0;f<num_freqs;++f){
      for(int r=0;r<num_rho;++r){
        double plan_opacity_low = plan_tab(f,temp_indices[i],r);
        double ross_opacity;

        //TEGAN: TO DO 

      }
    }
  }

  for(int t=0;t<num_temps;++t){

    
  }

  




}