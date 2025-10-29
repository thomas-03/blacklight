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
  num_arrays = 0;
  if (model_type == ModelType::simulation)
    num_arrays = slow_light_on ? slow_chunk_size : 1;

}

//--------------------------------------------------------------------------------------------------

// Simulation reader destructor
OpacityTableReader::~OpacityTableReader()
{
  if (num_dataset_names > 0)
    delete[] dataset_names;
  if (num_variable_names > 0)
    delete[] variable_names;
}

//--------------------------------------------------------------------------------------------------

// Simulation reader read and initialize function
// Inputs:
//   snapshot: index (starting at 0) of which snapshot is about to be prepared
// Outputs:
//   returned value: execution time in seconds
// Notes:
//   Does nothing if model does not need to be read from file.
//   The output file offset is always equal to snapshot; the input file offset is equal to snapshot
//       if slow_light_on == false.
//   Opens and closes stream for reading.
//   Initializes all member objects.
//   Implements a subset of the HDF5 standard:
//       portal.hdfgroup.org/display/HDF5/File+Format+Specification
double OpacityTableReader::Read(int snapshot)
{
  // Only proceed if needed
  if (model_type != ModelType::simulation)
    return 0.0;
  double time_start = omp_get_wtime();

  int num_read = 1;

  // Read new files

  // Determine file name
  std::string opacity_file_formatted = opacity_file;
  opacity_file_formatted = FormatFilename(-1);

  // Open input file
  if ( (opac_file=fopen("out_opacity_table_num_freqsq32.txt","r"))==NULL) {
    std::stringstream msg;
    msg << "FATAL ERROR: Could not open out_opacity_table_num_freqsq16.txt." << std::endl;
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
      }
    }
  }

  // planck mean for each frequency group
  for(int k=0; k<num_freqs; ++k) {
    for(int j=0; j<num_temps; ++j) {
      for(int i=0; i<num_rho; ++i) {
        fscanf(opac_file,"%lf",&(plan_tab(k,j,i)));
        plan_tab(k,j,i) *= rho_grid(i);
      }
    }
  }

  // Replaces plan_tab with free-free values
  /*Real dummy;
  for(int k=0; k<num_freqs; ++k) {
    Real ffnrm = 3.692146e8;
    Real heabund = 0.09; //hardcode for now (should be parameter)
    Real mp = 1.67262192369e-24;
    Real h = 6.62607015e-27;
    Real kb = 1.380649e-16;
    Real nu = freq_grid(k) / h;
    for(int j=0; j<num_temps; ++j) {
      Real tgas = temp_grid(j);
      Real ehnu = exp(-h*nu / (kb * tgas) );
      for(int i=0; i<num_rho; ++i) {
        Real nh = rho_grid(i) / (mp*(1.+4.*heabund));
        Real nhe = nh*heabund;
        Real ne = nh + 2.*nhe;
        fscanf(opac_file,"%lf",&(dummy));
        Real aff = ffnrm/sqrt(tgas)/pow(nu,3);
        Real opac = ne * (nh + 4. * nhe) * aff * (1. - ehnu);
        plan_tab(k,j,i) = opac;
      }
    }
    }*/

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

//--------------------------------------------------------------------------------------------------

// Function to construct filename formatted with file number
// Inputs:
//   file_number: number of simulation file to construct
// Outputs:
//   returned_value: string containing formatted filename
std::string OpacityTableReader::FormatFilename(int file_number)
{
  // Locate braces
  std::string::size_type opacity_pos_open = opacity_file.find_first_of('{');
  std::string::size_type opacity_pos_close = opacity_file.find_first_of('}', opacity_pos_open);

  // Parse integer format string
  int opacity_field_length = 0;
  if (opacity_pos_close - opacity_pos_open > 2)
    opacity_field_length = std::stoi(opacity_file.substr(opacity_pos_open + 1,
        opacity_pos_close - opacity_pos_open - 2));
  int file_number_length = std::snprintf(nullptr, 0, "%d", file_number);
  if (file_number_length < 0)
    throw BlacklightException("Could not format file name.");
  int num_zeros = 0;
  if (file_number_length < opacity_field_length)
    num_zeros = opacity_field_length - file_number_length;

  // Create filename
  std::ostringstream opacity_filename;
  opacity_filename << opacity_file.substr(0, opacity_pos_open);
  for (int n = 0; n < num_zeros; n++)
    opacity_filename << "0";
  opacity_filename << file_number;
  opacity_filename << opacity_file.substr(opacity_pos_close + 1);
  std::string opacity_file_formatted = opacity_filename.str();
  return opacity_file_formatted;
}

