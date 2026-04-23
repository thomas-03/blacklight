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
  mc_file_name = p_input_reader->mc_file.value();
  mc_freq_file_name = p_input_reader->mc_freq_file.value();
  simulation_coord = p_input_reader->simulation_coord.value();
  simulation_m_msun = p_input_reader->simulation_m_msun.value();
  simulation_all_cgs = p_input_reader->simulation_all_cgs.value();
    if(simulation_all_cgs){
      simulation_r_rg = Physics::c*Physics::c/(simulation_m_msun*Physics::gg_msun);
    }else{
      simulation_r_rg = p_input_reader->simulation_r_rg.value();
    }
  simulation_format = p_input_reader->simulation_format.value();
  scattering_source_terms = new Array<float>[1];

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
  std::cout<<"read freq file"<<std::endl;

  // Read new files

  // Determine file name

  // Open input file
  /*if ( (mc_file=fopen(mc_file_name.c_str(),"r"))==NULL) {
    std::stringstream msg;
    msg << "FATAL ERROR: Could not open " << mc_file_name << "." << std::endl;
    throw BlacklightException(msg.str().c_str());
  }*/

  data_stream =
        std::ifstream(mc_file_name, std::ios_base::in | std::ios_base::binary);
  if (not data_stream.is_open())
      throw BlacklightException("Could not open MC file for reading.");

  if (simulation_format == SimulationFormat::athenak or simulation_format == SimulationFormat::iharm3d)
    {
      throw BlacklightException("Reading MC files is only implemented for Athena++ formats.");
    }
  
  // Read block layout
   /* if (first_time)
    {
      if (simulation_format == SimulationFormat::athena)
      {
        SimulationReader::ReadHDF5IntArray("Levels", levels);
        SimulationReader::ReadHDF5IntArray("LogicalLocations", locations);
      }
      else if (simulation_format == SimulationFormat::iharm3d
          or simulation_format == SimulationFormat::harm3d)
      {
        levels.Allocate(1);
        levels(0) = 0;
        locations.Allocate(1, 3);
        locations(0,0) = 0;
        locations(0,1) = 0;
        locations(0,2) = 0;
      }
    }*/

    // Read coordinates
    if (first_time)
    {
      std::cout<<"in first time reading"<<std::endl;
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
    
    if (first_time)
      {
        //literally when i change the n1's here to n2's all of a sudden I get a "could not read MC file" from above
        //std::cout<<num_freqs<<","<<levels.n1<<","<<x3v.n1<<","<<x2v.n1<<","<<x1v.n1<<std::endl;
        //want it to be frequency, b, k ,j, i
        int n5 = num_freqs;
        //int n4 = levels.n1; // 256
        int n3 = x3v.n1; 
        int n2 = x2v.n1; 
        int n1 = x1v.n1;
        int n4 = levels.n1;
        //int n5 = phiBins;
        scattering_source_terms[0].Allocate(n5,n4, n3, n2, n1);
      }
      //Array<float> source_terms(scattering_source_terms[0]);

      Array<float> shallow_scatter(scattering_source_terms[0]);
       //std::cout<<"number of frequencies: "<<num_freqs<<std::endl;

      //ReadHDF5FloatArray("prim", hydro);
      //scattering_source_terms[0].Slice(5, 0, num_freqs-1);
      //ReadHDF5FloatArray("mcscat",scattering_source_terms[0]);
      std::cout<<"before read flaot array"<<std::endl;
      ReadHDF5FloatArray("mcscat",shallow_scatter);
      std::cout<<"after read it"<<std::endl;
      //(20, 256, 8, 8, 32)
      /*for(int i=0;i<num_freqs;i++){
          std::cout<<scattering_source_terms[0](i,0,0,0,0)<<", ";
      }*/

    // Close input file
    data_stream.close();

    // Update first time flag
    first_time = false;

  std::cout<<"finished MC"<<std::endl;
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