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
  //std::cout << "Formatted opacity file name: " << opacity_file_formatted << std::endl;

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
  freq_grid.Allocate(num_freqs);
  temp_grid.Allocate(num_temps);
  rho_grid.Allocate(num_rho);
  ross_tab.Allocate(num_freqs,num_temps,num_rho);
  plan_tab.Allocate(num_freqs,num_temps,num_rho);

  for(int i=0; i<num_freqs; ++i){
    fscanf(opac_file,"%lf",&(freq_grid(i)));
  }

  // convert to erg
  double keverg = 1.602176634e-9;
  for(int i=0; i<num_freqs; ++i)
    freq_grid(i) *= keverg;
  fmin = std::log10(freq_grid(0));
  fmax = std::log10(freq_grid(num_freqs-1));
  dlf = (fmax-fmin)/(num_freqs-1);

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


  std::cout<<"Opacity table frequency range from "<<freq_grid(0)/Physics::h<<" to "<<freq_grid(num_freqs-1)/Physics::h<< " Hz with "<<num_freqs<<" frequencies in total."<<std::endl;
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
  //std::ofstream kTFile;
  //kTFile.open("./opacityTableOverwritten.csv", std::ios_base::app);
  //kTFile<<"table value"<<","<<"free-free value"<<","<<"ratio"<<","<<"frequency(Hz)"<<","<<"temperature(K)"<<","<<"density(g/cm3)"<<std::endl;
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
        for(int k=0; k<num_freqs; ++k) {
          double partA = 4*pow(Physics::e,6.)/(3*Physics::m_e*Physics::c*Physics::h);
          double partB = std::sqrt(2.0*Math::pi/(3.0*Physics::k_b*temp_grid(j)*Physics::m_e));
          double gaunt_factor = 1.0; //approximate it as this because shouldn't impact too much
          double nu = freq_grid(k) / Physics::h;
          //note that I hardcoded in the 0.5 mean molecular weight here
          double coefficient = partA*partB*(rho_grid(i)/(0.5*Physics::m_p))*(rho_grid(i)/(0.5*Physics::m_p))*(1.0 - std::exp(-Physics::h*nu/ (Physics::k_b * temp_grid(j))))*gaunt_factor/(nu*nu*nu);
          //kTFile<<plan_tab(k,j,i)<<","<<coefficient<<","<<plan_tab(k,j,i)/coefficient<< ","<<nu<<","<<temp_grid(j)<<","<<rho_grid(i)<<std::endl;
          plan_tab(k,j,i) = coefficient;
        }
      }/*else{
        //manually check if the tables always deviated from free free or if it's a new thing

        for(int k=0; k<num_freqs; ++k) {
          double partA = 4*pow(Physics::e,6.)/(3*Physics::m_e*Physics::c*Physics::h);
          double partB = std::sqrt(2.0*Math::pi/(3.0*Physics::k_b*temp_grid(j)*Physics::m_e));
          double gaunt_factor = 1.0; //approximate it as this because shouldn't impact too much
          double nu = freq_grid(k) / Physics::h;
          double coefficient = partA*partB*(rho_grid(i)/(0.6*Physics::m_p))*(rho_grid(i)/(0.6*Physics::m_p))*(1.0 - std::exp(-Physics::h*nu/ (Physics::k_b * temp_grid(j))))*gaunt_factor/(nu*nu*nu);
          kTFile<<plan_tab(k,j,i)<<","<<coefficient<<","<<plan_tab(k,j,i)/coefficient<< ","<<nu<<","<<temp_grid(j)<<","<<rho_grid(i)<<std::endl;
        }
      }*/
     //my version of the free-free opacity calculation
     /*for(int k=0; k<num_freqs; ++k) {
           double partA = 4*pow(Physics::e,6.)/(3*Physics::m_e*Physics::c*Physics::h);
           double partB = std::sqrt(2.0*Math::pi/(3.0*Physics::k_b*temp_grid(j)*Physics::m_e));
           double gaunt_factor = 1.0; //approximate it as this because shouldn't impact too much
           //double n_cgs = rho_cgs / (plasma_mu * Physics::m_p);
           double nu = freq_grid(k) / Physics::h;
           double coefficient = partA*partB*(rho_grid(i)/(0.6*Physics::m_p))*(rho_grid(i)/(0.6*Physics::m_p))*(1.0 - std::exp(-Physics::h*nu/ (Physics::k_b * temp_grid(j))))*gaunt_factor/(nu*nu*nu);
           plan_tab(k,j,i) = coefficient;
      }*/
    }
  }
  //kTFile.close();
  
  //interpolate to uniform log T grid
  InterpolateToUniformTLog();
  
  // Close input file
  fclose(opac_file);

  // Calculate elapsed time
  return omp_get_wtime() - time_start;

}

void OpacityTableReader::InterpolateToUniformTLog(){
  //have a holder for the original, non-uniform temp grid so that we can refer to it
  Array<double> temp_grid_holder;
  temp_grid_holder.Allocate(num_temps);
  temp_grid_holder.CopyFrom(temp_grid,0,0,num_temps);
  Array<double> plan_tab_holder;
  plan_tab_holder.Allocate(num_freqs,num_temps,num_rho);
  plan_tab_holder.CopyFrom(plan_tab,0,0,plan_tab.GetNumBytes()/sizeof(double));
  Array<double> ross_tab_holder;
  ross_tab_holder.Allocate(num_freqs,num_temps,num_rho);
  ross_tab_holder.CopyFrom(ross_tab,0,0,ross_tab.GetNumBytes()/sizeof(double));
  int temp_indices[num_temps];
  int num_found = 0;

  temp_indices[0] = 0;

  //we start at the second index because we know the first and the last must be the same btwn both grids
  for(int i=1; i<num_temps-1; ++i){
    //update the true temp_grid to be uniform in log
    temp_grid(i) = pow(10.0, tmin + i * dlt);

    //identify which two points in the non-uniform grid you are between by performing binary searc
    int low = temp_indices[i-1];
    temp_indices[i] = low;
    int high = num_temps - 2;
    
    while (low <= high) {
        int mid = low + std::floor((high - low) / 2);
        
        if (temp_grid_holder(mid)< temp_grid(i) && temp_grid_holder(mid+1) >= temp_grid(i)) {
            temp_indices[i] = mid;
            break;
        }
        else if (temp_grid_holder(mid+1) <= temp_grid(i)) {
            low = mid + 1;
        }
        else {
            high = mid - 1;
        }
    }

    //now that you have your upper and lower t bounds, interpolate the opacities to the uniform log T grid
    for(int f=0;f<num_freqs;++f){
      for(int r=0;r<num_rho;++r){
        double plan_opacity_low = plan_tab_holder(f,temp_indices[i],r);
        double ross_opacity_low = ross_tab_holder(f,temp_indices[i],r);
        double plan_opacity_high = plan_tab_holder(f,temp_indices[i]+1,r);
        double ross_opacity_high = ross_tab_holder(f,temp_indices[i]+1,r);

        double t_low = temp_grid_holder(temp_indices[i]);
        double t_high = temp_grid_holder(temp_indices[i]+1);
        double t_target = temp_grid(i);

        plan_tab(f,i,r) = plan_opacity_low + (plan_opacity_high - plan_opacity_low)*(t_target - t_low)/(t_high - t_low);
        ross_tab(f,i,r) = ross_opacity_low + (ross_opacity_high - ross_opacity_low)*(t_target - t_low)/(t_high - t_low);

        if(plan_tab(f,i,r)<0){
          std::cout<<"negative opacity. prev value "<<plan_opacity_low<< " next value "<<plan_opacity_high<<" indice "<<temp_indices[i]<<std::endl;
        }
      }
    }
  }

  temp_grid_holder.Deallocate();
  plan_tab_holder.Deallocate();
  ross_tab_holder.Deallocate();
}