// Blacklight simulation reader header

#ifndef OPACITY_TABLE_READER_H_
#define OPACITY_TABLE_READER_H_

// C++ headers
#include <fstream>  // ifstream
#include <iosfwd>   // streampos
#include <string>   // string

// Blacklight headers
#include "../blacklight.hpp"                 // enums
#include "../input_reader/input_reader.hpp"  // InputReader
#include "../utils/array.hpp"                // Array

//--------------------------------------------------------------------------------------------------

// OpacityTable reader
struct OpacityTableReader
{
  // Constructors and destructor
  OpacityTableReader(const InputReader *p_input_reader_);
  OpacityTableReader(const OpacityTableReader &source) = delete;
  OpacityTableReader &operator=(const OpacityTableReader &source) = delete;
  ~OpacityTableReader();

  // Pointers to other objects
  const InputReader *p_input_reader;

  // Input data - general
  ModelType model_type;

  // Input data - opacity parameters
  std::string opacity_file;
  bool opacity_table;

  
  FILE  *opac_file;

  // Data arrays
  int num_freqs;
  int num_temps;
  int num_rho;

  Array<double> freq_grid;
  Array<double> temp_grid;
  Array<double> rho_grid;

  Array<double> ross_tab;
  Array<double> plan_tab;

  double fmin, fmax, dlf;
  double tmin, tmax, dlt;
  double rmin, rmax, dlr;
  
  // External function
  double Read(int snapshot);
  void InterpolateToUniformTLog();

  // Internal functions - opacity_reader.cpp
  //std::string FormatFilename(int file_number);

};

#endif
