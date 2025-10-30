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
  bool use_opacity_table;

  // Input data - slow-light parameters
  bool slow_light_on;
  int slow_chunk_size;
  double slow_t_start;
  double slow_dt;

  // Input data - plasma parameters
  double plasma_mu;
  PlasmaModel plasma_model;
  bool plasma_use_p;
  double plasma_gamma;
  double plasma_gamma_i;
  double plasma_gamma_e;

  // Flags for tracking function calls
  bool first_time = true;
  bool first_time_root_object_header = true;

  // Metadata


  std::string metric;
  double metric_a, metric_h, metric_r_in;
  double metric_poly_xt, metric_poly_alpha, metric_mks_smooth, metric_derived_poly_norm;
  std::ifstream::pos_type cell_data_address;
  std::string *dataset_names;
  int num_dataset_names = 0;
  std::string *variable_names;
  int num_variable_names = 0;
  Array<int> num_variables;
  int ind_hydro;
  int ind_bb;
  int ind_rho, ind_pgas, ind_kappa;
  int ind_u0, ind_uu1, ind_uu2, ind_uu3;
  int ind_b0, ind_bb1, ind_bb2, ind_bb3;
  int ind_rad;
  int num_arrays;
  int latest_file_number;
  const double extrapolation_tolerance = 1.0;
  const double angular_domain_tolerance = 0.1;
  bool gamma_set = false;
  bool gamma_i_set = false;
  bool gamma_e_set = false;

  // Coordinate interpolation data
  double sks_map_r_in, sks_map_r_out, sks_map_dr, sks_map_dtheta;
  Array<double> opacity_bounds;
  Array<double> sks_map;
  const int sks_map_n1 = 2048;
  const int sks_map_n2 = 2048;
  const int sks_map_max_iter = 1000;
  const double sks_map_tol = 1.0e-8;


  FILE  *opac_file;

  // Data


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

  // Internal functions - opacity_reader.cpp
  //std::string FormatFilename(int file_number);

  // Internal functions - hdf5_format_arrays.cpp
  void ReadHDF5StringArray(const char *name, bool allocate, std::string **p_string_array,
      int *p_array_length);
  void ReadHDF5IntArray(const char *name, Array<int> &int_array);
  void ReadHDF5FloatArray(const char *name, Array<float> &float_array);
  void ReadHDF5FloatArray(const char *name, Array<double> &double_array);
  void ReadHDF5DoubleArray(const char *name, Array<double> &double_array);
  static void SetHDF5StringArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, bool allocate,
      std::string **string_array, int *p_array_length);
  static void SetHDF5IntArray(const unsigned char *datatype_raw, const unsigned char *dataspace_raw,
      const unsigned char *data_raw, Array<int> &int_array);
  static void SetHDF5FloatArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, Array<float> &float_array);
  static void SetHDF5FloatArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw,
      Array<double> &double_array);
  static void SetHDF5DoubleArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw,
      Array<double> &double_array);
};

#endif
