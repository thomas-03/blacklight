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

// Function to read HDF5 superblock
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Indirectly sets root_object_header_address, root_btree_address, and root_name_heap_address.
//   Changes stream pointer.
//   Does not allow for superblock to be anywhere but beginning of file.
//   Must have superblock version 0.
//   Must have size of offsets 8.
//   Must have size of lengths 8.
void MCReader::ReadHDF5Superblock()
{
  // Check format signature
  data_stream.seekg(0);
  const unsigned char expected_signature[] = {0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a};
  for (int n = 0; n < 8; n++)
    if (data_stream.get() != expected_signature[n])
      throw BlacklightException("Unexpected HDF5 format signature.");

  // Check superblock version
  int superblock_vers = data_stream.get();
  if (superblock_vers != 0 && superblock_vers != 1){
    std::printf("superblock_vers: %d ",superblock_vers);
    ReadHDF5SuperblockVers2();
    return;
  }
  //  throw BlacklightException("Unexpected HDF5 superblock version.");

  // Check other version numbers
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 file free space storage version.");
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 root group symbol table entry version.");
  data_stream.ignore(1);
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 shared header message format version.");

  // Check sizes
  if (data_stream.get() != 8)
    throw BlacklightException("Unexpected HDF5 size of offsets.");
  if (data_stream.get() != 8)
    throw BlacklightException("Unexpected HDF5 size of lengths.");
  data_stream.ignore(1);

  // Skip checking tree parameters and consistency flags
  data_stream.ignore(2 * 2 + 4);

  // Skip Indexed Storage Internal Node K for superblock version 1
  if(superblock_vers==1){
    data_stream.ignore(2*2);
  }

  // Skip checking addresses
  data_stream.ignore(4 * 8);

  // Read root group symbol table entry
  ReadHDF5RootGroupSymbolTableEntry();
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 superblock
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Indirectly sets root_object_header_address, root_btree_address, and root_name_heap_address.
//   Changes stream pointer.
//   Does not allow for superblock to be anywhere but beginning of file.
//   Must have superblock version 0.
//   Must have size of offsets 8.
//   Must have size of lengths 8.
void MCReader::ReadHDF5SuperblockVers2()
{
  // Check sizes
  if (data_stream.get() != 8)
    throw BlacklightException("Unexpected HDF5 size of offsets.");
  if (data_stream.get() != 8)
    throw BlacklightException("Unexpected HDF5 size of lengths.");
  data_stream.ignore(1);

  // Skip checking addresses
  data_stream.ignore(3 * 8);

  // Read object header address
  data_stream.read(reinterpret_cast<char *>(&root_object_header_address), 8);
  ReadHDF5RootObjectHeader();
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 root group symbol table entry
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Sets root_object_header_address, root_btree_address, and root_name_heap_address.
//   Assumes stream pointer is already set.
//   Must have size of offsets 8.
//   Must be run on little-endian machine.
void MCReader::ReadHDF5RootGroupSymbolTableEntry()
{
  // Skip reading link name offset
  data_stream.ignore(8);

  // Read object header address
  data_stream.read(reinterpret_cast<char *>(&root_object_header_address), 8);

  // Check cache type
  unsigned int cache_type;
  data_stream.read(reinterpret_cast<char *>(&cache_type), 4);
  if (cache_type != 1)
    throw BlacklightException("Unexpected HDF5 root group symbol table entry cache type.");
  data_stream.ignore(4);

  // Read B-tree and name heap addresses
  data_stream.read(reinterpret_cast<char *>(&root_btree_address), 8);
  data_stream.read(reinterpret_cast<char *>(&root_name_heap_address), 8);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 local heap
// Inputs:
//   heap_address: global address of start of heap
// Outputs:
//   returned value: global address of start of heap data segment
// Notes:
//   Sets root_data_segment_address.
//   Changes stream pointer.
//   Must have size of offsets 8.
unsigned long int MCReader::ReadHDF5Heap(unsigned long int heap_address)
{
  // Check local heap signature and version
  data_stream.seekg(static_cast<std::streamoff>(heap_address));
  const unsigned char expected_signature[] = {'H', 'E', 'A', 'P'};
  for (int n = 0; n < 4; n++)
    if (data_stream.get() != expected_signature[n])
      throw BlacklightException("Unexpected HDF5 heap signature.");
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 heap version.");
  data_stream.ignore(3);

  // Skip data segement size and offset to head of free list
  data_stream.ignore(16);

  // Read address of data segment
  unsigned long int data_segment_address;
  data_stream.read(reinterpret_cast<char *>(&data_segment_address), 8);
  return data_segment_address;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 root object header
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Allocates and sets dataset_names, variable_names, nums_variables.
//   Sets num_dataset_names and num_variable_names.
//   Assumes root_object_header_address set.
//   Changes stream pointer.
//   Must have object header version 1.
//   Must not have shared header messages.
//   Must have attribute message version 1.
//   Must have size of offsets 8.
//   Must be run on little-endian machine.
void MCReader::ReadHDF5RootObjectHeader()
{
  if(!first_time_root_object_header){
    return;
  }

  // Check object header version
  data_stream.seekg(static_cast<std::streamoff>(root_object_header_address));
  if (data_stream.get() != 1){
    ReadHDF5RootObjectHeaderVers2();
    return;
  }
  data_stream.ignore(1);

  // Read number of header messages
  unsigned short int num_messages;
  data_stream.read(reinterpret_cast<char *>(&num_messages), 2);

  // Skip reading object reference count and object header size
  data_stream.ignore(8);

  // Align to 8 bytes within header (location of padding not documented)
  data_stream.ignore(4);

  // Go through messages
  bool root_grid_size_found = false;
  bool dataset_names_found = false;
  bool variable_names_found = false;
  bool num_variables_found = false;
  for (int n = 0; n < num_messages; n++)
  {
    // Read message type and size
    unsigned short int message_type, message_size;
    data_stream.read(reinterpret_cast<char *>(&message_type), 2);
    data_stream.read(reinterpret_cast<char *>(&message_size), 2);

    // Check message flags
    unsigned char message_flags;
    data_stream.read(reinterpret_cast<char *>(&message_flags), 1);
    data_stream.ignore(3);
    if (message_flags & 0b00000010)
      throw BlacklightException("Unexpected HDF5 header message flag.");

    // Read message data
    unsigned char *message_data = new unsigned char[message_size];
    data_stream.read(reinterpret_cast<char *>(message_data), message_size);

    // Follow any continuation messages
    if (message_type == 16)
    {
      unsigned long int new_offset;
      std::memcpy(&new_offset, message_data, 8);
      data_stream.seekg(static_cast<std::streamoff>(new_offset));
      continue;
    }

    // Inspect any attribute messages
    else if (message_type == 12)
    {

      // Check attribute message version
      int offset = 0;
      if (message_data[offset] != 1)
        throw BlacklightException("Unexpected HDF5 attribute message version.");
      offset += 2;

      // Read attribute message metadata
      unsigned short int name_size, datatype_size, dataspace_size;
      std::memcpy(&name_size, message_data + offset, 2);
      offset += 2;
      std::memcpy(&datatype_size, message_data + offset, 2);
      offset += 2;
      std::memcpy(&dataspace_size, message_data + offset, 2);
      offset += 2;
      unsigned short int name_size_pad = static_cast<unsigned short int>((8 - name_size % 8) % 8);
      unsigned short int datatype_size_pad =
          static_cast<unsigned short int>((8 - datatype_size % 8) % 8);
      unsigned short int dataspace_size_pad =
          static_cast<unsigned short int>((8 - dataspace_size % 8) % 8);

      // Read attribute message data
      unsigned char *name_raw = new unsigned char[name_size];
      unsigned char *datatype_raw = new unsigned char[datatype_size];
      unsigned char *dataspace_raw = new unsigned char[dataspace_size];
      std::memcpy(name_raw, message_data + offset, name_size);
      offset += name_size + name_size_pad;
      std::memcpy(datatype_raw, message_data + offset, datatype_size);
      offset += datatype_size + datatype_size_pad;
      std::memcpy(dataspace_raw, message_data + offset, dataspace_size);
      offset += dataspace_size + dataspace_size_pad;
      std::string name(reinterpret_cast<char *>(name_raw),
          static_cast<std::string::size_type>(static_cast<int>(name_size) - 1));

      // Read and set desired attributes
      if (name == "RootGridSize")
      {
        root_grid_size_found = true;
        Array<int> root_grid_size;
        SetHDF5IntArray(datatype_raw, dataspace_raw, message_data + offset, root_grid_size);
        n_3_root = root_grid_size(2);
      }
      else if (name == "DatasetNames")
      {
        dataset_names_found = true;
        SetHDF5StringArray(datatype_raw, dataspace_raw, message_data + offset,
            first_time_root_object_header, &dataset_names, &num_dataset_names);
      }
      else if (name == "VariableNames")
      {
        variable_names_found = true;
        SetHDF5StringArray(datatype_raw, dataspace_raw, message_data + offset,
            first_time_root_object_header, &variable_names, &num_variable_names);

      }
      else if (name == "NumVariables")
      {
        num_variables_found = true;
        SetHDF5IntArray(datatype_raw, dataspace_raw, message_data + offset, num_variables);
      }

      // Free raw buffers
      delete[] name_raw;
      delete[] datatype_raw;
      delete[] dataspace_raw;
    }

    else if(message_type == 17){
      data_stream.read(reinterpret_cast<char *>(&root_btree_address), 8);
      data_stream.read(reinterpret_cast<char *>(&root_name_heap_address), 8);
    }

    // Free raw buffer
    delete[] message_data;

    // Break when required information found
    if (root_grid_size_found and dataset_names_found and variable_names_found
        and num_variables_found)
      break;
  }

  // Check that appropriate messages were found
  if (not (root_grid_size_found and dataset_names_found and variable_names_found
      and num_variables_found)){
        std::printf("root_grid_size_found:%d  dataset_names_found: %d  variable_names_found: %d num_variables_found: %d \n ",root_grid_size_found,dataset_names_found,variable_names_found,num_variables_found);
        throw BlacklightException("Could not find needed file-level attributes.");
      }
  if (num_variables.n1 != num_dataset_names)
    throw BlacklightException("DatasetNames and NumVariables file-level attribute mismatch.");

  // Update first time flag
  first_time_root_object_header = false;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 root object header
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Allocates and sets dataset_names, variable_names, nums_variables.
//   Sets num_dataset_names and num_variable_names.
//   Assumes root_object_header_address set.
//   Changes stream pointer.
//   Must have object header version 1.
//   Must not have shared header messages.
//   Must have attribute message version 1.
//   Must have size of offsets 8.
//   Must be run on little-endian machine.
void MCReader::ReadHDF5RootObjectHeaderVers2()
{
  std::printf("root object header version 2 ");
  //the true signature is 4 bytes long but 1 is already read when trying to check version
  unsigned char *signature = new unsigned char[3];
    data_stream.read(reinterpret_cast<char *>(signature), 3);
  if(signature[0]!=72 || signature[1]!=68 || signature[2]!=82)
    throw BlacklightException("Unexpected HDF5 Root Object Version 2 Header");

  // Check header version
  if (data_stream.get() != 2)
    throw BlacklightException("Unexpected HDF5 Root Object Header version");
  
  // Read Flags!
  char flags;
  data_stream.read(&flags, 1);
  int chunk_size;
  int size_flag = flags & 0x03;
  if (size_flag == 0) chunk_size = 1;
  else if (size_flag == 1) chunk_size = 2;
  else if (size_flag == 2) chunk_size = 4;
  else chunk_size = 8;

  bool attr_stored = (flags>>4)&1;
  bool times_stored = (flags>>5)&1;

  if(times_stored){
    data_stream.ignore(4*4);
  }

  if(attr_stored){
    data_stream.ignore(2*2);
  }

  //data_stream.ignore(chunk_size);
  unsigned long long actual_chunk_size = 0;
  data_stream.read(reinterpret_cast<char *>(&actual_chunk_size), chunk_size);

  long long bytes_remaining = static_cast<long long>(actual_chunk_size) - 4;


  //have link inof message, group info message, attribute info message, and nill message
  // Go through messages
  root_grid_size_found = false;
  dataset_names_found = false;
  variable_names_found = false;
  num_variables_found = false;
  tree_n_heap_address_found = false;
  while(!(root_grid_size_found && dataset_names_found && variable_names_found && num_variables_found 
    && tree_n_heap_address_found)&&bytes_remaining > 0){
    
    // Read message type and size
    unsigned char raw_type;
    unsigned short int message_size;
    data_stream.read(reinterpret_cast<char *>(&raw_type), 1);
    unsigned short int message_type = raw_type; // Safe upcast

    std::printf("message: %d ",message_type);
    data_stream.read(reinterpret_cast<char *>(&message_size), 2);

    // Check message flags
    unsigned char message_flags;
    data_stream.read(reinterpret_cast<char *>(&message_flags), 1);
    bytes_remaining -= 4;

    //deleting this check and hoping it's fine because I don't think this is strictly bad
    if (message_flags & 0b00000010)
      throw BlacklightException("Unexpected HDF5 header message flag.");
      

    // Read message data
    unsigned char *message_data = new unsigned char[message_size];
    data_stream.read(reinterpret_cast<char *>(message_data), message_size);

    // Follow any continuation messages
    if (message_type == 16)
    {
      unsigned long int new_offset;
      std::memcpy(&new_offset, message_data, 8);
      data_stream.seekg(static_cast<std::streamoff>(new_offset));
    }

    // Inspect any attribute messages
    else if (message_type == 12)
    {

      // Check attribute message version
      int offset = 0;
      if (message_data[offset] != 1)
        throw BlacklightException("Unexpected HDF5 attribute message version.");
      offset += 2;

      // Read attribute message metadata
      unsigned short int name_size, datatype_size, dataspace_size;
      std::memcpy(&name_size, message_data + offset, 2);
      offset += 2;
      std::memcpy(&datatype_size, message_data + offset, 2);
      offset += 2;
      std::memcpy(&dataspace_size, message_data + offset, 2);
      offset += 2;
      unsigned short int name_size_pad = static_cast<unsigned short int>((8 - name_size % 8) % 8);
      unsigned short int datatype_size_pad =
          static_cast<unsigned short int>((8 - datatype_size % 8) % 8);
      unsigned short int dataspace_size_pad =
          static_cast<unsigned short int>((8 - dataspace_size % 8) % 8);

      // Read attribute message data
      unsigned char *name_raw = new unsigned char[name_size];
      unsigned char *datatype_raw = new unsigned char[datatype_size];
      unsigned char *dataspace_raw = new unsigned char[dataspace_size];
      std::memcpy(name_raw, message_data + offset, name_size);
      offset += name_size + name_size_pad;
      std::memcpy(datatype_raw, message_data + offset, datatype_size);
      offset += datatype_size + datatype_size_pad;
      std::memcpy(dataspace_raw, message_data + offset, dataspace_size);
      offset += dataspace_size + dataspace_size_pad;
      std::string name(reinterpret_cast<char *>(name_raw),
          static_cast<std::string::size_type>(static_cast<int>(name_size) - 1));

      // Read and set desired attributes
      if (name == "RootGridSize")
      {
        root_grid_size_found = true;
        Array<int> root_grid_size;
        SetHDF5IntArray(datatype_raw, dataspace_raw, message_data + offset, root_grid_size);
        n_3_root = root_grid_size(2);
      }
      else if (name == "DatasetNames")
      {
        dataset_names_found = true;
        SetHDF5StringArray(datatype_raw, dataspace_raw, message_data + offset,
            first_time_root_object_header, &dataset_names, &num_dataset_names);
      }
      else if (name == "VariableNames")
      {
        variable_names_found = true;
        SetHDF5StringArray(datatype_raw, dataspace_raw, message_data + offset,
            first_time_root_object_header, &variable_names, &num_variable_names);

      }
      else if (name == "NumVariables")
      {
        num_variables_found = true;
        SetHDF5IntArray(datatype_raw, dataspace_raw, message_data + offset, num_variables);
      }

      // Free raw buffers
      delete[] name_raw;
      delete[] datatype_raw;
      delete[] dataspace_raw;
    }

    else if(message_type == 17){
      std::printf("Found Symbol Table Message!\n");
      // Safely extract from the data buffer we already pulled
      std::memcpy(&root_btree_address, message_data, 8);
      std::memcpy(&root_name_heap_address, message_data + 8, 8);
      tree_n_heap_address_found = true;
    }

    else if(message_type == 2){
     std::printf("link info message! \n");

      // Byte 0 is version. Byte 1 is the flags byte.
      unsigned char link_flags = message_data[1];
      
      // Bit 0 of the flags determines if the 8-byte Creation Index is present
      bool creation_index_present = (link_flags & 0x01);

      int offset = 2; // Start right after Version and Flags
      if (creation_index_present) {
        offset += 8; // Skip the 8-byte Maximum Creation Index
      }

      // Read Fractal Heap Address (8 bytes)
      std::memcpy(&root_name_heap_address, message_data + offset, 8);
      offset += 8;

      // Read v2 B-tree Address (8 bytes)
      std::memcpy(&root_btree_address, message_data + offset, 8);

      tree_n_heap_address_found = true;
      
      // Now that the pointer is fully accurate, traverse the B-tree
      ReadHDF5AttributeBTreeVers2(root_btree_address);
    }

    // Free raw buffer
    delete[] message_data;
    
    bytes_remaining -= message_size;
    
  }
  // Skip the remaining bytes (like the 4-byte chunk checksum) to leave the stream clean
  if (bytes_remaining >= 0) {
      data_stream.ignore(bytes_remaining + 4); 
  }

  // Check that appropriate messages were found
  if (not (root_grid_size_found and dataset_names_found and variable_names_found
      and num_variables_found and tree_n_heap_address_found)){
        std::printf("root_grid_size_found:%d  dataset_names_found: %d  variable_names_found: %d num_variables_found: %d tree_n_head_address_found: %d \n ",root_grid_size_found,dataset_names_found,variable_names_found,num_variables_found,tree_n_heap_address_found);
        throw BlacklightException("Could not find needed file-level attributes.");
      }
  if (num_variables.n1 != num_dataset_names)
    throw BlacklightException("DatasetNames and NumVariables file-level attribute mismatch.");

  // Update first time flag
  first_time_root_object_header = false;
  return;
}

// Function to traverse an HDF5 Version 2 B-tree looking for attributes
// Inputs: 
//   btree_address: Global file offset to the B-tree header block
// Outputs: (none)
void MCReader::ReadHDF5AttributeBTreeVers2(unsigned long int btree_address)
{
  if (btree_address == 0 || btree_address == 0xFFFFFFFFFFFFFFFFULL) return;

  data_stream.seekg(static_cast<std::streamoff>(btree_address));

  // Verify B-tree v2 Signature: 'B','T','H','D'
  const unsigned char expected_sig[4] = {'B', 'T', 'H', 'D'};
  //when printed out I got -1 in the tree signature
  for (int i = 0; i < 4; i++) {
    //std::printf(" %d ",data_stream.get());
    if (data_stream.get() != expected_sig[i])
      throw BlacklightException("Unexpected HDF5 v2 B-tree signature.");
  }

  // Version must be 0
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 v2 B-tree version.");

  // Type 5 is typically used for Attribute Name Indexes
  unsigned char btree_type = data_stream.get();
  //data_stream.ignore(2); // Node size (2 bytes), Number of records (2 bytes) processed below

  unsigned short int node_size;
  unsigned short int num_records;
  data_stream.read(reinterpret_cast<char*>(&node_size), 4);
  data_stream.read(reinterpret_cast<char*>(&num_records), 2);
  
  // 6. Depth (1 byte)
  unsigned char depth = data_stream.get();

  // 7. Split & Merge Percentages + Reserved (5 bytes total)
  data_stream.ignore(5);

  // 8. Number of records (2 bytes)
  //unsigned short int num_records;
  data_stream.read(reinterpret_cast<char*>(&num_records), 2);

  // Track state of found variables globally or via class properties
  bool root_grid_size_found = false;
  bool dataset_names_found = false;
  bool variable_names_found = false;
  bool num_variables_found = false;

  // Each record in an Attribute Name Index contains a Name Hash (4 bytes) 
  // and a Fractal Heap ID (size depends on file configuration, usually 6 or 8 bytes)
  // For most v2 configurations, the Heap ID size is 8 bytes.
  int heap_id_size = 8; 

  for (int n = 0; n < num_records; n++)
  {
    unsigned int name_hash;
    data_stream.read(reinterpret_cast<char*>(&name_hash), 4);

    unsigned char* heap_id = new unsigned char[heap_id_size];
    data_stream.read(reinterpret_cast<char*>(heap_id), heap_id_size);

    // Save current stream location before jumping into the Fractal Heap
    std::streampos return_pos = data_stream.tellg();

    // Pass the Heap ID to the reader to resolve the attribute data
    std::string attr_name = "";
    unsigned char* datatype_raw = nullptr;
    unsigned char* dataspace_raw = nullptr;
    unsigned char* data_raw = nullptr;

    // Call the secondary helper to lookup this object within the Fractal Heap
    ReadHDF5FractalHeapObject(root_name_heap_address, heap_id, heap_id_size, 
                              attr_name, &datatype_raw, &dataspace_raw, &data_raw);

    // Reconstruct your specific logic mapping matching metadata
    if (!attr_name.empty())
    {
      if (attr_name == "RootGridSize")
      {
        root_grid_size_found = true;
        Array<int> root_grid_size;
        SetHDF5IntArray(datatype_raw, dataspace_raw, data_raw, root_grid_size);
        n_3_root = root_grid_size(2);
      }
      else if (attr_name == "DatasetNames")
      {
        dataset_names_found = true;
        SetHDF5StringArray(datatype_raw, dataspace_raw, data_raw,
            first_time_root_object_header, &dataset_names, &num_dataset_names);
      }
      else if (attr_name == "VariableNames")
      {
        variable_names_found = true;
        SetHDF5StringArray(datatype_raw, dataspace_raw, data_raw,
            first_time_root_object_header, &dataset_names, &num_variable_names);
      }
      else if (attr_name == "NumVariables")
      {
        num_variables_found = true;
        SetHDF5IntArray(datatype_raw, dataspace_raw, data_raw, num_variables);
      }
      
      //delete[] attr_data;
    }
    
    delete[] datatype_raw;
    delete[] dataspace_raw;
    delete[] data_raw;
    delete[] heap_id;

    // Restore stream to parse the next B-tree record row
    data_stream.seekg(return_pos);
  }
}

// Function to resolve a dynamic object out of the HDF5 Fractal Heap
// Inputs:
//   heap_address: Global address targeting the fractal heap root block
//   heap_id: Raw buffer holding the unique ID extracted from the record
//   heap_id_len: The configuration length byte sizing of the active ID
// Outputs:
//   out_name: String variable populated with the target attribute descriptor 
//   out_data: Binary array returned holding payload values
void MCReader::ReadHDF5FractalHeapObject(unsigned long int heap_address, 
                                        unsigned char* heap_id, 
                                        int heap_id_len,
                                        std::string &out_name, 
                                        unsigned char** out_datatype_raw,
                                        unsigned char** out_dataspace_raw,
                                        unsigned char** out_data_raw)
{
  // Initialize pointers to safe null state
  *out_datatype_raw = nullptr;
  *out_dataspace_raw = nullptr;
  *out_data_raw = nullptr;

  data_stream.seekg(static_cast<std::streamoff>(heap_address));

  unsigned char expected_heap_sig[4] = {'F', 'R', 'H', 'P'};
  for (int i = 0; i < 4; i++) {
    if (data_stream.get() != expected_heap_sig[i])
      throw BlacklightException("Unexpected HDF5 Fractal Heap signature.");
  }
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 Fractal Heap version.");

  data_stream.ignore(2 + 2 + 1); 
  
  unsigned long int data_block_address;
  data_stream.read(reinterpret_cast<char*>(&data_block_address), 8);
  if (data_block_address == 0 || data_block_address == 0xFFFFFFFFFFFFFFFFULL) return;

  unsigned long int internal_offset = 0;
  std::memcpy(&internal_offset, heap_id + 1, (heap_id_len - 1 < 8) ? heap_id_len - 1 : 8);

  // Jump to actual attribute message storage position
  data_stream.seekg(static_cast<std::streamoff>(data_block_address + internal_offset));

  // Check Version Flag (Must be 3 for modern v2 object tracking configurations)
  unsigned char attribute_version = data_stream.get();
  if (attribute_version != 3) return;

  // Flags field (1 byte)
  data_stream.ignore(1);

  // Read message segments layout sizes
  unsigned short int name_size, datatype_size, dataspace_size;
  data_stream.read(reinterpret_cast<char*>(&name_size), 2);
  data_stream.read(reinterpret_cast<char*>(&datatype_size), 2);
  data_stream.read(reinterpret_cast<char*>(&dataspace_size), 2);

  // Character set byte (1 byte)
  data_stream.ignore(1);

  // 1. Read the Attribute Name String
  char* raw_string = new char[name_size + 1];
  data_stream.read(raw_string, name_size);
  raw_string[name_size] = '\0';
  out_name = std::string(raw_string);
  delete[] raw_string;

  // 2. Read raw Datatype structure
  *out_datatype_raw = new unsigned char[datatype_size];
  data_stream.read(reinterpret_cast<char*>(*out_datatype_raw), datatype_size);

  // 3. Read raw Dataspace structure
  *out_dataspace_raw = new unsigned char[dataspace_size];
  data_stream.read(reinterpret_cast<char*>(*out_dataspace_raw), dataspace_size);

  // Determine payload data size dynamically. The actual variable attribute data values 
  // occupy whatever space is left over at the end of the message.
  // Note: For certain complex types, HDF5 encodes a data size here, but since it is 
  // standard raw data, we can calculate its size based on the specific datatypes.
  // For safety with your integer/string layouts, calculate the value payload space needed:
  
  unsigned int data_size = 0;
  if (out_name == "RootGridSize") {
    data_size = 3 * 4; // 3 element dimension array of 4-byte integers
  } else {
    // For variable strings or arrays, extract data size descriptor directly out of datatype structure metadata
    std::memcpy(&data_size, (*out_datatype_raw) + 4, 4);
    
    // For array dimensions, we check if multiple elements are packed together:
    // DatasetNames and VariableNames string arrays are sized to num_dataset_names / num_variable_names
    if (out_name == "DatasetNames") data_size *= num_dataset_names;
    if (out_name == "VariableNames") data_size *= num_variable_names;
  }

  // 4. Read the Raw Data Payload Block
  if (data_size > 0) {
    *out_data_raw = new unsigned char[data_size];
    data_stream.read(reinterpret_cast<char*>(*out_data_raw), data_size);
  }
}
//--------------------------------------------------------------------------------------------------

// Function to read a scalar, single-precision file attribute
// Inputs:
//   attribute_name: name of attribute
// Outputs:
//   *p_val: value set
// Notes:
//   Assumes root_object_header_address set.
//   Changes stream pointer.
//   Must have object header version 1.
//   Must not have shared header messages.
//   Must have attribute message version 1.
//   Must have size of offsets 8.
//   Must be run on little-endian machine.
void MCReader::ReadHDF5FloatAttribute(const char *attribute_name, float *p_val)
{
  // Check object header version
  data_stream.seekg(static_cast<std::streamoff>(root_object_header_address));
  if (data_stream.get() != 1)
    throw BlacklightException("Unexpected HDF5 object header version.");
  data_stream.ignore(1);

  // Read number of header messages
  unsigned short int num_messages;
  data_stream.read(reinterpret_cast<char *>(&num_messages), 2);

  // Skip reading object reference count and object header size
  data_stream.ignore(8);

  // Align to 8 bytes within header (location of padding not documented)
  data_stream.ignore(4);

  // Go through messages
  bool attribute_found = false;
  for (int n = 0; n < num_messages; n++)
  {
    // Read message type and size
    unsigned short int message_type, message_size;
    data_stream.read(reinterpret_cast<char *>(&message_type), 2);
    data_stream.read(reinterpret_cast<char *>(&message_size), 2);

    // Check message flags
    unsigned char message_flags;
    data_stream.read(reinterpret_cast<char *>(&message_flags), 1);
    data_stream.ignore(3);
    if (message_flags & 0b00000010)
      throw BlacklightException("Unexpected HDF5 header message flag.");

    // Read message data
    unsigned char *message_data = new unsigned char[message_size];
    data_stream.read(reinterpret_cast<char *>(message_data), message_size);

    // Follow any continuation messages
    if (message_type == 16)
    {
      unsigned long int new_offset;
      std::memcpy(&new_offset, message_data, 8);
      data_stream.seekg(static_cast<std::streamoff>(new_offset));
      continue;
    }

    // Inspect any attribute messages
    else if (message_type == 12)
    {

      // Check attribute message version
      int offset = 0;
      if (message_data[offset] != 1)
        throw BlacklightException("Unexpected HDF5 attribute message version.");
      offset += 2;

      // Read attribute message metadata
      unsigned short int name_size, datatype_size, dataspace_size;
      std::memcpy(&name_size, message_data + offset, 2);
      offset += 2;
      std::memcpy(&datatype_size, message_data + offset, 2);
      offset += 2;
      std::memcpy(&dataspace_size, message_data + offset, 2);
      offset += 2;
      unsigned short int name_size_pad = static_cast<unsigned short int>((8 - name_size % 8) % 8);
      unsigned short int datatype_size_pad =
          static_cast<unsigned short int>((8 - datatype_size % 8) % 8);
      unsigned short int dataspace_size_pad =
          static_cast<unsigned short int>((8 - dataspace_size % 8) % 8);

      // Read attribute message data
      unsigned char *name_raw = new unsigned char[name_size];
      unsigned char *datatype_raw = new unsigned char[datatype_size];
      unsigned char *dataspace_raw = new unsigned char[dataspace_size];
      std::memcpy(name_raw, message_data + offset, name_size);
      offset += name_size + name_size_pad;
      std::memcpy(datatype_raw, message_data + offset, datatype_size);
      offset += datatype_size + datatype_size_pad;
      std::memcpy(dataspace_raw, message_data + offset, dataspace_size);
      offset += dataspace_size + dataspace_size_pad;
      std::string name(reinterpret_cast<char *>(name_raw),
          static_cast<std::string::size_type>(static_cast<int>(name_size) - 1));

      // Read and set desired attributes
      if (name == attribute_name)
      {
        std::printf("Found attribute %s\n",attribute_name);
        attribute_found = true;
        Array<float> attribute;
        SetHDF5FloatArray(datatype_raw, dataspace_raw, message_data + offset, attribute);
        *p_val = attribute(0);
      }

      // Free raw buffers
      delete[] name_raw;
      delete[] datatype_raw;
      delete[] dataspace_raw;
    }

    // Free raw buffer
    delete[] message_data;

    // Break when required information found
    if (attribute_found)
      break;
  }

  // Check that appropriate message was found
  if (not attribute_found)
    throw BlacklightException("Could not find needed file-level attributes.");
  return;
}