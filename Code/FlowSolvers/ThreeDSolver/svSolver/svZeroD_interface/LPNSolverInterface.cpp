#include "LPNSolverInterface.h"
#include <dlfcn.h>
#include <iostream>
#include <string>

//--------------------
// LPNSolverInterface
//--------------------
//
LPNSolverInterface::LPNSolverInterface()
{
//// Set the default names of the LPN interface functions.
  lpn_initialize_name_ = "initialize";
  lpn_increment_time_name_ = "increment_time";
  lpn_run_simulation_name_ = "run_simulation";
  lpn_update_block_params_name_ = "update_block_params";
  lpn_read_block_params_name_ = "read_block_params";
  lpn_get_block_node_IDs_name_ = "get_block_node_IDs";
  lpn_update_state_name_ = "update_state";
  lpn_return_ydot_name_ = "return_ydot";
  lpn_set_external_step_size_name_ = "set_external_step_size";
}

LPNSolverInterface::~LPNSolverInterface()
{
  dlclose(library_handle_);
}

//--------------
// load_library
//--------------
// Load the LPN shared library and get pointers to its interface functions.
//
void LPNSolverInterface::load_library(const std::string& interface_lib)
{
  library_handle_ = dlopen(interface_lib.c_str(), RTLD_LAZY);

  if (!library_handle_) {
    std::cerr << "Error loading shared library '" << interface_lib << "'  with error: " << dlerror() << std::endl;
    return;
  }

  // Get a pointer to the svzero 'initialize' function.
  *(void**)(&lpn_initialize_) = dlsym(library_handle_, "initialize");
  if (!lpn_initialize_) {
    std::cerr << "Error loading function 'initialize' with error: " << dlerror() << std::endl;
    dlclose(library_handle_);
    return;
  }

  // Get a pointer to the svzero 'increment_time' function.
  *(void**)(&lpn_increment_time_) = dlsym(library_handle_, "increment_time");
  if (!lpn_increment_time_) {
    std::cerr << "Error loading function 'increment_time' with error: " << dlerror() << std::endl;
    dlclose(library_handle_);
    return;
  }
  
  // Get a pointer to the svzero 'run_simulation' function.
  *(void**)(&lpn_run_simulation_) = dlsym(library_handle_, "run_simulation");
  if (!lpn_run_simulation_) {
    std::cerr << "Error loading function 'run_simulation' with error: " << dlerror() << std::endl;
    dlclose(library_handle_);
    return;
  }

  // Get a pointer to the svzero 'update_block_params' function.
  *(void**)(&lpn_update_block_params_) = dlsym(library_handle_, "update_block_params");
  if (!lpn_update_block_params_) {
    std::cerr << "Error loading function 'update_block_params' with error: " << dlerror() << std::endl;
    dlclose(library_handle_);
    return;
  }

  // Get a pointer to the svzero 'read_block_params' function.
  *(void**)(&lpn_read_block_params_) = dlsym(library_handle_, "read_block_params");
  if (!lpn_read_block_params_) {
    std::cerr << "Error loading function 'read_block_params' with error: " << dlerror() << std::endl;
    dlclose(library_handle_);
    return;
  }
  
  // Get a pointer to the svzero 'get_block_node_IDs' function.
  *(void**)(&lpn_get_block_node_IDs_) = dlsym(library_handle_, "get_block_node_IDs");
  if (!lpn_get_block_node_IDs_) {
    std::cerr << "Error loading function 'lpn_get_block_node_IDs' with error: " << dlerror() << std::endl;
    dlclose(library_handle_);
    return;
  }
  
  // Get a pointer to the svzero 'update_state' function.
  *(void**)(&lpn_update_state_) = dlsym(library_handle_, "update_state");
  if (!lpn_update_state_) {
    std::cerr << "Error loading function 'lpn_update_state' with error: " << dlerror() << std::endl;
    dlclose(library_handle_);
    return;
  }
  
  // Get a pointer to the svzero 'return_ydot' function.
  *(void**)(&lpn_return_ydot_) = dlsym(library_handle_, "return_ydot");
  if (!lpn_return_ydot_) {
    std::cerr << "Error loading function 'lpn_return_ydot' with error: " << dlerror() << std::endl;
    dlclose(library_handle_);
    return;
  }
  
  // Get a pointer to the svzero 'set_external_step_size' function.
  *(void**)(&lpn_set_external_step_size_) = dlsym(library_handle_, "set_external_step_size");
  if (!lpn_set_external_step_size_) {
    std::cerr << "Error loading function 'lpn_set_external_step_size' with error: " << dlerror() << std::endl;
    dlclose(library_handle_);
    return;
  }
}

//------------
// initialize
//------------
// Initialze the LPN solver.
//
// Parameters:
//
//   file_name: The name of the LPN configuration file (JSON).
//
//   time_step: The time step to use in the LPN simulation. 
//
void LPNSolverInterface::initialize(std::string file_name)
{
  lpn_initialize_(file_name, problem_id_, pts_per_cycle_, num_cycles_, num_output_steps_, block_names_, variable_names_);
  std::cout << "[LPNSolverInterface::initialize] Problem ID: " << problem_id_ << std::endl;
  system_size_ = variable_names_.size();
  std::cout << "[LPNSolverInterface::initialize] System size: " << system_size_ << std::endl;
  //solution_.resize(system_size_*num_output_steps_);
}

//----------------
// set_external_step_size
//----------------
void LPNSolverInterface::set_external_step_size(double step_size)
{
  lpn_set_external_step_size_(problem_id_, step_size);
}

//----------------
// increment_time
//----------------
// Increment the LPN solution in time.
//
// Parameters:
//
//   time: The solution time.
//
//   solution: The returned LPN solution.
//
void LPNSolverInterface::increment_time(const double time, std::vector<double>& solution)
{
  lpn_increment_time_(problem_id_, time, solution);
}

//----------------
// run_simulation
//----------------
void LPNSolverInterface::run_simulation(const double time, std::vector<double>& output_times, std::vector<double>& output_solutions, int& error_code)
{
  lpn_run_simulation_(problem_id_, time, output_times, output_solutions, error_code);
}

//----------------
// update_block_params
//----------------
void LPNSolverInterface::update_block_params(std::string block_name, std::vector<double>& new_params)
{
  lpn_update_block_params_(problem_id_, block_name, new_params);
}

//----------------
// read_block_params
//----------------
void LPNSolverInterface::read_block_params(std::string block_name, std::vector<double>& block_params)
{
  lpn_read_block_params_(problem_id_, block_name, block_params);
}

//----------------
// get_variable_ids
//----------------
void LPNSolverInterface::get_block_node_IDs(std::string block_name, std::vector<int>& IDs)
{
  lpn_get_block_node_IDs_(problem_id_, block_name, IDs);
}

//----------------
// update_state
//----------------
void LPNSolverInterface::update_state(std::vector<double> state_y, std::vector<double> state_ydot)
{
  lpn_update_state_(problem_id_, state_y, state_ydot);
}

//----------------
// get_variable_ids
//----------------
void LPNSolverInterface::return_ydot(std::vector<double>& ydot)
{
  lpn_return_ydot_(problem_id_, ydot);
}
