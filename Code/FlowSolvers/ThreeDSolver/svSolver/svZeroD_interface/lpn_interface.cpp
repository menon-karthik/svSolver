
#include "LPNSolverInterface.h"

#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <string>
//#include <boost/tokenizer.hpp>

//static std::string lpn_library;
static std::map<int,LPNSolverInterface*> interfaces;

////--------------------
//// lpn_interface_init
////--------------------
////
//extern "C" void lpn_interface_init_(const char* lpn_library_name)
//{
//  lpn_library = std::string(lpn_library_name);

//}

//-------------------------
// lpn_interface_add_block
//-------------------------
//
extern "C" void lpn_interface_add_model_(const char* lpn_library_name, int* lib_filename_len, const char* lpn_json_file, int* json_filename_len, int* interface_id, int* num_output_steps, int* system_size)
{
  //std::cout<<std::string(lpn_library_name)<<std::endl;
  // Load library
  //lpn_library = std::string(lpn_library_name);
  std::string lpn_library(lpn_library_name,0,*lib_filename_len);
  std::cout<<"lpn_library: "<<lpn_library<<std::endl;
  //std::cout<<lpn_library.c_str()<<std::endl;
  //lpn_library = "/home/users/kmenon13/svZeroDPlus/Release/src/interface/libsvzero_interface_library.so";
  auto interface = new LPNSolverInterface();
  interface->load_library(lpn_library);
  std::cout << "[lpn_interface_init_] Loaded library." << std::endl;

  // Initialize model
//std::cout<<"lpn_json_file: "<<std::string(lpn_json_file)<<std::endl;
//std::string lpn_file(lpn_json_file,0,json_filename_len);
//std::cout<<"lpn_json_file (test): "<<test<<std::endl;
//lpn_json_file = "with_coronaries_ext.json";
//std::cout<<"lpn_json_file: "<<std::string(lpn_json_file)<<std::endl;
//interface->initialize(std::string(lpn_json_file));
  std::string lpn_file(lpn_json_file,0,*json_filename_len);
  interface->initialize(std::string(lpn_file));
  *interface_id = interface->problem_id_;
  std::cout << "[lpn_interface_add_model_] interface->problem_id_:" <<interface->problem_id_<< std::endl;
  interfaces[*interface_id] = interface;
  
  //std::cout << "[lpn_interface_add_model_] 1" << std::endl;
  
  // Save model parameters
  *num_output_steps = interface->num_output_steps_;
  //std::cout << "[lpn_interface_add_model_] 4" << std::endl;
  *system_size = interface->system_size_;

  //std::cout << "[lpn_interface_add_model_] END" << std::endl;
}

//----------------------------
// lpn_interface_set_external_step_size
//----------------------------
//
extern "C" void lpn_interface_set_external_step_size_(const int* interface_id, double* step_size)
{
  auto interface = interfaces[*interface_id];
  std::cout<<"[lpn_interface_set_external_step_size] step_size = "<<*step_size<<std::endl;
  interface->set_external_step_size(*step_size);
}

//----------------------------
// lpn_interface_get_variables_ids_
//----------------------------
//
extern "C" void lpn_interface_get_variable_ids_(const int* interface_id, char* block_name, int* block_name_len, int* blk_ids)
{
  auto interface = interfaces[*interface_id];
  std::vector<std::string> variable_names;
  variable_names = interface->variable_names_;
  for (int i = 0; i < variable_names.size(); i++) {
    std::cout<<"[lpn_interface_get_variable_ids] Variable name: "<<variable_names[i]<<" , index: "<< i<<std::endl;
  }
  std::vector<int> IDs;
  std::string block_name_cpp(block_name,0,*block_name_len);
  interface->get_block_node_IDs(std::string(block_name_cpp), IDs);
  // IDs stores info in the following format:
  // {num inlet nodes, inlet flow[0], inlet pressure[0],..., num outlet nodes, outlet flow[0], outlet pressure[0],...}
  std::cout<<"IDs: ";
  for (int i = 0;i<IDs.size();i++) {
    std::cout<<IDs[i]<<", ";
  }
  std::cout<<std::endl;
  int num_inlet_nodes = IDs[0];
  int num_outlet_nodes = IDs[1+num_inlet_nodes*2];
  std::cout<<"Inlet and outlet nodes: "<<num_inlet_nodes<<", "<<num_outlet_nodes<<std::endl;
  if ((num_inlet_nodes == 0) && (num_outlet_nodes = 1)) {
    std::cout<<"Only outlet nodes"<<std::endl;
    blk_ids[0] = IDs[1+num_inlet_nodes*2+1]; // Outlet flow
    blk_ids[1] = IDs[1+num_inlet_nodes*2+2]; // Outlet pressure
  } else if ((num_inlet_nodes == 1) && (num_outlet_nodes == 0)) {
    std::cout<<"Only inlet nodes"<<std::endl;
    blk_ids[0] = IDs[1]; // Inlet flow
    blk_ids[1] = IDs[2]; // Inlet pressure
  } else {
    std::runtime_error("ERROR: [lpn_interface_get_variable_ids] Not a flow/pressure block.");
  }
  std::cout<<"blk_ids[0]: "<<blk_ids[0]<<std::endl;
  std::cout<<"blk_ids[1]: "<<blk_ids[1]<<std::endl;

//std::cout<<"[lpn_interface_get_variable_ids] Block name: "<<block_name<<std::endl;
//for (int i = 0; i < num_inlet_nodes; i++) {
//  std::cout<<"[lpn_interface_get_variable_ids] Inlet node IDs: "<<IDs[1+2*i]<<", "<<IDs[1+2*i+1];
//}
//std::cout<<std::endl;
//for (int i = 0; i < num_outlet_nodes; i++) {
//  std::cout<<"[lpn_interface_get_variable_ids] Outlet node IDs: "<<IDs[2+num_inlet_nodes*2+2*i]<<", "<<IDs[2+num_inlet_nodes*2+2*i+1];
//}
//std::cout<<std::endl;
}

//// Functions to trim a string
//static inline std::string &ltrim(std::string &s) {
//    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int c) {return !std::isspace(c);}));
//    return s;
//}
//static inline std::string &rtrim(std::string &s) {
//    s.erase(std::find_if(s.rbegin(), s.rend(), [](int c) {return !std::isspace(c);}).base(), s.end());
//    return s;
//}

//----------------------------
// lpn_interface_get_solution
//----------------------------
//
extern "C" void lpn_interface_get_solution_(const int* interface_id, const double* time, double* solution)
{
  auto interface = interfaces[*interface_id];
  std::vector<double> lpn_solution(interface->system_size_);
  interface->increment_time(*time, lpn_solution);

  for (int i = 0; i < lpn_solution.size(); i++) {
    solution[i] = lpn_solution[i];
  }
}

//----------------------------
// lpn_interface_update_state
//----------------------------
//
extern "C" void lpn_interface_update_state_(const int* interface_id, double* y, double* ydot)
{
  auto interface = interfaces[*interface_id];
  std::vector<double> state_y(interface->system_size_);
  std::vector<double> state_ydot(interface->system_size_);
  for (int i = 0; i < interface->system_size_; i++) {
    state_y[i] = y[i];
    state_ydot[i] = ydot[i];
  }
  interface->update_state(state_y, state_ydot);
}

//----------------------------
// lpn_interface_return_ydot
//----------------------------
//
extern "C" void lpn_interface_return_ydot_(const int* interface_id, double* ydot)
{
  auto interface = interfaces[*interface_id];
  std::vector<double> state_ydot(interface->system_size_);
  interface->return_ydot(state_ydot);
  for (int i = 0; i < interface->system_size_; i++) {
    ydot[i] = state_ydot[i];
    //std::cout<<"[lpn_interface_return_ydot] "<<state_ydot[i]<<std::endl;
  }
}

//----------------------------
// lpn_interface_get_solution
//----------------------------
//
extern "C" void lpn_interface_update_block_params_(const int* interface_id, char* block_name, int* block_name_len, double* time, double* params, int* num_time_pts)
{
  std::cout << "lpn_interface_update_block_params START" << std::endl;
  auto interface = interfaces[*interface_id];
  //std::cout << "lpn_interface_update_block_params 1" << std::endl;
  int param_len = *num_time_pts;
  std::vector<double> new_params(1+2*param_len);
  //std::cout << new_params.size() << std::endl;
  //std::cout << param_len << std::endl;
  new_params[0] = (double) param_len;
  //double temp;
  //temp = time[0];
  //std::cout << "lpn_interface_update_block_params 2" << std::endl;
  //temp = time[1];
  //std::cout << "lpn_interface_update_block_params 2.1" << std::endl;
  for (int i = 0; i < param_len; i++) {
    //std::cout<<i << std::endl;
    //std::cout << "lpn_interface_update_block_params time["<<i<<"]: "<<time[i]<< std::endl;
    new_params[1+i] = time[i];
    new_params[1+param_len+i] = params[i];
  }
  //std::string block = std::string(block_name);
  //std::cout << "lpn_interface_update_block_params block_name length: "<<block.length()<<std::endl;
  //block = "INFLOW";
  //std::cout << "lpn_interface_update_block_params block_name length: "<<block.length()<<std::endl;
  //std::cout << "lpn_interface_update_block_params 3" << std::endl;
  std::string block_name_cpp(block_name,0,*block_name_len);
  interface->update_block_params(std::string(block_name_cpp), new_params);
  std::cout << "lpn_interface_update_block_params END" << std::endl;
}

//----------------------------
// lpn_interface_run_simulation
//----------------------------
//
extern "C" void lpn_interface_run_simulation_(const int* interface_id, const double* time, double* lpn_times, double* lpn_solutions, int* error_code)
{
  auto interface = interfaces[*interface_id];
  int solutions_vec_size = interface->system_size_*interface->num_output_steps_;
  std::vector<double> lpn_solutions_vec(solutions_vec_size);
  std::vector<double> lpn_times_vec(interface->num_output_steps_);
  int error_code_ret = 0;
  interface->run_simulation(*time, lpn_times_vec, lpn_solutions_vec, error_code_ret);
  *error_code = error_code_ret;

  for (int i = 0; i < interface->num_output_steps_; i++) {
    lpn_times[i] = lpn_times_vec[i];
  }
  for (int i = 0; i < solutions_vec_size; i++) {
    lpn_solutions[i] = lpn_solutions_vec[i];
  }
}
