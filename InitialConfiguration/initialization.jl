function initialize_parameters(parameters)
    parameters["num_polymers_per_type"] = parameters["initial_values"]
    parameters = read_input_data(parameters) 
    parameters = basic_setup(parameters) 
    parameters = initialization_state(parameters)  
    return parameters
end