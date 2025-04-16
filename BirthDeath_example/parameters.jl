
# ### Parameters for Simulation
parameters = Dict(
    # Gillespie Algorithm Parameters
    "initial_values" => [100],                      # Gillespie initial values
    "reaction_matrix" => [1; -1],                   # Gillespie reaction
    "reaction_rates" => [0.1, 0.001],                # Gillespie generation rate, degradation rate (reaction_rates)
    
    "time_scale" => 0.1,                             # tau_D/tau_p
    "alpha" => 1,                                 # regulates the strength of the non-equilibrium correction

    # Initial Polymer Counts
    "droplet_polymers_count" => [30],               # Droplet polymers count
    "slab_polymers_count" => [200],                  # Slab polymers count
    "num_polymers_per_type" => [200],
    
    # Monte Carlo Move
    "single_bead_move_count" => 100,
    "paired_bead_move_count" => 100,
    "reconnect_move_count" => 200,
    "reptation_move_count" => 100,
    "chain_translation_move_count" => 0,
    "chain_formed_cluster_move_count" => 100,
    "non_largest_cluster_move_count" => 1,
    "double_pivot_move_count" => 0
    )


