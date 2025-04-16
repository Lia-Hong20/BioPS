function run_simulation(parameters, STEPS_PER_SIMULATION, TOTAL_SIMULATION_STEPS)
    prod_count, deg_count, time = 0, 0, 0.0
    largest_cluster = 0  # 初始化 largest_cluster

    state = parameters["num_polymers_per_type"]
    for step in 1:STEPS_PER_SIMULATION:TOTAL_SIMULATION_STEPS
        println("Step:", step)

        move_step = 0
        for substep in step:(step + STEPS_PER_SIMULATION - 1)
            # biochemical reaction process (Gillespie)
            time_interval, reaction_index, updated_state = gillespie(state, propensity_function, parameters["reaction_matrix"], parameters["reaction_rates"]) #判断发生的反应
            state = updated_state
            time += time_interval
            old_polymer_count, old_energy = parameters["num_polymers_per_type"][1], parameters["total_energy"]
            
            # biochemical reaction mmc
            if reaction_index == 1 # 生成
                if parameters["total_polymers_count"] < MAX_TOTAL_POLYMERS[1] # 蛋白质数没有超过总数
                    parameters = protein_production(parameters)
                    prod_count += 1
                end
            elseif reaction_index == 2 # 降解
                if parameters["total_polymers_count"] > 0 #系统中还在有蛋白质
                    parameters = protein_degradation(parameters)
                    deg_count += 1
                end
            end
            new_polymer_count, new_energy = parameters["num_polymers_per_type"][1], parameters["total_energy"]
            
            # 修正
            delta_chemical_potential = kBT * log(new_polymer_count / old_polymer_count) + (new_energy - old_energy)
            time_scale_factor = parameters["alpha"] * (parameters["time_scale"]/(parameters["time_scale"] + 1))
            parameters["adjustment_factor"] = exp(-time_scale_factor * delta_chemical_potential / kBT)

            # 在delta_t时间间隔内执行相分离的MMC步骤
            if parameters["num_polymers_per_type"][1] > 0
                move_step = round(Int, time_interval / 0.001)
                for _ ∈ 1:move_step
                    parameters = move_process_simulation(parameters)
                end
            end

            scatter(parameters["bead_positions_list"][:, 1], parameters["bead_positions_list"][:, 2], parameters["bead_positions_list"][:, 3], xlims=[0, 100], ylims=[0, 100], zlims=[0, 100])

        end
        largest_cluster = find_largest_cluster(parameters)[1]
        average_Rg = center_of_mass(parameters)
        write_simulation_results(step, STEPS_PER_SIMULATION, parameters, largest_cluster, time, move_step)   
    end

    return largest_cluster, parameters
end