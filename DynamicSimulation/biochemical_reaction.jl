using Distributions
using DataFrames
using StaticArrays
using Plots
using DelimitedFiles
using Setfield

# 用累积和优化 pfsample 函数
function pfsample(w::AbstractArray{Float64,1}, s::Float64, n::Int64)
    t = rand() * s
    i = 1
    cw = w[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += w[i]
    end
    return i
end

# Gillespie 算法实现（避免覆盖原始 state）
function gillespie(state::AbstractVector{Int64}, propensity_function::Base.Callable, stoichiometry::AbstractVector{Int64}, params::AbstractVector{Float64})
    # 当前状态的倾向函数
    propensities = propensity_function(state, params)
    
    # 所有反应的总发生率
    total_propensity = sum(propensities)
    if total_propensity == 0.0
        error("All propensities are zero. No reaction can occur.")
    end
    
    # 计算时间间隔 dt
    dt = rand(Exponential(1 / total_propensity))
    
    # 选择发生的反应事件
    event_idx = pfsample(propensities, total_propensity, size(stoichiometry, 1))
    
    # 返回一个新的 state，而不是修改原始 state
    new_state = state .+ stoichiometry[event_idx, :]

    return dt, event_idx, new_state
end

# propensity function F
function propensity_function(state, params)
	(λ, δ) = params
	return [λ, state[1] * δ]
end

function protein_degradation(parameters)
	bead_positions_list = parameters["bead_positions_list"]
	polymers_type = parameters["polymers_type"]
	max_beads = parameters["max_beads"]
	interaction_list = parameters["interaction_list"]
	bead_index_lattice = parameters["bead_index_lattice"]
	neighbor_matrix = parameters["neighbor_matrix"]
	energy = parameters["energy"]
	total_energy = parameters["total_energy"]
	num_polymers_per_type = parameters["num_polymers_per_type"]

	# 删除
	current_polymer_index = rand(1:sum(num_polymers_per_type)) # 随机选择一条链
	current_polymer_type = polymers_type[current_polymer_index]

	# 删除与当前链有相互作用的珠子的绑定，并建立新的连接（注意删除后，删除的链后的珠子索引要前移）
	bead_range_start = (current_polymer_index - 1) * max_beads + 1
	bead_range_end = bead_range_start + BEADS_PER_POLYMER[current_polymer_type] - 1

	if any(interaction_list[bead_range_start:bead_range_end] .!= 0)  # 如果链上有绑定，删除
		for bead_index in bead_range_start:bead_range_end # 遍历当前链上的所有珠子
			interaction_bead_global_index = interaction_list[bead_index]

			if interaction_bead_global_index != 0
				current_bead_index = bead_global_to_local_index(bead_index, parameters, polymers_type)
				interaction_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type)
				
				# 更新能量：断开链上所有的相互作用，并且不形成新的链接
				total_energy = calculate_energy(total_energy, .-energy, current_bead_index, interaction_bead_index, parameters, bead_index, interaction_bead_global_index)
				interaction_list[bead_index] = 0 
				interaction_list[interaction_bead_global_index] = 0 # 不形成新的连接

				#= 更新能量：断开链上所有的相互作用，并且不形成新的链接
				if bead_index_to_polymer_index(bead_index, parameters) != bead_index_to_polymer_index(interaction_bead_global_index, parameters)
					total_energy = calculate_energy(total_energy, .-energy, current_bead_index, interaction_bead_index, parameters, bead_index, interaction_bead_global_index)
				else # 绑定的珠子在删除的链上
					if interaction_bead_global_index > bead_index 
						total_energy = calculate_energy(total_energy, .-energy, current_bead_index, interaction_bead_index, parameters, bead_index, interaction_bead_global_index)
					end
				end
                
				if bead_index_to_polymer_index(bead_index, parameters) != bead_index_to_polymer_index(interaction_bead_global_index, parameters) # 绑定的珠子不在删除的链上，被绑定的珠子换绑
					interaction_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type)
					interaction_candidates_count = 1 # 可相互作用的珠子数
					interaction_candidates = zeros(Int, 27)
					for neighbor_index in 1:26 #遍历周围格点
						neighbor_position = round.(Int, periodic_boundary_3d(bead_positions_list[interaction_bead_global_index, :] .+ neighbor_matrix[neighbor_index, :])) # 邻居格点的坐标
						neighbor_bead_index = bead_index_lattice[neighbor_position[1], neighbor_position[2], neighbor_position[3]]

						if neighbor_bead_index != 0 # 邻居格点处有珠子
							if energy[interaction_bead_index, bead_global_to_local_index(neighbor_bead_index, parameters, polymers_type)] != 0 && interaction_list[neighbor_bead_index] == 0 # 两个珠子可以相互作用，且周围的珠子没有绑定
								interaction_candidates_count += 1
								if neighbor_bead_index in bead_range_start:bead_range_end # 如果邻居格点处的珠子也在要删除的链上，不会绑定
									interaction_candidates[interaction_candidates_count] = 0
								else
									interaction_candidates[interaction_candidates_count] = neighbor_bead_index # 储存可以绑定的珠子索引
								end
							end
						end
					end

					if interaction_candidates_count > 1
						selected_interaction = interaction_candidates[rand(1:interaction_candidates_count)] # 可绑定的珠子列表中随机选择一个珠子
						if selected_interaction > 0
							interaction_list[interaction_bead_global_index] = selected_interaction # 更新连接
							interaction_list[selected_interaction] = interaction_bead_global_index
							current_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type)
							interaction_bead_index = bead_global_to_local_index(selected_interaction, parameters, polymers_type)
							total_energy = calculate_energy(total_energy, energy, current_bead_index, interaction_bead_index, parameters, interaction_bead_global_index, selected_interaction)
						else
							interaction_list[interaction_bead_global_index] = 0
						end
					else
						interaction_list[interaction_bead_global_index] = 0
					end
				end
				=#
			end
		end
	end
	# 更新变量
	interaction_list = interaction_list[[1:(bead_range_start-1); (bead_range_end+1):end], :] # 在相互作用列表interaction_list中删除该链
	logical_index = interaction_list .> bead_range_end
	interaction_list[logical_index] .-= BEADS_PER_POLYMER[current_polymer_type] # 后面的珠子索引前移

	for delete_bead ∈ bead_range_start:bead_range_end # 删除空间中的珠子
		bead_index_lattice[bead_positions_list[delete_bead, 1], bead_positions_list[delete_bead, 2], bead_positions_list[delete_bead, 3]] = 0
	end
	logical_index = bead_index_lattice .> bead_range_end
	bead_index_lattice[logical_index] .-= BEADS_PER_POLYMER[current_polymer_type]

	bead_positions_list = bead_positions_list[[1:(bead_range_start-1); (bead_range_end+1):end], :]

	deleteat!(polymers_type, current_polymer_index) # 删除该链对应的类型
	num_polymers_per_type[current_polymer_type] -= 1

	parameters["interaction_list"] = interaction_list
	parameters["bead_positions_list"] = bead_positions_list
	parameters["bead_index_lattice"] = bead_index_lattice
	parameters["num_polymers_per_type"] = num_polymers_per_type
	parameters["polymers_type"] = polymers_type
	parameters["total_energy"] = total_energy
	parameters["total_polymers_count"] = sum(parameters["num_polymers_per_type"])
	return parameters
end

function protein_production(parameters)
	bead_positions_list = parameters["bead_positions_list"]
	polymers_type = parameters["polymers_type"]
	max_beads = parameters["max_beads"]
	interaction_list = parameters["interaction_list"]
	bead_index_lattice = parameters["bead_index_lattice"]
	neighbor_matrix = parameters["neighbor_matrix"]
	energy = parameters["energy"]
	total_energy = parameters["total_energy"]
	num_polymers_per_type = parameters["num_polymers_per_type"]

	current_polymer_type = rand(1:NUM_POLYMER_TYPES) # 随机选择要添加的蛋白质的类型
	flag = false
	success = false  # 标志变量，记录是否成功更新
	add_bead_positions = []

	for _ in 1:100000
		bead_reference_position = rand_int_array(LATTICE_SIZE)
		bead_position_local = add_protein(current_polymer_type, parameters)
		flag, add_bead_positions = check_overlap(bead_position_local, bead_reference_position, bead_index_lattice, current_polymer_type)

		if flag # 放置成功，更新
			# 将加入的珠串信息添加到bead_positions_list, polymers_type, interaction_list后面, 该类型的链数＋1
			bead_positions_list = vcat(bead_positions_list, add_bead_positions)
			polymers_type = vcat(polymers_type, current_polymer_type)
			for bead_index in 1:BEADS_PER_POLYMER[current_polymer_type] # 新加入的链上的所有珠子放入空间bead_index_lattice
				current_bead_global_index = num_polymers_per_type[current_polymer_type] * max_beads + bead_index
				bead_index_lattice[add_bead_positions[bead_index, 1], add_bead_positions[bead_index, 2], add_bead_positions[bead_index, 3]] = current_bead_global_index
			end
			success = true  # 更新标志变量
			break
		end
	end
	if !success
		error("Consider adding more room or trying with a different seed")
	end

	interaction_list = vcat(interaction_list, zeros(Int, BEADS_PER_POLYMER[current_polymer_type]))
	# 新生成的蛋白质链与已有链的相互作用
	for bead_index in 1:BEADS_PER_POLYMER[current_polymer_type] # 遍历新加入的链上的所有珠子
		current_bead_global_index = num_polymers_per_type[current_polymer_type] * max_beads + bead_index
		if interaction_list[current_bead_global_index] == 0 # 当前珠子没有连接
			interaction_candidates = zeros(Int, 27)
			interaction_candidates_count = 1 # 可发生相互作用的珠子的数量
			for neighbor_index in 1:26
				neighbor_global_position = periodic_boundary_3d(add_bead_positions[bead_index, :] .+ neighbor_matrix[neighbor_index, :]) # 当前位置的周围坐标
				neighbor_bead_global_index = bead_index_lattice[neighbor_global_position[1], neighbor_global_position[2], neighbor_global_position[3]]
				if neighbor_bead_global_index != 0
					if energy[bead_index, bead_global_to_local_index(neighbor_bead_global_index, parameters, polymers_type)] != 0 && interaction_list[neighbor_bead_global_index] == 0 # 可相互作用，且与周围的珠子没有连接
						interaction_candidates_count += 1
						interaction_candidates[interaction_candidates_count] = neighbor_bead_global_index # 储存可以旋转的珠子索引
					end
				end
			end

			if interaction_candidates_count > 1
				interaction_bead_global_index = interaction_candidates[rand(1:interaction_candidates_count)] # 随机选择一个可以相互作用的珠子
				if interaction_bead_global_index > 0
					interaction_list[current_bead_global_index] = interaction_bead_global_index # 更新变量
					interaction_list[interaction_bead_global_index] = current_bead_global_index
					interaction_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type)
					total_energy = calculate_energy(total_energy, energy, bead_index, interaction_bead_index, parameters, current_bead_global_index, interaction_bead_global_index)
				end
			end
		end
	end
	num_polymers_per_type[current_polymer_type] += 1

	parameters["bead_positions_list"] = bead_positions_list
	parameters["polymers_type"] = polymers_type
	parameters["bead_index_lattice"] = bead_index_lattice
	parameters["interaction_list"] = interaction_list
	parameters["total_energy"] = total_energy
	parameters["num_polymers_per_type"] = num_polymers_per_type
	parameters["total_polymers_count"] = sum(num_polymers_per_type)
	return parameters
end
