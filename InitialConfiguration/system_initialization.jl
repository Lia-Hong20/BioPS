function initialization_state(parameters) #在系统中放置蛋白质链
	total_polymers_count = parameters["total_polymers_count"]
	energy = parameters["energy"]
	droplet_polymers_count = parameters["droplet_polymers_count"]
	slab_polymers_count = parameters["slab_polymers_count"]
	neighbor_matrix = parameters["neighbor_matrix"]
	max_beads = parameters["max_beads"]
	num_polymers_per_type = parameters["num_polymers_per_type"]
   
	bead_index_lattice = zeros(Int, LATTICE_SIZE[1], LATTICE_SIZE[2], LATTICE_SIZE[3]) #lattice[x, y, z]=k:在坐标(x, y, z)处存储了一个蛋白质链的珠子，该珠子在List中的索引值为k
	bead_positions_list = zeros(Int, max_beads * total_polymers_count, 3) #每行代表一个珠子的三维坐标
	polymers_type = zeros(Int, total_polymers_count)
	interaction_list = zeros(Int, max_beads * total_polymers_count)
	total_energy = 0.0

	# 放置蛋白质链
	for current_polymers_type in 1:NUM_POLYMER_TYPES
		for current_polymers_index in 1:num_polymers_per_type[current_polymers_type]
			place_protein(current_polymers_type, current_polymers_index, bead_index_lattice, bead_positions_list, polymers_type, max_beads, parameters)
		end
	end

	interaction_list, total_energy = form_interactions(NUM_POLYMER_TYPES, num_polymers_per_type, BEADS_PER_POLYMER, total_energy, energy, polymers_type, interaction_list, bead_index_lattice, bead_positions_list, neighbor_matrix, parameters)

	parameters["bead_index_lattice"] = bead_index_lattice
	parameters["bead_positions_list"] = bead_positions_list
	parameters["polymers_type"] = polymers_type
	parameters["interaction_list"] = interaction_list
	parameters["total_energy"] = total_energy
	parameters["propose_moves"] = zeros(Int, 8)
	parameters["rejection_moves"] = zeros(Int, 8)
	parameters["moves_counts"] = cumsum([
        parameters["single_bead_move_count"],
        parameters["paired_bead_move_count"],
        parameters["reconnect_move_count"],
        parameters["reptation_move_count"],
        parameters["chain_translation_move_count"],
        parameters["chain_formed_cluster_move_count"],
        parameters["non_largest_cluster_move_count"],
        parameters["double_pivot_move_count"],
    ])

	return parameters
end

function place_protein(current_polymers_type, current_polymers_index, bead_index_lattice, bead_positions_list, polymers_type, max_beads, parameters)
	for _ in 1:10000
        bead_reference_position = rand_int_array(LATTICE_SIZE) # 第一个珠子的位置
        if START_FROM_SLAB
            if current_polymers_index <= parameters["slab_polymers_count"][current_polymers_type]
                bead_reference_position[3] = rand(SLAB_HEIGHT) + (LATTICE_SIZE[3] - SLAB_HEIGHT) / 2
            end
        else
            if current_polymers_index <= parameters["droplet_polymers_count"][current_polymers_type]
                bead_reference_position = generate_droplet_position()
            end
        end

        bead_position_local = add_protein(current_polymers_type, parameters)
        flag, bead_absolute_position = check_overlap(bead_position_local, bead_reference_position, bead_index_lattice, current_polymers_type)

		if flag
			add_to_lattice(bead_absolute_position, current_polymers_type, current_polymers_index, bead_index_lattice, bead_positions_list, polymers_type, max_beads)
			return
		end
	end
	error("Consider adding more room or trying with a different seed")
end

function add_protein(polymer_type, parameters)
    total_connections = parameters["total_connections"]
    connection_matrix = parameters["connection_matrix"]
    connection_lengths = parameters["connection_lengths"]
    max_beads = parameters["max_beads"]
    
    protein_positions = zeros(Int, max_beads, 3) # 行：链所含的最大珠子数，列：坐标
    protein_positions[1, :] .= [0, 0, 0] # 第一个珠子的位置为原点(0, 0, 0)
    
    for bead_index in 2:BEADS_PER_POLYMER[polymer_type] # 对蛋白质链上的所有珠子遍历，除第一个珠子
        for _ in 1:100   # 对每个珠子，进行最多100次的尝试（防止无限循环）
            bead_global_index = sum(BEADS_PER_POLYMER[1:polymer_type-1]) + bead_index # 当前珠子在整个蛋白质链中的索引
            neighbor_count = total_connections[bead_global_index] # 与当前珠子相邻的珠子数量
            flag = true
            x1, x2, x3 = zeros(Int, 3), zeros(Int, 3), zeros(Int, 3)
            
            for neighbor_index in 1:neighbor_count # 遍历与当前珠子相连的所有珠子
                if connection_matrix[bead_global_index, neighbor_index] < bead_global_index # 只考虑当前珠子之前的相邻珠子
                    conn_index = connection_matrix[bead_global_index, neighbor_index] - sum(BEADS_PER_POLYMER[1:polymer_type-1])
                    if flag # 当前珠子的第一个连接：根据其位置和连接关系计算出当前珠子可能的位置范围
                        x2 = protein_positions[conn_index, :] .- connection_lengths[bead_global_index, neighbor_index]
                        x3 = protein_positions[conn_index, :] .+ connection_lengths[bead_global_index, neighbor_index]
                        flag = false
                    else # 不是当前珠子的第一个连接：使用当前珠子和之前相邻珠子的位置范围，更新当前珠子的位置范围
                        x1 = max.(x2, protein_positions[conn_index, :] .- connection_lengths[bead_global_index, neighbor_index])
                        x2 = min.(x3, protein_positions[conn_index, :] .+ connection_lengths[bead_global_index, neighbor_index])
                    end
                end
            end
            
            for dim in 1:3
                protein_positions[bead_index, dim] = rand(1:(x3[dim] - x2[dim] + 1)) + x2[dim] - 1
            end

            if all(!all(protein_positions[bead_index, :] .== protein_positions[prev, :]) for prev in 1:bead_index-1)
                break
            end
        end
    end

    return protein_positions
end

function check_overlap(bead_position_local, bead_reference_position, bead_index_lattice, current_polymers_type) # 检查给定位置上的蛋白质链是否与已存在的蛋白质链或者自身重叠
    bead_absolute_position = zeros(Int, size(bead_position_local))
    for current_bead_index in 1:BEADS_PER_POLYMER[current_polymers_type]
        bead_absolute_position[current_bead_index, :] .= periodic_boundary_3d(bead_position_local[current_bead_index, :] .+ bead_reference_position)
        if bead_index_lattice[bead_absolute_position[current_bead_index, 1], bead_absolute_position[current_bead_index, 2], bead_absolute_position[current_bead_index, 3]] != 0
            return false, bead_absolute_position
        end
        for previous_bead_index in 1:current_bead_index-1
            if all(bead_absolute_position[current_bead_index, :] .== bead_absolute_position[previous_bead_index, :])
                return false, bead_absolute_position
            end
        end
    end
    return true, bead_absolute_position
end

function generate_droplet_position() # 生成液滴中蛋白质链的位置
    # 生成球坐标
	θ, φ, r = rand() * 2π, acos(rand() * 2 - 1), (rand()^(1/3)) * DROPLET_RADIUS
	cartesian_coords = [r * sin(φ) * cos(θ), r * sin(φ) * sin(θ), r * cos(φ)] .+ (LATTICE_SIZE .+ 1) / 2
    
    # 将坐标值取整并返回整数数组
    return floor.(Int, cartesian_coords)
end

function add_to_lattice(bead_absolute_position, current_polymers_type, current_polymers_index, bead_index_lattice, bead_positions_list, polymers_type, max_beads) # 将蛋白质链添加到 lattice 中，并更新相关的数据结构
	num_polymers_per_type = parameters["num_polymers_per_type"]
	polymer_start_index = sum(num_polymers_per_type[1:current_polymers_type-1]) # 计算当前类型之前的所有聚合物数量总和，作为当前类型的起始索引
	for current_bead_index in 1:BEADS_PER_POLYMER[current_polymers_type]
		idx = (current_polymers_index + polymer_start_index - 1) * max_beads + current_bead_index # 当前珠子的全局索引
		bead_index_lattice[bead_absolute_position[current_bead_index, 1], bead_absolute_position[current_bead_index, 2], bead_absolute_position[current_bead_index, 3]] = idx
		bead_positions_list[idx, :] .= bead_absolute_position[current_bead_index, :]
		polymers_type[current_polymers_index+polymer_start_index] = current_polymers_type # 聚合物类型信息
	end
    return bead_index_lattice, bead_positions_list, polymers_type
end

function form_interactions(NUM_POLYMER_TYPES, num_polymers_per_type, BEADS_PER_POLYMER, total_energy, energy, polymers_type, interaction_list, bead_index_lattice, bead_positions_list, neighbor_matrix, parameters) # 蛋白质链间的随机键形成
	# 实现形成随机键
	current_bead_global_index = 0 # 当前珠子的索引
	for current_polymer_type in 1:NUM_POLYMER_TYPES # 遍历所有类型的聚合物
		for current_polymer_index in 1:num_polymers_per_type[current_polymer_type] # 遍历当前类型的所有链
			for bead_index_in_polymer in 1:BEADS_PER_POLYMER[current_polymer_type] # 遍历当前链上的所有珠子
				current_bead_global_index += 1
				current_bead_index = bead_global_to_local_index(current_bead_global_index, parameters, polymers_type) # 珠子的类型
				if interaction_list[current_bead_global_index] == 0 # 当前珠子没有连接
					forward_moves_count = 1 # 前向旋转移动的数量
					interaction_beads_list = zeros(Int, 27) # 存储可以旋转的珠子索引
					for neighbor_index in 1:26 # 查看周围的珠子以确定可能的相互作用
						neighbor_global_position = periodic_boundary_3d(bead_positions_list[current_bead_global_index, :] .+ neighbor_matrix[neighbor_index, :]) # 当前位置的周围坐标
						if bead_index_lattice[neighbor_global_position[1], neighbor_global_position[2], neighbor_global_position[3]] != 0 # 位置 neighbor_global_position 处有珠子
							neighbor_bead_global_index = bead_index_lattice[neighbor_global_position[1], neighbor_global_position[2], neighbor_global_position[3]]
							if energy[current_bead_index, bead_global_to_local_index(neighbor_bead_global_index, parameters, polymers_type)] != 0 && interaction_list[neighbor_bead_global_index] == 0 # 两个珠子有连接，且周围的珠子没有连接
								forward_moves_count += 1
								interaction_beads_list[forward_moves_count] = neighbor_bead_global_index 
							end
						end
					end

					if forward_moves_count > 1
						interaction_bead_global_index = interaction_beads_list[rand(1:forward_moves_count)] # 随机选择一个可以旋转的珠子
						if interaction_bead_global_index > 0
							interaction_list[current_bead_global_index] = interaction_bead_global_index # current_bead_global_index 与随机选择的可以旋转的珠子 interaction_bead_global_index 进行连接
							interaction_list[interaction_bead_global_index] = current_bead_global_index
							
                            current_bead_index = bead_global_to_local_index(current_bead_global_index, parameters, polymers_type)
							interaction_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type)
							total_energy = calculate_energy(total_energy, energy, current_bead_index, interaction_bead_index, parameters, current_bead_global_index, interaction_bead_global_index)
						end
					end
				end
			end
		end
	end
    return interaction_list, total_energy
end
