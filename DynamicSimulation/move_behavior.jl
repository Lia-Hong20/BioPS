function single_bead_move(parameters) #随机选择一个珠子，计算出珠子在x, y, z每个维度上的最小位移和最大位移，分别移动（会改变链的形状）
	bead_positions_list = parameters["bead_positions_list"]
	polymers_type = parameters["polymers_type"]
	total_connections = parameters["total_connections"]
	max_beads = parameters["max_beads"]
	connection_matrix = parameters["connection_matrix"]
	connection_lengths = parameters["connection_lengths"]
	interaction_list = parameters["interaction_list"]
	bead_index_lattice = parameters["bead_index_lattice"]
	neighbor_matrix = parameters["neighbor_matrix"]
	energy = parameters["energy"]
	total_energy = parameters["total_energy"]
	num_polymers_per_type = parameters["num_polymers_per_type"]

	parameters["propose_moves"][1] += 1
	current_polymer_index = rand(1:sum(num_polymers_per_type))  # 随机选择一条链
	current_polymer_type = polymers_type[current_polymer_index]  # 链的类型
	current_bead_global_index = (current_polymer_index - 1) * max_beads + rand(1:BEADS_PER_POLYMER[current_polymer_type])  # 随机选择该链上的一个珠子，珠子索引包括前面所有珠子
	current_bead_index = bead_global_to_local_index(current_bead_global_index, parameters, polymers_type)  # 该珠子在链上的顺序，参照z_Connection
	current_bead_connection_count = total_connections[current_bead_index]  # 该珠子固定键的连接数
	current_bead_position = bead_positions_list[current_bead_global_index, :] # 获取当前珠子的坐标

	connected_bead_displacements = [] # 存储当前珠子连接的其他珠子的位置位移的数组
	connected_bead_lengths = [] # 存储与当前珠子连接的所有珠子的距离
	for connected_bead_index in 1:current_bead_connection_count # 遍历与当前珠子连接的所有珠子
		connected_bead_position = bead_positions_list[current_bead_global_index-current_bead_index+connection_matrix[current_bead_index, connected_bead_index], :] # 获取与当前珠子连接的珠子的坐标
		push!(connected_bead_displacements, distance_periodic_boundary(connected_bead_position, current_bead_position)) # 与要移动的珠子连接的其他珠子的位移
		push!(connected_bead_lengths, connection_lengths[current_bead_index, connected_bead_index])
	end
	connected_bead_displacements = hcat(connected_bead_displacements...)'
	relative_move_position = calculate_relative_move_position(connected_bead_displacements, connected_bead_lengths, current_bead_connection_count)
	new_bead_position = round.(Int, periodic_boundary_3d(current_bead_position + relative_move_position))

	flag = bead_index_lattice[new_bead_position[1], new_bead_position[2], new_bead_position[3]] == 0 || bead_index_lattice[new_bead_position[1], new_bead_position[2], new_bead_position[3]] == current_bead_global_index

	if flag
		energy_before_move = 0.0
		interaction_bead_global_index = interaction_list[current_bead_global_index]
		if interaction_bead_global_index > 0  # 当前珠子已经与其他珠子连接
			interaction_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type)
			energy_before_move = calculate_energy(energy_before_move, energy, current_bead_index, interaction_bead_index, parameters, current_bead_global_index, interaction_bead_global_index)
		end

		# 移动前相互作用的邻居数
		pre_move_interacting_neighbors_count, _ = count_interacting_neighbors(current_bead_position, bead_index_lattice, energy, interaction_list, current_bead_index, current_bead_global_index, parameters, polymers_type, neighbor_matrix, false)
		# 移动后相互作用的邻居数
		post_move_interacting_neighbors_count, interaction_beads_list =
			count_interacting_neighbors(new_bead_position, bead_index_lattice, energy, interaction_list, current_bead_index, current_bead_global_index, parameters, polymers_type, neighbor_matrix, true)

		flag = true
	else
		flag = false
	end

	# Accept or Reject the move
	if flag
		interaction_bead_global_index = interaction_beads_list[rand(1:post_move_interacting_neighbors_count)] # 随机选择一个移动后可以相互作用的珠子
		energy_after_move = 0.0
		if interaction_bead_global_index != 0  # 选择的位置有效
			interaction_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type) # 移动后相邻珠子的类型
			energy_after_move = calculate_energy(energy_after_move, energy, current_bead_index, interaction_bead_index, parameters, current_bead_global_index, interaction_bead_global_index)
		end

		if rand() < exp(-(energy_after_move - energy_before_move) / kBT) * (post_move_interacting_neighbors_count / pre_move_interacting_neighbors_count) * parameters["adjustment_factor"] # Metropolis准则
			total_energy += (energy_after_move - energy_before_move)  # 更新系统的总能量

			bead_index_lattice[bead_positions_list[current_bead_global_index, 1], bead_positions_list[current_bead_global_index, 2], bead_positions_list[current_bead_global_index, 3]] = 0
			bead_index_lattice[new_bead_position[1], new_bead_position[2], new_bead_position[3]] = current_bead_global_index
			bead_positions_list[current_bead_global_index, :] .= new_bead_position  # 将当前珠子从旧位置移除并放置到新位置
			if interaction_list[current_bead_global_index] > 0
				interaction_list[interaction_list[current_bead_global_index]] = 0
			end
			interaction_list[current_bead_global_index] = interaction_bead_global_index
			if interaction_bead_global_index != 0
				interaction_list[interaction_bead_global_index] = current_bead_global_index
			end
		else
			parameters["rejection_moves"][1] += 1 #成功计算Metropolis准则但拒绝移动的次数
		end
	else
		parameters["rejection_moves"][1] += 1 #在计算Metropolis准则之前就已经判断无法移动的次数
	end
	parameters["total_energy"] = total_energy
	parameters["bead_index_lattice"] = bead_index_lattice
	parameters["bead_positions_list"] = bead_positions_list
	parameters["interaction_list"] = interaction_list
	return parameters
end

function paired_bead_move(parameters) # 随机选择一个珠子，得到与该珠子相互作用的珠子，即移动一对相互作用的珠子
	bead_positions_list = parameters["bead_positions_list"]
	polymers_type = parameters["polymers_type"]
	total_connections = parameters["total_connections"]
	max_beads = parameters["max_beads"]
	connection_matrix = parameters["connection_matrix"]
	connection_lengths = parameters["connection_lengths"]
	interaction_list = parameters["interaction_list"]
	bead_index_lattice = parameters["bead_index_lattice"]
	num_polymers_per_type = parameters["num_polymers_per_type"]

	parameters["propose_moves"][2] += 1 # 记录尝试移动的次数

	current_polymer_index = rand(1:sum(num_polymers_per_type)) # 随机选择一条蛋白质链
	current_polymer_type = polymers_type[current_polymer_index] # 获取蛋白质链的类型
	current_bead_global_index = (current_polymer_index - 1) * max_beads + rand(1:BEADS_PER_POLYMER[current_polymer_type])  # 随机选择该链上的一个珠子，珠子索引包括前面所有珠子
	current_bead_index = bead_global_to_local_index(current_bead_global_index, parameters, polymers_type) # 当前珠子在该链上的序号
	current_bead_position = bead_positions_list[current_bead_global_index, :] # 当前珠子的位置
	interaction_bead_global_index = interaction_list[current_bead_global_index] # 与当前珠子相互作用的珠子的编号

	# 如果没有相互作用的珠子，则拒绝移动
	if interaction_bead_global_index == 0
		parameters["rejection_moves"][2] += 1  # 记录拒绝的移动次数
		return parameters
	end

	# 如果有相互作用的珠子，继续计算
	connection_count = 0  # 链接的珠子数
	flag = true

	# 遍历与当前珠子连接的所有固定键，并计算这些键对应的其他珠子的位置和位移
	# 处理与 current_bead 相关的连接
	connected_bead_displacements, connected_bead_lengths, connection_count_current =
		process_connections(current_bead_global_index, current_bead_index, interaction_bead_global_index, total_connections, connection_matrix, bead_positions_list, current_bead_position, connection_lengths)
	connection_count += connection_count_current

	# 当前珠子有连接
	interaction_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type)
	interaction_bead_position = bead_positions_list[interaction_bead_global_index, :]
	connected_interaction_bead_displacements, connected_interaction_bead_lengths, connection_count_interaction =
		process_connections(interaction_bead_global_index, interaction_bead_index, current_bead_global_index, total_connections, connection_matrix, bead_positions_list, interaction_bead_position, connection_lengths)
	append!(connected_bead_displacements, connected_interaction_bead_displacements)
	append!(connected_bead_lengths, connected_interaction_bead_lengths)
	connection_count += connection_count_interaction

	# 当前珠子的目标位置向量
	connected_bead_displacements = hcat(connected_bead_displacements...)'
	relative_move_position = calculate_relative_move_position(connected_bead_displacements, connected_bead_lengths, connection_count)
	new_bead_position = round.(Int, periodic_boundary_3d(current_bead_position + relative_move_position))

	# 检查新位置是否发生碰撞
	if bead_index_lattice[new_bead_position...] != 0 && bead_index_lattice[new_bead_position...] != interaction_bead_global_index
		flag = false
	end
	new_interaction_bead_position = round.(Int, periodic_boundary_3d(interaction_bead_position + relative_move_position))
	if bead_index_lattice[new_interaction_bead_position...] != 0 && bead_index_lattice[new_interaction_bead_position...] != current_bead_global_index
		flag = false
	end

	if flag #移动可行,更新变量
		bead_index_lattice[bead_positions_list[current_bead_global_index, 1], bead_positions_list[current_bead_global_index, 2], bead_positions_list[current_bead_global_index, 3]] = 0
		bead_positions_list[current_bead_global_index, :] = new_bead_position
		bead_index_lattice[new_bead_position[1], new_bead_position[2], new_bead_position[3]] = current_bead_global_index

		bead_index_lattice[bead_positions_list[interaction_bead_global_index, 1], bead_positions_list[interaction_bead_global_index, 2], bead_positions_list[interaction_bead_global_index, 3]] = 0
		bead_positions_list[interaction_bead_global_index, :] = new_interaction_bead_position
		bead_index_lattice[new_interaction_bead_position[1], new_interaction_bead_position[2], new_interaction_bead_position[3]] = interaction_bead_global_index
	else #如果移动不可行，则拒绝+1
		parameters["rejection_moves"][2] += 1
	end
	parameters["bead_index_lattice"] = bead_index_lattice
	parameters["bead_positions_list"] = bead_positions_list
	return parameters
end

function reconnect_move(parameters) # 更改与珠子MS相互作用的珠子
	bead_positions_list = parameters["bead_positions_list"]
	polymers_type = parameters["polymers_type"]
	max_beads = parameters["max_beads"]
	interaction_list = parameters["interaction_list"]
	bead_index_lattice = parameters["bead_index_lattice"]
	energy = parameters["energy"]
	neighbor_matrix = parameters["neighbor_matrix"]
	total_energy = parameters["total_energy"]
	rotatable_beads_list = parameters["rotatable_beads_list"]
	num_polymers_per_type = parameters["num_polymers_per_type"]

	parameters["propose_moves"][3] += 1
	energy_after_move = 0.0
	energy_before_move = 0.0

	current_polymer_index = rand(1:sum(num_polymers_per_type)) # 随机选择一条蛋白质链
	current_polymer_type = polymers_type[current_polymer_index]
	current_bead_global_index = (current_polymer_index - 1) * max_beads + rand(1:BEADS_PER_POLYMER[current_polymer_type])  # 随机选择链上的一个珠子
	current_bead_index = bead_global_to_local_index(current_bead_global_index, parameters, polymers_type)
	current_bead_position = bead_positions_list[current_bead_global_index, :]

	if rotatable_beads_list[current_bead_index] != 0
		if interaction_list[current_bead_global_index] > 0 # 当前珠子有一个与之连接的珠子
			interaction_bead_global_index = interaction_list[current_bead_global_index]
			interaction_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type)
			energy_before_move = calculate_energy(energy_before_move, energy, current_bead_index, interaction_bead_index, parameters, current_bead_global_index, interaction_bead_global_index)
		end

		post_move_interacting_neighbors, interaction_beads_list =
			count_interacting_neighbors(current_bead_position, bead_index_lattice, energy, interaction_list, current_bead_index, current_bead_global_index, parameters, polymers_type, neighbor_matrix, false)

		interaction_bead_global_index = interaction_beads_list[rand(1:post_move_interacting_neighbors)] # 从可选的旋转移动中随机选择一个珠子编号
		if interaction_bead_global_index != 0
			interaction_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type)
			energy_after_move = calculate_energy(energy_after_move, energy, current_bead_index, interaction_bead_index, parameters, current_bead_global_index, interaction_bead_global_index)
		end

		if rand() < exp(-(energy_after_move - energy_before_move) / kBT) * parameters["adjustment_factor"]
			total_energy += energy_after_move - energy_before_move
			if interaction_list[current_bead_global_index] > 0
				interaction_list[interaction_list[current_bead_global_index]] = 0
			end
			interaction_list[current_bead_global_index] = interaction_bead_global_index
			if interaction_bead_global_index != 0
				interaction_list[interaction_bead_global_index] = current_bead_global_index
			end
		else
			parameters["rejection_moves"][3] += 1
		end
	end
	parameters["total_energy"] = total_energy
	parameters["interaction_list"] = interaction_list
	return parameters
end

function reptation_move(parameters) #随机选择一个移动长度SlitherTimes，聚合物向前（向后）滑动SlitherTimes个珠子，后面的珠子继承前面珠子的坐标，前面的珠子在移动范围内随机选择一个坐标
	bead_positions_list = parameters["bead_positions_list"]
	polymers_type = parameters["polymers_type"]
	max_beads = parameters["max_beads"]
	interaction_list = parameters["interaction_list"]
	bead_index_lattice = parameters["bead_index_lattice"]
	energy = parameters["energy"]
	neighbor_matrix = parameters["neighbor_matrix"]
	connection_lengths = parameters["connection_lengths"]
	total_energy = parameters["total_energy"]
	num_polymers_per_type = parameters["num_polymers_per_type"]

	parameters["propose_moves"][5] += 1
	slither_position_list = zeros(Int, max_beads, 3) #存储滑移链的位置

	# 找到特定类型的蛋白质链
	current_polymer_type = rand(1:NUM_POLYMER_TYPES) # 选择蛋白质类型
	current_polymer_index = rand(1:num_polymers_per_type[current_polymer_type]) # 随机选择一条该类型的聚合物链
	current_bead_global_index_reference = (current_polymer_index - 1) * max_beads  # 全局
	current_bead_local_index_reference = sum(BEADS_PER_POLYMER[1:current_polymer_type-1]) # 当前链

	flag = false
	if BEADS_PER_POLYMER[current_polymer_type] > 5 #确保滑动前移步数不超过5个珠子
		slither_length = rand(1:5)  #如果蛋白质链的长度超过5个珠子，SlitherTimes将在1到5之间随机选择
	else
		slither_length = rand(1:BEADS_PER_POLYMER[current_polymer_type]-1) #SlitherTimes将在1到Vk[current_polymer_type]-1之间随机选择
	end

	if rand() < 0.5 #尾带着走
		slither_forward = 1
		for current_bead_index in 1:slither_forward:BEADS_PER_POLYMER[current_polymer_type]-slither_length
			slither_position_list[current_bead_index, :] = bead_positions_list[current_bead_global_index_reference+slither_length+current_bead_index, :] #滑移操作的潜在终点位置的坐标
		end
		criterion_star = current_bead_global_index_reference
		criterion_end = current_bead_global_index_reference + slither_length
		for current_bead_index in BEADS_PER_POLYMER[current_polymer_type]-slither_length+1:slither_forward:BEADS_PER_POLYMER[current_polymer_type] #从上一个循环结束的位置遍历直到滑移操作的最后一个珠子
			check_range = 1:current_bead_index-1
			flag = attempt_bead_move(
				slither_position_list[current_bead_index-1, :],
				connection_lengths[current_bead_local_index_reference+current_bead_index, 1],
				bead_index_lattice,
				slither_position_list,
				current_bead_index,
				criterion_star,
				criterion_end,
				check_range,
				flag,
			)
			if flag
				break
			end
		end
	else  #头带着走
		slither_forward = -1
		for current_bead_index in BEADS_PER_POLYMER[current_polymer_type]:slither_forward:slither_length+1 #从滑移操作的末尾位置开始，逆序遍历直到滑移操作的起始位置
			slither_position_list[current_bead_index, :] = bead_positions_list[current_bead_global_index_reference-slither_length+current_bead_index, :]
		end
		for current_bead_index in slither_length:slither_forward:1 #逆序遍历滑移操作的每个珠子
			check_range = current_bead_index+1:BEADS_PER_POLYMER[current_polymer_type]
			criterion_star = current_bead_global_index_reference + BEADS_PER_POLYMER[current_polymer_type] - slither_length
			criterion_end = current_bead_global_index_reference + BEADS_PER_POLYMER[current_polymer_type]
			flag = attempt_bead_move(
				slither_position_list[current_bead_index+1, :],
				connection_lengths[current_bead_local_index_reference+current_bead_index, 1],
				bead_index_lattice,
				slither_position_list,
				current_bead_index,
				criterion_star,
				criterion_end,
				check_range,
				flag,
			)
			if flag
				break
			end
		end
	end

	if !flag
		rows = [Tuple(slither_position_list[i, :]) for i in 1:size(slither_position_list, 1)] # 全局检查是否有位置重叠
		has_duplicate_rows = length(rows) != length(unique(rows))
		if has_duplicate_rows
			parameters["rejection_moves"][5] += 1
		else
			before_log_bind_count, after_log_bind_count, energy_before_move, energy_after_move, slither_rotatable_list =
				calculate_energy_and_Rosenbluth(current_polymer_type, slither_position_list, interaction_list, current_polymer_index, max_beads, energy, parameters, neighbor_matrix, bead_positions_list, polymers_type, bead_index_lattice)

			# Accept or reject the move
			acceptance_prob = exp(after_log_bind_count - before_log_bind_count - (energy_after_move - energy_before_move) / kBT) * parameters["adjustment_factor"]
			if rand() < acceptance_prob
				total_energy += energy_after_move - energy_before_move
				for current_bead_index in 1:BEADS_PER_POLYMER[current_polymer_type]  # 清空旧绑定
					global_index = current_bead_global_index_reference + current_bead_index
					bound_bead = interaction_list[global_index]
					if bound_bead != 0
						interaction_list[bound_bead] = 0
					end
					bead_index_lattice[bead_positions_list[global_index, 1], bead_positions_list[global_index, 2], bead_positions_list[global_index, 3]] = 0  # 清空老坐标
				end
				for current_bead_index ∈ 1:BEADS_PER_POLYMER[current_polymer_type]  #新的绑定
					global_index = current_bead_global_index_reference + current_bead_index
					slither_bead = slither_rotatable_list[current_bead_index]
					interaction_list[global_index] = slither_bead
					if slither_bead != 0
						interaction_list[slither_bead] = global_index
					end
					new_position = slither_position_list[current_bead_index, :]
					bead_index_lattice[new_position[1], new_position[2], new_position[3]] = global_index
					bead_positions_list[global_index, :] = new_position
				end
			else
				parameters["rejection_moves"][5] += 1
			end
		end
	end
	parameters["interaction_list"] = interaction_list
	parameters["total_energy"] = total_energy
	parameters["bead_index_lattice"] = bead_index_lattice
	parameters["bead_positions_list"] = bead_positions_list
	return parameters
end

function chain_translation_move(parameters)
	bead_positions_list = parameters["bead_positions_list"]
	polymers_type = parameters["polymers_type"]
	max_beads = parameters["max_beads"]
	interaction_list = parameters["interaction_list"]
	bead_index_lattice = parameters["bead_index_lattice"]
	total_energy = parameters["total_energy"]
	energy = parameters["energy"]
	neighbor_matrix = parameters["neighbor_matrix"]
	num_polymers_per_type = parameters["num_polymers_per_type"]

	parameters["propose_moves"][8] += 1
	current_polymer_index = rand(1:sum(num_polymers_per_type))  # 随机选择一条链
	current_polymer_type = polymers_type[current_polymer_index]  # 链的类型

	# 生成一个在 L/2 半径内的随机位移
	move_vector = [rand(-LATTICE_SIZE[i]÷2:LATTICE_SIZE[i]÷2) for i in 1:3]
	if all(move_vector .== 0)
		parameters["rejection_moves"][8] += 1
		return parameters
	end

	# 检查该位移下的链是否发生碰撞
	collision_flag = false
	bead_range_start = (current_polymer_index - 1) * BEADS_PER_POLYMER[current_polymer_type] + 1
	bead_range_end = bead_range_start + BEADS_PER_POLYMER[current_polymer_type] - 1

	translation_position_list = bead_positions_list[bead_range_start:bead_range_end, :]
	for bead_index in bead_range_start:bead_range_end
		new_position = periodic_boundary_3d(bead_positions_list[bead_index, :] .+ move_vector)
		if bead_index_lattice[new_position[1], new_position[2], new_position[3]] == 0
			translation_position_list[bead_index-bead_range_start+1, :] = new_position
		else
			collision_flag = true
			parameters["rejection_moves"][8] += 1
			break
		end
	end

	if !collision_flag # 如果没有碰撞，则继续计算能量和接受概率
		before_log_bind_count, after_log_bind_count, energy_before_move, energy_after_move, translation_rotatable_list =
			calculate_energy_and_Rosenbluth(current_polymer_type, translation_position_list, interaction_list, current_polymer_index, max_beads, energy, parameters, neighbor_matrix, bead_positions_list, polymers_type, bead_index_lattice)

		# Accept or reject the move
		acceptance_prob = exp(after_log_bind_count - before_log_bind_count + (-(energy_after_move - energy_before_move)) / kBT) * parameters["adjustment_factor"]
		if rand() < acceptance_prob
			total_energy += energy_after_move - energy_before_move
			for current_bead_index in 1:BEADS_PER_POLYMER[current_polymer_type]  # 清空旧绑定
				global_index = bead_range_start - 1 + current_bead_index
				bound_bead = interaction_list[global_index]
				if bound_bead != 0
					interaction_list[bound_bead] = 0
				end
				bead_index_lattice[bead_positions_list[global_index, 1], bead_positions_list[global_index, 2], bead_positions_list[global_index, 3]] = 0  # 清空老坐标
			end
			for current_bead_index ∈ 1:BEADS_PER_POLYMER[current_polymer_type]  #新的绑定
				global_index = bead_range_start - 1 + current_bead_index
				translation_bead = translation_rotatable_list[current_bead_index]
				interaction_list[global_index] = translation_bead
				if translation_bead != 0
					interaction_list[translation_bead] = global_index
				end
				new_position = translation_position_list[current_bead_index, :]
				bead_index_lattice[new_position[1], new_position[2], new_position[3]] = global_index
				bead_positions_list[global_index, :] = new_position
			end
		else
			parameters["rejection_moves"][8] += 1
		end
	end
	parameters["interaction_list"] = interaction_list
	parameters["total_energy"] = total_energy
	parameters["bead_index_lattice"] = bead_index_lattice
	parameters["bead_positions_list"] = bead_positions_list
	return parameters
end

function chain_formed_cluster_move(parameters)  # 随机选择一条链，检查共价/可逆键，在蛋白质链的数量小于六的团簇中，将团簇移动到晶格中的一个随机新位置
	bead_positions_list = parameters["bead_positions_list"]
	polymers_type = parameters["polymers_type"]
	max_beads = parameters["max_beads"]
	bead_index_lattice = parameters["bead_index_lattice"]
	interaction_list = parameters["interaction_list"]
	num_polymers_per_type = parameters["num_polymers_per_type"]

	parameters["propose_moves"][7] += 1
	current_polymer_index = rand(1:sum(num_polymers_per_type)) # 随机选择一条链
	current_polymer_type = polymers_type[current_polymer_index] # 链的类型

	cluster_list = [current_polymer_index]  # 包含选定蛋白质链的集群列表
	is_polymer_in_cluster = trues(sum(num_polymers_per_type))  # 标记蛋白质链不属于任何集群
	is_polymer_in_cluster[current_polymer_index] = false  # 标记蛋白质链已属于当前集群

	cluster_index = 1
	while cluster_index <= length(cluster_list) && length(cluster_list) <= 5 # 集群中蛋白质链的数量不超过5，当前正在处理的集群没有超出集群总数
		current_polymer_index = cluster_list[cluster_index]
		current_polymer_type = polymers_type[current_polymer_index]
		bead_start_index = (current_polymer_index - 1) * max_beads + 1 # 当前蛋白质链的第一个珠子的全局索引
		bead_end_index = bead_start_index + BEADS_PER_POLYMER[current_polymer_type] - 1 # 当前蛋白质链的最后一个珠子的全局索引

		for current_bead_index in bead_start_index:bead_end_index # 遍历当前链上的所有珠子
			interaction_index = interaction_list[current_bead_index] # 当前珠子相互作用的珠子
			if interaction_index != 0
				interacting_polymer_index = bead_index_to_polymer_index(interaction_index, parameters) # 相互作用的珠子所在的链
				if is_polymer_in_cluster[interacting_polymer_index] # 该链不属于任何集群
					push!(cluster_list, interacting_polymer_index)
					is_polymer_in_cluster[interacting_polymer_index] = false
				end
			end
		end
		cluster_index += 1
	end

	if length(cluster_list) <= 5
		relative_move_position = rand_int_array(LATTICE_SIZE) # 生成随机移动向量
		if all(relative_move_position .== 0) # 发生了移动
			parameters["rejection_moves"][7] += 1
		else
			is_move = true
			new_positions = Vector{Tuple{Int, Int, Int}}() # 创建一个空的向量

			for polymer_index in cluster_list # 遍历集群中的所有蛋白质链
				polymer_type = polymers_type[polymer_index]
				bead_start_index = (polymer_index - 1) * BEADS_PER_POLYMER[polymer_type] + 1
				bead_end_index = bead_start_index + BEADS_PER_POLYMER[polymer_type] - 1

				for current_bead_index in bead_start_index:bead_end_index # 遍历当前链上的所有珠子
					new_bead_position = periodic_boundary_3d(bead_positions_list[current_bead_index, :] .+ relative_move_position)
					if bead_index_lattice[new_bead_position[1], new_bead_position[2], new_bead_position[3]] != 0 # 如果新位置有其他珠子存在，跳出循环
						is_move = false
						break
					else
						push!(new_positions, (new_bead_position[1], new_bead_position[2], new_bead_position[3]))
					end
				end
				if !is_move
					break
				end
			end

			if is_move # 更新变量
				for (i, polymer_index) in enumerate(cluster_list) # 遍历集群中的所有蛋白质链
					polymer_type = polymers_type[polymer_index]
					bead_start_index = (polymer_index - 1) * BEADS_PER_POLYMER[polymer_type] + 1
					bead_end_index = bead_start_index + BEADS_PER_POLYMER[polymer_type] - 1

					for (current_bead_index, new_bead_position) in zip(bead_start_index:bead_end_index, new_positions[(i-1)*BEADS_PER_POLYMER[polymer_type]+1:i*BEADS_PER_POLYMER[polymer_type]])
						bead_index_lattice[bead_positions_list[current_bead_index, 1], bead_positions_list[current_bead_index, 2], bead_positions_list[current_bead_index, 3]] = 0
						bead_positions_list[current_bead_index, :] = collect(new_bead_position)
						bead_index_lattice[new_bead_position[1], new_bead_position[2], new_bead_position[3]] = current_bead_index
					end
				end
			else
				parameters["rejection_moves"][7] += 1
			end
		end
	else
		parameters["rejection_moves"][7] += 1
	end

	parameters["bead_index_lattice"] = bead_index_lattice
	parameters["bead_positions_list"] = bead_positions_list
	return parameters
end

function non_largest_cluster_move(parameters) # 随机选择一个不是最大的集群，随机得到移动向量，更新集群中的每个珠子的坐标。每次移动，从选择随机集群开始重复10次
	bead_positions_list = parameters["bead_positions_list"]
	polymers_type = parameters["polymers_type"]
	max_beads = parameters["max_beads"]
	bead_index_lattice = parameters["bead_index_lattice"]

	parameters["propose_moves"][6] += 1

	num_clusters, cluster_ranges, clustered_polymers = identify_protein_clusters(parameters)

	if num_clusters > 1 # 集群数量大于1
		move_attempts = min(num_clusters - 1, 10)  # 移动次数
		available_clusters = num_clusters - 1
		max_cluster_index = argmax(cluster_ranges[1:num_clusters, 2] .- cluster_ranges[1:num_clusters, 1] .+ 1) # 获取最大集群的索引

		Fisher = [i for i in 1:num_clusters if i != max_cluster_index] # 包含除最大集群外的所有集群索引的列表
		for attempt in 1:move_attempts
			selected_cluster_idx = rand(1:available_clusters)
			current_cluster_index = Fisher[selected_cluster_idx] # 从列表中随机选择一个集群
			Fisher[selected_cluster_idx] = Fisher[available_clusters]
			available_clusters -= 1

			collision_flag = false
			move_vector = rand_int_array(LATTICE_SIZE)  # 生成一个随机的移动向量，指示蛋白质的移动方向和距离

			if all(move_vector .== 0)
				parameters["rejection_moves"][6] += 1
				continue
			end

			cluster_start, cluster_end = cluster_ranges[current_cluster_index, :]
			cluster_polymers = clustered_polymers[cluster_start:cluster_end]

			for cluster_polymer in cluster_polymers  # 检查移动是否有重叠，遍历集群中的所有蛋白质
				polymer_type = polymers_type[cluster_polymer]
				bead_range_start = (cluster_polymer - 1) * max_beads + 1
				bead_range_end = bead_range_start + BEADS_PER_POLYMER[polymer_type] - 1

				for bead_index in bead_range_start:bead_range_end #遍历所有的珠子
					new_position = round.(Int, periodic_boundary_3d(bead_positions_list[bead_index, :] .+ move_vector))
					if bead_index_lattice[new_position[1], new_position[2], new_position[3]] != 0
						collision_flag = true
						break
					end
				end
				if collision_flag
					break
				end
			end

			if !collision_flag #更新珠子位置
				for cluster_polymer in cluster_polymers # 清空老位置
					polymer_type = polymers_type[cluster_polymer]
					bead_range_start = (cluster_polymer - 1) * BEADS_PER_POLYMER[polymer_type] + 1
					bead_range_end = bead_range_start + BEADS_PER_POLYMER[polymer_type] - 1

					for bead_index in bead_range_start:bead_range_end
						bead_index_lattice[bead_positions_list[bead_index, :][1], bead_positions_list[bead_index, :][2], bead_positions_list[bead_index, :][3]] = 0
					end
				end

				for cluster_polymer in cluster_polymers # 在新位置上放入珠子
					polymer_type = polymers_type[cluster_polymer]
					bead_range_start = (cluster_polymer - 1) * BEADS_PER_POLYMER[polymer_type] + 1
					bead_range_end = bead_range_start + BEADS_PER_POLYMER[polymer_type] - 1

					for bead_index in bead_range_start:bead_range_end
						new_position = periodic_boundary_3d(bead_positions_list[bead_index, :] .+ move_vector)
						bead_positions_list[bead_index, :] = new_position
						bead_index_lattice[new_position[1], new_position[2], new_position[3]] = bead_index
					end
				end
			else
				parameters["rejection_moves"][6] += 1
			end
		end
	else
		parameters["rejection_moves"][6] += 10
	end

	parameters["bead_index_lattice"] = bead_index_lattice
	parameters["bead_positions_list"] = bead_positions_list
	return parameters
end

function double_pivot_move(parameters)
	bead_positions_list = parameters["bead_positions_list"]
	polymers_type = parameters["polymers_type"]
	max_beads = parameters["max_beads"]
	connection_lengths = parameters["connection_lengths"]
	interaction_list = parameters["interaction_list"]
	bead_index_lattice = parameters["bead_index_lattice"]
	total_energy = parameters["total_energy"]
	energy = parameters["energy"]
	num_polymers_per_type = parameters["num_polymers_per_type"]

	parameters["propose_moves"][4] += 1

	current_polymer_index = rand(1:sum(num_polymers_per_type))  # 随机选择一条链
	current_polymer_type = polymers_type[current_polymer_index]  # 链的类型
	current_bead_global_index = (current_polymer_index - 1) * max_beads + rand(1:BEADS_PER_POLYMER[current_polymer_type])  # 随机选择该链上的一个珠子，珠子索引包括前面所有珠子

	# 珠子是否在链的开头或末端
	bead_start_index = (current_polymer_index - 1) * BEADS_PER_POLYMER[current_polymer_type] + 1
	bead_end_index = bead_start_index + BEADS_PER_POLYMER[current_polymer_type] - 1
	if bead_start_index == current_bead_global_index || bead_end_index == current_bead_global_index
		parameters["rejection_moves"][4] += 1
		return parameters
	end

	# 当前珠子和下一个珠子在链中的位置
	current_bead_index = bead_global_to_local_index(current_bead_global_index, parameters, polymers_type)  # 该珠子在链上的顺序，参照z_Connection
	next_bead_global_index = current_bead_global_index + 1
	next_bead_index = current_bead_index + 1

	# 邻居半径
	search_radius = connection_lengths[current_bead_index, 2]
	current_bead_position = bead_positions_list[current_bead_global_index, :] # 当前珠子的坐标

	available_beads_list = []
	available_beads_position_list = []
	# 找到邻居中不在当前链条的珠子
	candidate_beads_list = double_pivot_neighbor_search(current_bead_global_index, current_bead_position, search_radius)
	candidate_beads_list = unique(candidate_beads_list)
	neighbor_num = length(candidate_beads_list)
	for i in 1:neighbor_num
		candidate_beads_polymer_index = bead_index_to_polymer_index(candidate_beads_list[i], parameters)
		candidate_bead_index = bead_global_to_local_index(candidate_beads_list[i], parameters, polymers_type)
		if candidate_beads_polymer_index != current_polymer_index && candidate_bead_index == next_bead_index # 与当前珠子不在同一条链，且索引为当前珠子+1
			append!(available_beads_list, candidate_beads_list[i])
		end
	end

	# 如果没有可行的前进移动，则拒绝
	pre_move_available_beads_count = length(available_beads_list)
	if pre_move_available_beads_count == 0
		parameters["rejection_moves"][4] += 1
		return parameters
	end

	# 随机选择一个候选者
	selected_next_bead_global_index = available_beads_list[rand(1:pre_move_available_beads_count)]
	selected_next_bead_index = bead_global_to_local_index(selected_next_bead_global_index, parameters, polymers_type)
	selected_this_bead_global_index = selected_next_bead_global_index - 1
	selected_this_bead_polymer_index = bead_index_to_polymer_index(selected_this_bead_global_index, parameters)
	selected_this_bead_index = bead_global_to_local_index(selected_this_bead_global_index, parameters, polymers_type)
	selected_this_bead_position = bead_positions_list[selected_this_bead_global_index, :] # 当前珠子的坐标
	search_radius = connection_lengths[selected_this_bead_index, 2]

	available_beads_list = []
	available_beads_position_list = []
	# 找到selected_this_bead_global_index邻居中不在当前链条的珠子，计算反向因子
	candidate_beads_list = double_pivot_neighbor_search(selected_this_bead_global_index, selected_this_bead_position, search_radius)
	candidate_beads_list = unique(candidate_beads_list)
	neighbor_num = length(candidate_beads_list)
	for i in 1:neighbor_num
		candidate_beads_polymer_index = bead_index_to_polymer_index(candidate_beads_list[i], parameters)
		candidate_bead_index = bead_global_to_local_index(candidate_beads_list[i], parameters, polymers_type)
		if candidate_beads_polymer_index != selected_this_bead_polymer_index && candidate_bead_index == selected_next_bead_index # 与当前珠子不在同一条链，且索引为当前珠子+1
			append!(available_beads_list, candidate_beads_list[i])
		end
	end
	post_move_available_beads_count = length(available_beads_list)

	# 如果反向候选者数量为0，则进行链交换
	if post_move_available_beads_count == 0
		for j in 0:(bead_end_index-next_bead_global_index)
			temporary_bead_positions_list, temporary_interaction_list, temporary_bead_index_lattice = double_pivot_swap_beads(next_bead_global_index + j, selected_next_bead_global_index + j, bead_positions_list, interaction_list, bead_index_lattice)
		end
	else
		# 计算能量
		energy_before_move = 0.0
		interaction_bead_global_index = interaction_list[current_bead_global_index]
		if interaction_bead_global_index > 0  # 当前珠子已经与其他珠子连接
			interaction_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type)
			energy_before_move = calculate_energy(energy_before_move, energy, current_bead_index, interaction_bead_index, parameters, current_bead_global_index, interaction_bead_global_index)
		end

		# 链交换
		temporary_bead_positions_list = bead_positions_list
		temporary_interaction_list = interaction_list
		temporary_bead_index_lattice = bead_index_lattice
		for j in 0:(bead_end_index-next_bead_global_index)
			temporary_bead_positions_list, temporary_interaction_list, temporary_bead_index_lattice = double_pivot_swap_beads(next_bead_global_index + j, selected_next_bead_global_index + j, bead_positions_list, interaction_list, bead_index_lattice)
		end

		interaction_bead_global_index = temporary_interaction_list[selected_this_bead_global_index]
		energy_after_move = 0.0
		if interaction_bead_global_index != 0  # 选择的位置有效
			interaction_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type) # 移动后相邻珠子的类型
			energy_after_move = calculate_energy(energy_after_move, energy, selected_this_bead_index, interaction_bead_index, parameters, selected_this_bead_global_index, interaction_bead_global_index)
		end

		# Metropolis-Hastings acceptance
		if rand() < exp(-(energy_after_move - energy_before_move) / kBT) * (pre_move_available_beads_count / post_move_available_beads_count) * parameters["adjustment_factor"] # Metropolis准则
			parameters["total_energy"] = total_energy + (energy_after_move - energy_before_move)  # 更新系统的总能量
			parameters["bead_index_lattice"] = temporary_bead_index_lattice
			parameters["bead_positions_list"] = temporary_bead_positions_list
			parameters["interaction_list"] = temporary_interaction_list
		else
			parameters["rejection_moves"][4] += 1
		end
	end
	return parameters
end