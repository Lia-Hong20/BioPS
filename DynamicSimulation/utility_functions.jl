function rand_int_array(j::Vector{Int})
	return [rand(1:j[i]) for i in 1:length(j)]
end

function periodic_boundary_1d(j1::Int, i1::Int) # 一维周期性边界条件
	return mod(j1 + LATTICE_SIZE[i1] - 1, LATTICE_SIZE[i1]) + 1
end

function periodic_boundary_3d(j3::Vector{Int}) # 三维周期性边界条件
	return mod.(j3 .+ LATTICE_SIZE .- 1, LATTICE_SIZE) .+ 1
end

function distance_periodic_boundary(x::Vector{Int}, y::Vector{Int}) # 带有周期性边界条件的两点之间的距离
	distance = mod.(x .- y .+ (LATTICE_SIZE .÷ 2) .+ LATTICE_SIZE, LATTICE_SIZE) .- LATTICE_SIZE .÷ 2
	return Int.(distance)
end

function real_distance_periodic_boundary(x::Vector{Int}, y::Vector{Int}) # 带有周期性边界条件的真实距离
	return mod.(x .- y .+ LATTICE_SIZE .÷ 2 .+ LATTICE_SIZE, Float64.(LATTICE_SIZE)) .- LATTICE_SIZE .÷ 2
end

function bead_index_to_polymer_index(j1::Int, parameters::Dict{String, Any}) # 根据珠子索引确定珠子所属的蛋白质链的索引
	return div(j1 - 1, parameters["max_beads"]) + 1
end

function bead_global_to_local_index(j1::Int, parameters::Dict{String, Any}, polymers_type::Vector{Int}) # 确定给定珠子的索引
	j2 = polymers_type[bead_index_to_polymer_index(j1, parameters)]
	return sum(BEADS_PER_POLYMER[1:j2-1]) + mod(j1 - 1, parameters["max_beads"]) + 1
end

# 记录可以相互作用的相邻珠子
function count_interacting_neighbors(position, bead_index_lattice, energy, interaction_list, current_bead_index, current_bead_global_index, parameters, polymers_type, neighbor_matrix, exclude_current_bead)
	interacting_neighbors = 1
	interaction_beads_list = zeros(Int, 27)

	for interaction_bead_index in 1:26
		neighbor_position = round.(Int, periodic_boundary_3d(position .+ neighbor_matrix[interaction_bead_index, :]))
		neighbor_bead_index = bead_index_lattice[neighbor_position[1], neighbor_position[2], neighbor_position[3]]

		if neighbor_bead_index != 0 && (!exclude_current_bead || neighbor_bead_index != current_bead_global_index)
			if energy[current_bead_index, bead_global_to_local_index(neighbor_bead_index, parameters, polymers_type)] != 0 &&
			   (interaction_list[neighbor_bead_index] == 0 || interaction_list[neighbor_bead_index] == current_bead_global_index)
				interacting_neighbors += 1
				interaction_beads_list[interacting_neighbors] = neighbor_bead_index
			end
		end
	end

	return interacting_neighbors, interaction_beads_list
end

# 计算能量
function calculate_energy(energy_value, energy, current_bead_index, interaction_bead_index, parameters, current_bead_global_index, interaction_bead_global_index)
	energy_value += energy[current_bead_index, interaction_bead_index]
	if bead_index_to_polymer_index(current_bead_global_index, parameters) == bead_index_to_polymer_index(interaction_bead_global_index, parameters)
		energy_value = energy_value + floor(SELF_INTERACTION_ENERGY * exp(real(-(abs(current_bead_index - interaction_bead_index) - 1) / INTERACTION_LENGTH)))
	end
	return energy_value
end

# 计算相对移动位置
function calculate_relative_move_position(connected_bead_displacements, connected_bead_lengths, connection_count)
	relative_move_position = rand(-2:2, 3)  # 初始化为随机移动
	if connection_count > 0
		min_displacement = maximum(connected_bead_displacements .- connected_bead_lengths, dims = 1)
		max_displacement = minimum(connected_bead_displacements .+ connected_bead_lengths, dims = 1)
		diff_displacement = max_displacement .- min_displacement
		relative_move_position = [rand(0:diff_displacement[i]) for i in 1:length(diff_displacement)]' .+ min_displacement
		relative_move_position = vec(relative_move_position)  # 转换为向量
	end
	return relative_move_position
end

# move_pair_beads：遍历并处理与某个珠子连接的所有固定键，计算这些键对应的其他珠子的位置和位移
function process_connections(global_index_1, local_index_1, global_index_2, total_connections, connection_matrix, bead_positions_list, current_bead_position, connection_lengths)
	connected_bead_displacements = []
	connected_bead_lengths = []
	connection_count = 0

	for connected_index in 1:total_connections[local_index_1]
		if global_index_2 - global_index_1 != connection_matrix[local_index_1, connected_index] - local_index_1 # 可逆相互作用较弱
			connection_count += 1
			connected_bead_position = bead_positions_list[global_index_1-local_index_1+connection_matrix[local_index_1, connected_index], :]
			push!(connected_bead_displacements, distance_periodic_boundary(connected_bead_position, current_bead_position))
			push!(connected_bead_lengths, connection_lengths[local_index_1, connected_index])
		end
	end

	return connected_bead_displacements, connected_bead_lengths, connection_count
end

# 滑动移动
function attempt_bead_move(start_position, connection_length, bead_index_lattice, slither_position_list, bead_index, criterion_star, criterion_end, check_range, flag)
	relative_move_position = rand(1:2*connection_length+1, 3) .- connection_length .- 1  # goes from -LL to LL
	new_bead_position = round.(Int, periodic_boundary_3d(start_position .+ relative_move_position))

	occupied_condition = ((bead_index_lattice[new_bead_position[1], new_bead_position[2], new_bead_position[3]] > criterion_star) && (bead_index_lattice[new_bead_position[1], new_bead_position[2], new_bead_position[3]] <= criterion_end))
	if bead_index_lattice[new_bead_position[1], new_bead_position[2], new_bead_position[3]] == 0 || occupied_condition
		rows = [Tuple(slither_position_list[check_range, :][i, :]) for i in 1:size(slither_position_list[check_range, :], 1)] #检查新增珠子是否和已有珠子重叠
		row_tuple = Tuple(new_bead_position)
		has_duplicate_rows = !(row_tuple in rows)
		if has_duplicate_rows
			slither_position_list[bead_index, :] = new_bead_position
		else
			flag = true
		end
	else
		flag = true
	end
	return flag
end

function double_pivot_neighbor_search(bead_index::Int, bead_position::Vector{Int}, neighbor_radius::Int)
	neighbor_list = Int[]  # 初始化空的邻居列表
	search_bead_position = fill(0, 3)  # 用于存储搜索位置

	for dx in -neighbor_radius:neighbor_radius
		for dy in -neighbor_radius:neighbor_radius
			for dz in -neighbor_radius:neighbor_radius

				search_bead_position = round.(Int, periodic_boundary_3d(bead_position .+ [dx, dy, dz]))
				neighbor_bead_index = parameters["bead_index_lattice"][search_bead_position[1], search_bead_position[2], search_bead_position[3]]

				if neighbor_bead_index != 0 && neighbor_bead_index != bead_index
					append!(neighbor_list, neighbor_bead_index)
				end
			end
		end
	end
	return neighbor_list
end

function double_pivot_swap_beads(bead1::Int, bead2::Int, bead_positions_list::Matrix{Int}, interaction_list::Vector{Int}, bead_index_lattice)
	# 临时存储粒子坐标的向量
	temporary_pos1 = Vector{Int}(undef, 3)
	temporary_pos2 = Vector{Int}(undef, 3)

	# 交换粒子的坐标
	temporary_pos1 = bead_positions_list[bead1, :]
	temporary_pos2 = bead_positions_list[bead2, :]
	bead_positions_list[bead1, :] = temporary_pos2
	bead_positions_list[bead2, :] = temporary_pos1

	# 交换绑定珠子
	bead1_interaction, bead2_interaction = interaction_list[bead1], interaction_list[bead2]

	if bead1_interaction != bead2
		interaction_list[bead1], interaction_list[bead2] = bead2_interaction, bead1_interaction
		if bead1_interaction != 0
			interaction_list[bead1_interaction] = bead2
		end
		if bead2_interaction != 0
			interaction_list[bead2_interaction] = bead1
		end
	end

	# 交换晶格上位置
	bead_index_lattice[temporary_pos1[1], temporary_pos1[2], temporary_pos1[3]] = bead2
	bead_index_lattice[temporary_pos2[1], temporary_pos2[2], temporary_pos2[3]] = bead1

	return bead_positions_list, interaction_list, bead_index_lattice
end

# 计算移动前、后的能量和Rosenbluth因子
function calculate_energy_and_Rosenbluth(current_polymer_type, slither_position_list, interaction_list, current_polymer_index, max_beads, energy, parameters, neighbor_matrix, bead_positions_list, polymers_type, bead_index_lattice)
	rotatable_beads_list = parameters["rotatable_beads_list"]

	energy_before_move = 0.0
	energy_after_move = 0.0
	before_log_bind_count = 0.0
	after_log_bind_count = 0.0
	
	slither_rotatable_list = zeros(Int, max_beads) # 存储滑移后链的绑定状态

	current_bead_global_index_reference = (current_polymer_index - 1) * max_beads  # 当前蛋白质链的第一个珠子的全局索引
	current_bead_local_index_reference = sum(BEADS_PER_POLYMER[1:current_polymer_type-1])
	for bead_index in 1:BEADS_PER_POLYMER[current_polymer_type]  # 移动蛋白质链中的每个珠子
		current_bead_index = current_bead_local_index_reference + bead_index
		current_bead_global_index = current_bead_global_index_reference + bead_index

		if rotatable_beads_list[current_bead_index]  # 该珠子能够发生相互作用
			bound_count_before_slither_move = 1
			bound_count_after_slither_move = 1

			# 计算移动前的能量
			if (interaction_list[current_bead_global_index] != 0) &&
			   ((interaction_list[current_bead_global_index] <= current_bead_global_index_reference) ||
				(interaction_list[current_bead_global_index] > current_bead_global_index)) # 珠子已经被绑定，且该绑定要不在索引低的其他链上，要不在索引更高的珠子上

				interaction_bead_global_index = interaction_list[current_bead_global_index]
				interaction_bead_index = bead_global_to_local_index(interaction_bead_global_index, parameters, polymers_type)
				energy_before_move = calculate_energy(energy_before_move, energy, current_bead_index, interaction_bead_index, parameters, current_bead_global_index, interaction_bead_global_index)
			end

			# 计算移动前的Rosenbluth因子
			if (interaction_list[current_bead_global_index] <= current_bead_global_index_reference) || (interaction_list[current_bead_global_index] > current_bead_global_index)
				for neighbor_index in 1:26
					neighbor_position = round.(Int, periodic_boundary_3d(bead_positions_list[current_bead_global_index, :] .+ neighbor_matrix[neighbor_index, :])) # 周围格子坐标
					neighbor_bead_index = bead_index_lattice[neighbor_position[1], neighbor_position[2], neighbor_position[3]]

					if neighbor_bead_index != 0 # 存在邻居

						has_interaction_energy = energy[current_bead_index, bead_global_to_local_index(neighbor_bead_index, parameters, polymers_type)] != 0
						if has_interaction_energy # 该珠子和邻居间存在相互作用能量

							is_different_chain = (neighbor_bead_index <= current_bead_global_index_reference) ||
												 (neighbor_bead_index > current_bead_global_index_reference + BEADS_PER_POLYMER[current_polymer_type])
							if is_different_chain # 邻居在不同链上

								is_unbound_or_bound_in_current_chain =
									(interaction_list[neighbor_bead_index] == 0) ||
									((interaction_list[neighbor_bead_index] >= current_bead_global_index) &&
									 (interaction_list[neighbor_bead_index] <= current_bead_global_index_reference + BEADS_PER_POLYMER[current_polymer_type]))
								if is_unbound_or_bound_in_current_chain # 邻居未绑定或者邻居绑定到了当前链（当前链的所有可逆键要断开，因此该珠子后续可以绑定）
									bound_count_before_slither_move += 1
								end

							else # 邻居在同一条链上
								is_later_in_chain = neighbor_bead_index > current_bead_global_index # 邻居的索引更高，避免重复计算
								is_not_within_current_chain =
									(interaction_list[neighbor_bead_index] == 0) ||
									(interaction_list[neighbor_bead_index] <= current_bead_global_index_reference) ||
									(interaction_list[neighbor_bead_index] >= current_bead_global_index) # 邻居未绑定，或者绑定到了其他链，或者绑定到索引更高的珠子上（当前链的所有可逆键要断开，因此该珠子后续可以绑定）
								if is_later_in_chain && is_not_within_current_chain
									bound_count_before_slither_move += 1
								end
							end
						end
					end
				end
			end

            # 计算移动后的Rosenbluth因子，为当前珠子寻找一个新的绑定对象，slither_rotatable_list为移动后的绑定列表
			if slither_rotatable_list[bead_index] == 0 
				current_position = slither_position_list[bead_index, :] # 移动后珠子的坐标
				interaction_candidates = zeros(Int, 27) 
                
				# 不同链上的绑定对象
				for neighbor_index in 1:26 # 检查邻居珠子
					neighbor_position = round.(Int, periodic_boundary_3d(current_position .+ neighbor_matrix[neighbor_index, :]))
					neighbor_bead_index = bead_index_lattice[neighbor_position[1], neighbor_position[2], neighbor_position[3]]

					is_different_chain = (neighbor_bead_index <= current_bead_global_index_reference) ||
										 (neighbor_bead_index > current_bead_global_index_reference + BEADS_PER_POLYMER[current_polymer_type]) # 邻居和当前珠子不在同一条链上

					if neighbor_bead_index != 0 && is_different_chain # 有邻居且不在同一条链上
						has_interaction_energy = energy[current_bead_index, bead_global_to_local_index(neighbor_bead_index, parameters, polymers_type)] != 0
						is_unbound_or_bound_in_current_chain =
							(interaction_list[neighbor_bead_index] == 0) ||
							((interaction_list[neighbor_bead_index] > current_bead_global_index_reference) &&
							 (interaction_list[neighbor_bead_index] <= current_bead_global_index_reference + BEADS_PER_POLYMER[current_polymer_type])) # 邻居未绑定或者绑定到了当前链（当前链的所有可逆键要断开，因此该珠子后续可以绑定）

						if has_interaction_energy && is_unbound_or_bound_in_current_chain
							if all(slither_rotatable_list[1:bead_index-1] .!= neighbor_bead_index) # 该邻居没有被前面的珠子绑定
								bound_count_after_slither_move += 1
								interaction_candidates[bound_count_after_slither_move] = neighbor_bead_index # interaction_candidates的第一个位置为不绑定状态
							end
						end
					end
				end
                # 同一条链上的绑定对象
				for future_bead_index in bead_index+1:BEADS_PER_POLYMER[current_polymer_type] 
					if energy[current_bead_index, current_bead_local_index_reference+future_bead_index] != 0
						future_position = slither_position_list[future_bead_index, :]
						if (slither_rotatable_list[future_bead_index] == 0) && all(abs.(distance_periodic_boundary(current_position, future_position)) .<= 1)
							bound_count_after_slither_move += 1
							interaction_candidates[bound_count_after_slither_move] = current_bead_global_index_reference + future_bead_index 
						end
					end
				end

				if bound_count_after_slither_move > 0
					selected_interaction = interaction_candidates[rand(1:bound_count_after_slither_move)] # 选择该珠子的绑定
					slither_rotatable_list[bead_index] = selected_interaction

					if selected_interaction != 0
						if (selected_interaction > current_bead_global_index) &&
						   (selected_interaction <= current_bead_global_index_reference + BEADS_PER_POLYMER[current_polymer_type]) # 绑定的珠子在同一条链上，且索引更高
							slither_rotatable_list[selected_interaction-current_bead_global_index_reference] = current_bead_global_index
						elseif (selected_interaction > current_bead_global_index_reference) &&
							   (selected_interaction < current_bead_global_index) # 绑定的珠子在同一条链上，索引低则报错
							println("Slither is binding to something it should not be able to!")
							return
						elseif selected_interaction == current_bead_global_index
							println("Slither is binding to itself!")
							return
						end

						interaction_bead_index = bead_global_to_local_index(selected_interaction, parameters, polymers_type)
						energy_after_move = calculate_energy(energy_after_move, energy, current_bead_index, interaction_bead_index, parameters, current_bead_global_index, selected_interaction)
					end
				end
			end

			before_log_bind_count += log(bound_count_before_slither_move)
			after_log_bind_count += log(bound_count_after_slither_move)
		end
	end
	return before_log_bind_count, after_log_bind_count, energy_before_move, energy_after_move, slither_rotatable_list
end
