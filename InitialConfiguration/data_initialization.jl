using Random
using DelimitedFiles
using CSV
using DataFrames

function read_energy_matrix(file_path::String, beads_count::Int)
	energy = zeros(Float64, beads_count, beads_count)
	open(file_path, "r") do file
		for (i, line) in enumerate(eachline(file))
			values = parse.(Float64, split(line))
			@assert length(values) == beads_count "Invalid energy matrix dimension"
			energy[i, :] = values
		end
	end
	validate_symmetry(energy)
	return energy
end

function read_connection_matrix(file_path::String, beads_count::Int)
    # 初始化连接计数器
    total_connections = zeros(Int, beads_count)  # 每个节点的所有连接数
    connections_with_higher_indices = zeros(Int, beads_count)  # 只计数 i < j 的连接数

    # 首次遍历文件，计算连接数
    open(file_path, "r") do file
        for line in eachline(file)
            i, j, len = parse.(Int, split(line, ","))
            if i >= j
                error("Error in Connection File, index i must be smaller than j")
            end
            connections_with_higher_indices[i] += 1
            total_connections[i] += 1
            total_connections[j] += 1
        end
    end

    # 初始化连接矩阵和长度矩阵
    max_connections = maximum(total_connections)
    connection_matrix = zeros(Int, beads_count, max_connections)
    connection_lengths = zeros(Int, beads_count, max_connections)
    connections_done = zeros(Int, beads_count)  # 重置连接计数器

    # 重新遍历文件，填充连接矩阵
    open(file_path, "r") do file
        for line in eachline(file)
            i, j, len = parse.(Int, split(line, ","))
            if i >= j
                error("Error in Connection File, index i must be smaller than j")
            end
            connections_done[i] += 1
            connections_done[j] += 1
            connection_matrix[i, connections_done[i]] = j
            connection_lengths[i, connections_done[i]] = len
            connection_matrix[j, connections_done[j]] = i
            connection_lengths[j, connections_done[j]] = len
        end
    end

    return connection_matrix, connection_lengths, total_connections
end

function validate_symmetry(matrix::Array{Float64, 2})
	for i in 1:size(matrix, 1)
		for j in 1:size(matrix, 2)
			@assert matrix[i, j] == matrix[j, i] "Matrix is not symmetric"
		end
	end
end

function read_input_data(parameters::Dict{String, Any})
	parameters["droplet_polymers_count"] = START_FROM_DROPLET ? parameters["droplet_polymers_count"] : 0
	parameters["slab_polymers_count"] = START_FROM_SLAB ? parameters["slab_polymers_count"] : 0
	parameters["total_polymers_count"] = sum(parameters["num_polymers_per_type"])

	# 能量矩阵
	beads_count = sum(BEADS_PER_POLYMER)
	parameters["energy"] = read_energy_matrix("interaction_matrix.txt", beads_count)

	# 连接矩阵
	parameters["connection_matrix"], parameters["connection_lengths"], parameters["total_connections"] = read_connection_matrix("polymer_topology.txt", beads_count)

	return parameters
end

function basic_setup(parameters)
	# 检查是否存在旋转移动的能量状态
	if !any(parameters["energy"] .!= 0) && parameters["reconnect_move_count"] > 0
		println("Rotation moves are not allowed without energy states.")
		parameters["reconnect_move_count"] = 0
	end

	# 初始化 reptation_move 移动
	if parameters["reptation_move_count"] > 0
		is_slither_move_few_beads = true
		valid_slither_few_indices = filter(i -> begin
				if BEADS_PER_POLYMER[i] > 1
					first, last = sum(BEADS_PER_POLYMER[1:i-1]) + 1, sum(BEADS_PER_POLYMER[1:i])
					con = parameters["total_connections"]
					lengths = parameters["connection_lengths"]
					return con[first] == 1 && con[last] == 1 &&
						   all(con[first+1:last-1] .== 2) &&
						   all(lengths[first+1:last-1, 1] .== lengths[first+1:last-1, 2])
				end
				return false
			end, 1:NUM_POLYMER_TYPES)
		parameters["reptation_move_indices"] = valid_slither_few_indices
		
		if length(valid_slither_few_indices) == 0
			is_slither_move_few_beads = false
			println("Slither_move_few_beads aren't appropriate for this system, turning them off")
			parameters["reptation_move_count"] = 0
		end
	else
		is_slither_move_few_beads = false
	end

	# 检查初始条件设置
	if START_FROM_DROPLET
		if sum(parameters["droplet_polymers_count"] .* BEADS_PER_POLYMER) > 4 * PI * DROPLET_RADIUS^3/3
			error("Not enough room for all the proteins you want in the dense phase!")
		end
		for i in 1:NUM_POLYMER_TYPES
			parameters["droplet_polymers_count"][i] = min(parameters["droplet_polymers_count"][i], parameters["num_polymers_per_type"][i])
		end
	end

	if sum(parameters["num_polymers_per_type"]) > parameters["total_polymers_count"]
		error("Not enough room to fit all the proteins you've asked for in MaxSumNk!")
	end

	if START_FROM_SLAB
		if START_FROM_IMPORTED
			error("Starting from imported is not compatible with starting from slab.")
		end
		if LATTICE_SIZE[3] < 2 * LATTICE_SIZE[1] || LATTICE_SIZE[3] < 2 * LATTICE_SIZE[2]
			error("Slab direction size is not long enough; make it longer.")
		end
		if sum(parameters["slab_polymers_count"] .* BEADS_PER_POLYMER) > 0.7 * LATTICE_SIZE[1] * LATTICE_SIZE[2] * SLAB_HEIGHT
			error("Not enough room for all the proteins you want in the slab dense phase!")
		end
		for i in 1:NUM_POLYMER_TYPES
			parameters["slab_polymers_count"][i] = min(parameters["slab_polymers_count"][i], parameters["num_polymers_per_type"][i])
		end
	elseif !(LATTICE_SIZE[1] == LATTICE_SIZE[2] == LATTICE_SIZE[3])
		println("Your box is not a cube, assuming you are starting Slab.")
	end

	parameters["max_beads"] = maximum(BEADS_PER_POLYMER)

	# 邻居列表
	neighbors = collect(Iterators.product(-1:1, -1:1, -1:1))
	neighbors = filter(x -> x != (0, 0, 0), neighbors)
	neighbor_matrix = zeros(Int, length(neighbors), 3)
	for (i, neighbor) in enumerate(neighbors)
		neighbor_matrix[i, :] = collect(neighbor)
	end
	parameters["neighbor_matrix"] = neighbor_matrix

	# 可旋转的珠子
	parameters["rotatable_beads_list"] = any(parameters["energy"] .!= 0, dims = 2)

	return parameters
end


