function move_process_simulation(parameters)
	# 定义移动类型及其对应的函数
	move_functions = [
		single_bead_move,
		paired_bead_move,
		reconnect_move,
		reptation_move,
		chain_translation_move,
		chain_formed_cluster_move,
		non_largest_cluster_move,
		double_pivot_move,
	]

	# 根据 move_choice 选择对应的移动函数
	move_choice = rand(1:parameters["moves_counts"][8]) - 1
	for (i, move_count) in enumerate(parameters["moves_counts"])
		if move_choice < move_count
			parameters = move_functions[i](parameters)
			break
		end
	end

	return parameters
end