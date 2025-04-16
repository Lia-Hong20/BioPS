function identify_protein_clusters(parameters) #网络分析，构建蛋白质的网络：遍历所有蛋白质链，遍历当前蛋白质链上的所有珠子的相互作用，得到所有存在相互作用的蛋白质的集合
    total_polymers_count = parameters["total_polymers_count"]
    polymers_type = parameters["polymers_type"]
    interaction_list = parameters["interaction_list"]
    num_polymers_per_type = parameters["num_polymers_per_type"]

    cluster_ranges = zeros(Int, total_polymers_count, 2) # 每个群集的起始和结束位置
    clustered_polymers_list = zeros(Int,total_polymers_count) # 存储所有群集中的蛋白质索引
    current_cluster = Vector{Int}() # 临时存储当前处理的群集中的蛋白质
    is_unclustered = fill(true, total_polymers_count) # 标记每个蛋白质是否未被处理
    total_clusters = 0 # 蛋白质网络中的群集数量
    first_cluster = true

    total_polymer_count = sum(num_polymers_per_type)
    for polymer_index in 1:total_polymer_count # 遍历所有蛋白质链
 
        if is_unclustered[polymer_index] # 当前蛋白质未被处理，初始化一个新的集群
            total_clusters += 1 # 群集的数量计数

            empty!(current_cluster) # 储存新集群中的链的索引
            push!(current_cluster, polymer_index)  # 将当前蛋白质添加到群集中
            
            is_unclustered[polymer_index] = false # 将当前蛋白质标记为已处理

            current_cluster_index = 1
            while current_cluster_index <= length(current_cluster) # 循环直到所有与群集中的蛋白质有相互作用的蛋白质都被添加到群集中
                current_polymer = current_cluster[current_cluster_index]
                current_polymer_type = polymers_type[current_polymer] # 与当前蛋白质相连的蛋白质类型
                bead_start_index = (current_polymer-1) * BEADS_PER_POLYMER[current_polymer_type] + 1
                bead_end_index = bead_start_index + BEADS_PER_POLYMER[current_polymer_type] - 1

                for bead_index in bead_start_index:bead_end_index # 遍历当前蛋白质链上的所有珠子
                    interacting_index = interaction_list[bead_index]
                    
                    if interacting_index != 0 # 珠子与其他蛋白质相连
                        interacting_polymer = bead_index_to_polymer_index(interacting_index, parameters)
                        
                        if is_unclustered[interacting_polymer] # 相互作用的蛋白质链未被处理
                            push!(current_cluster, interacting_polymer)
                            is_unclustered[interacting_polymer] = false
                        end
                    end
                end
                current_cluster_index += 1
            end

            if first_cluster
                cluster_ranges[total_clusters, :] = [1, length(current_cluster)] # 第一个群集的起始位置和结束位置
                first_cluster = false
            else
                cluster_ranges[total_clusters, :] = [cluster_ranges[total_clusters-1, 2] + 1, cluster_ranges[total_clusters-1, 2] + length(current_cluster)] # 记录新群集的起始和结束位置
            end
            clustered_polymers_list[cluster_ranges[total_clusters, 1]:cluster_ranges[total_clusters, 2]] .= current_cluster[1:length(current_cluster)]
        end
    end

    return total_clusters, cluster_ranges, clustered_polymers_list
end

