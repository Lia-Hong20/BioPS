function find_largest_cluster(parameters)
    total_clusters, cluster_ranges, clustered_polymers_list = identify_protein_clusters(parameters)
    cluster_sizes = cluster_ranges[1:total_clusters, 2].- cluster_ranges[1:total_clusters, 1].+ 1
    largest_cluster_size = maximum(cluster_sizes)
    largest_cluster_index = findfirst(isequal(largest_cluster_size), cluster_sizes)
    largest_cluster_range = cluster_ranges[largest_cluster_index, :]
    return largest_cluster_size, largest_cluster_range, clustered_polymers_list
end

# 计算配对分布函数 P(2)(r)
function pair_distribution_function(bead_positions_list::Matrix{Int}, bins::Int, N::Int)
    # 计算最大距离和 bin 宽度
    max_distance = sqrt(3) * LATTICE_SIZE[1] / 2.0  # 在周期性立方晶格中，两个单体的最大可能距离
    bin_width = max_distance / bins  # 每个 bin 的宽度

    # 初始化配对分布函数数组
    pair_distribution = zeros(bins)  # 初始化配对分布函数数组
    bin_volumes = [(4.0/3.0) * π * (((bin+1) * bin_width)^3 - (bin * bin_width)^3) for bin in 0:(bins-1)]  # 每个bin的体积

    # 全局平均密度
    rho0 = N / (LATTICE_SIZE[1]^3)

    # 计算所有珠子对的距离，并统计到对应的 bin 中
    for i in 1:N-1
        for j in i+1:N
            delta_dis_vec = distance_periodic_boundary(bead_positions_list[i, :], bead_positions_list[j, :])  # 计算周期边界下的距离
            delta_dis = sqrt(sum(delta_dis_vec.^2))
            bin_index = Int(floor(delta_dis / bin_width)) + 1  # 确定距离所在的bin
            if bin_index <= bins
                pair_distribution[bin_index] += 2  # 每对粒子贡献两次（i-j 和 j-i）
            end
        end
    end

    # 计算径向分布函数
    N_pairs = N * (N - 1)  # 粒子对的总数（每对粒子贡献两次）
    radial_distribution = pair_distribution ./ (bin_volumes .* rho0 .* N_pairs)

    return pair_distribution, radial_distribution, bin_width
end
