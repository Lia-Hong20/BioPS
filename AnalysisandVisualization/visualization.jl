using Plots
using XLSX
using LinearAlgebra
using Statistics

# Set PyPlot backend to support .eps output
pyplot()
default(fontfamily="Arial", titlefontsize=16, guidefontsize=14, tickfontsize=12, legendfontsize=12, grid=false)

# Function to read simulation data
function read_cluster_data(file_path)
    try
        xlsx_data = XLSX.readxlsx(file_path)
        data_values = XLSX.getdata(xlsx_data["Sheet1"])
        chain_data = data_values[:, 1]  # 总聚合物链
        cluster_data = data_values[:, 2]  # 最大簇
        time_data = data_values[:, 3]   # 时间
        energy_data = data_values[:, 4]  # 总能量
        return time_data, chain_data, cluster_data, energy_data
    catch e
        @error "读取数据失败：$e"
        return nothing, nothing, nothing, nothing
    end
end

# Function to read bead position data
function read_bead_positions(file_path)
    try
        xlsx_data = XLSX.readxlsx(file_path)
        bead_positions = XLSX.getdata(xlsx_data["Bead Positions"])
        valid_rows = findall(row -> !any(ismissing, row), eachrow(bead_positions))
        return bead_positions[valid_rows, :]
    catch e
        @error "Failed to read bead position data: $e"
        return nothing
    end
end

# Files to read
file_paths = ["simulation_results.xlsx"]
bins = 500
num_steps = 27  # Only plot the most recent num_steps data points

# Define plots
p1 = plot(xlabel="Time", ylabel="Largest Cluster Size", title="Cluster Size Over Time", legend=:topright)
p2 = plot(xlabel="Time", ylabel="Number of Polymer Chains", title="Polymer Chains Over Time", legend=:topright)
p3 = plot(xlabel="Time", ylabel="Total Energy", title="Energy Over Time", legend=:topright)

# 3D scatter plots for bead positions
p4_list = [plot(xlabel="X", ylabel="Y", zlabel="Z", title="Bead Positions $(i)", legend=false) for i in 1:length(file_paths)]

colors = palette(:tab10)  

for (i, file_path) in enumerate(file_paths)
    time, chain, cluster, energy = read_cluster_data(file_path)
    if time === nothing
        continue
    end

    # 3D scatter plots for bead positions
    idx_range = max(1, length(time) - num_steps + 1):length(time)

    plot!(p1, time, cluster, label="File $(i)", color=colors[i], lw=1)
    plot!(p2, time, chain, label="File $(i)", color=colors[i], lw=2, xlim = (0, time[end]))
    plot!(p3, time, energy, label="File $(i)", color=colors[i], lw=2, xlim = (0, time[end]))

    # Plot bead positions if available
    bead_positions = read_bead_positions(file_path)
    if bead_positions !== nothing
        scatter!(p4_list[i], bead_positions[:, 1], bead_positions[:, 2], bead_positions[:, 3], markersize=2, color=colors[i])
    end
end

# Generate final figure layout and save as .eps
plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 800))
savefig("simulation_results.eps")
