using XLSX  # 在模块内部加载 XLSX 包

function init_results_file()
    XLSX.openxlsx("simulation_results.xlsx", mode="w") do xf
        sheet1 = xf[1]  # 第一个工作表
        sheet1[1, 1:5] = ["Protein count", "Largest Cluster", "Time", "Total Energy", "Move Steps"]

        sheet2 = XLSX.addsheet!(xf, "Bead Positions")  # 添加第二个工作表
    
        sheet3 = XLSX.addsheet!(xf, "Interaction List")  # 添加第三个工作表
    end
end

function write_simulation_results(step, STEPS_PER_SIMULATION, parameters, largest_cluster, time, move_step)
    XLSX.openxlsx("simulation_results.xlsx", mode="rw") do xf
        # 处理第一个工作表
        sheet1 = xf[1]  # 第一个工作表
        row = step ÷ STEPS_PER_SIMULATION + 1  # 计算当前写入的行数
        sheet1[row, 1] = parameters["num_polymers_per_type"][1]
        sheet1[row, 2] = largest_cluster
        sheet1[row, 3] = time
        sheet1[row, 4] = parameters["total_energy"]
        sheet1[row, 5] = move_step

        # 清空并更新第二个工作表（只清空所有行）
        sheet2 = xf["Bead Positions"]
        rows2 = length(sheet2[:, 1])  # 获取工作表的行数，通过第一列来估算
        cols2 = length(sheet2[1, :])  # 获取工作表的列数，通过第一行来估算
        for i in 1:rows2
            for j in 1:cols2
                sheet2[i, j] = nothing  # 清空每一行的每一列
            end
        end
        # 填充新数据
        for i in 1:size(parameters["bead_positions_list"], 1)
            sheet2[i, 1:3] = parameters["bead_positions_list"][i, :]
        end

        # 清空并更新第三个工作表（只清空所有行）
        sheet3 = xf["Interaction List"]
        rows3 = length(sheet3[:, 1])  # 获取工作表的行数，通过第一列来估算
        for i in 1:rows3
            sheet3[i, 1] = nothing  # 清空每一行的第1列
        end
        # 填充新数据
        for i in 1:length(parameters["interaction_list"])
            sheet3[i, 1] = parameters["interaction_list"][i, 1]
        end
    end
end