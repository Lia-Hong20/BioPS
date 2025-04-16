using Plots
using XLSX

include("constants.jl")
include("parameters.jl")

include("initialization.jl")
include("data_initialization.jl")
include("system_initialization.jl")

include("utility_functions.jl")
include("simulation.jl")
include("biochemical_reaction.jl")
include("phase_separation.jl")
include("move_behavior.jl")

include("analysis.jl")
include("network_analysis.jl")
include("results_writer.jl")

# 主函数
function main(parameters)
    # 调用 ResultsWriter 模块中的函数
    init_results_file()  # 使用模块名调用函数
    parameters = initialize_parameters(parameters)

    # 运行模拟
    @time begin
        largest_cluster, parameters = run_simulation(parameters, STEPS_PER_SIMULATION, TOTAL_SIMULATION_STEPS)
    end
end

main(parameters)