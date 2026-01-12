#RUN FILE

# Calls the Packages used for the optimization problem
using JuMP
using Printf
using Gurobi
#using CPLEX
using MathOptInterface
using JLD
using TimerOutputs
using DataFrames
using XLSX
using Parameters
using Dates

using CSV

# Calls the other Julia files
include("Structures.jl")
include("SetInputParameters.jl")
include("solveOptimizationAlgorithm.jl")        #solveOptimizationAlgorithm_3cuts
include("ProblemFormulationInequalities.jl")      #ProblemFormulationCutsTaylor_3
include("Saving in xlsx.jl")

date = string(today())

# PREPARE INPUT DATA
to = TimerOutput()

@timeit to "Set input data" begin

  #Set run case - indirizzi delle cartelle di input ed output
  case = set_runCase()

  @unpack (DataPath,InputPath,ResultPath,CaseName) = case;

  # Set run mode (how and what to run) and Input parameters
  runMode = read_runMode_file()
  InputParameters = set_parameters(runMode, case)
  @unpack (NYears, NMonths, NHoursStep, NStages, Big, conv, bin, Hours_rolling, Hours_saved)= InputParameters;    #NSteps, NHoursStage

  # Upload battery's characteristics
  Battery = set_battery_system(runMode, case)
  @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull, cost, fix) = Battery; 

  # Set solver parameters (Gurobi etc)
  SolverParameters = set_solverParameters()

  a,b,c,disc = calculate_coefficients(min_SOC,max_SOC,bin)

  # Read power prices from a file [€/MWh]
  #Steps_stages = [0 1132 2259 3397 4465 5650 6723 7806 8843 9953 10946 12124 13171 14283 15273 16385 17375 18460 19470 20496 21505]
  Steps_stages = [0 2259 4465 3723 8843 10946 13171 15273 17375 19470 21505] # yearly revamping
  Steps_stop = [128 127 120 111 127 115 115 123 120]                      # 3 weeks dowtime 
  NSteps = Steps_stages[NStages+1]

  Battery_price_purchase = read_csv("Mid-cost projections 2026-2036.csv",case.DataPath)
  Battery_price_sale = set_price(Battery_price_purchase,cost);
  
  Power_prices = read_csv("Prezzi_2026_2035_filtered.csv", case.DataPath);    

  # Where and how to save the results
  FinalResPath= set_run_name(case, ResultPath, NSteps)

end


@timeit to "Solve optimization problem" begin
  if bin ==3
  ResultsOpt = solveOptimizationProblem_3(InputParameters,SolverParameters,Battery);
  else
    ResultsOpt = solveOptimizationProblem_4(InputParameters,SolverParameters,Battery);
  end
end

# SAVE DATA IN EXCEL FILES
if runMode.excel_savings
  cartella = "C:\\GitHub_Asja\\OCRA_unitary\\Results"
  cd(cartella)
  Saving = data_saving(InputParameters,ResultsOpt)
else
  println("Solved without saving results in xlsx format.")
end

