# STRUCTURES USED IN THE PROBLEM

# Input data
#-----------------------------------------------

# Input parameters 
@with_kw struct InputParam{F<:Float64,I<:Int}
    NYears::F                                     # Number of years
    NMonths::I
    NStages::I                                    # Number of stages of N months in the problem FORMULATION-- calcolato come NYears/NMonths*12
    NHoursStep::F                                 # Number of hours in each time step 
    Hours_rolling::I
    Hours_saved::I
    Big::F                                        # A big number
    conv::F                                       # A small number for degradation convergence
    bin::I                                       # Number of binaries for discretization of SE
end

# Battery's characteristics
@with_kw struct BatteryParam{F<:Float64,I<:Int}
    min_SOC::F                                     # Battery's minimum energy storage capacity
    max_SOC::F                                     # Batter's maximum energy storage capacity
    Eff_charge::F                                  # Battery's efficiency for charging operation
    Eff_discharge::F                               # Battery's efficiency for discharging operation
    min_P::F
    max_P::F
    max_SOH::F                                     # Maximum SOH that can be achieved because of volume issues
    min_SOH::F                                     # Minimum SOH to be respected by contract
    Nfull::I                                       # Maximum number of full cycles for DoD=100%     
    fix::F                                         # fixed costs for battery replacement                                   
    cost::F
end
  
# solver parameters
@with_kw struct SolverParam{F<:Float64,I<:Int}
    MIPGap::F 
    MIPFocus::I
    Method::F
    Cuts::F
    Heuristics::F
end
  
# Indirizzi cartelle
@with_kw struct caseData{S<:String}
    DataPath::S
    InputPath::S
    ResultPath::S
    CaseName::S
end

# runMode Parameters
@with_kw mutable struct runModeParam{B<:Bool}

    # Solver settings
    solveMIP::B     #If using SOS2

    batterySystemFromFile::B 

    #runMode self defined reading of input 
    setInputParameters::B             #from .in file
 
    excel_savings::B 
    endlife::B
    endlife_rolling::B
    rolling_newDeg::B

end

# Optimization problem
struct BuildStageProblem_3
    M::Any
    soc::Any
    soc_quad::Any
    e_charge::Any 
    e_discharge::Any
    #bin_op::Any
    deg::Any
    x::Any
    y::Any
    z::Any
    h_x::Any
    h_y::Any
    h_z::Any
    w_xx::Any
    w_yy::Any
    w_zz::Any
    w_xy::Any
    w_xz::Any
    w_zy::Any
    capacity::Any
    revamping::Any
    e::Any
    rev_vendita::Any
    rev_acquisto::Any
    h_xx::Any
    h_xy::Any
    h_xz::Any
    h_yy::Any
    h_zz::Any
    h_yz::Any
    aux_deg::Any
end

struct Results_3
    objective::Any
    gain_stage::Any
    cost_rev::Any
    deg_stage::Any
    soc::Any
    e_charge::Any
    e_discharge::Any
    #bin_op::Any
    deg::Any
    soc_quad::Any
    x::Any
    y::Any
    z::Any
    w_xx::Any
    w_yy::Any
    w_zz::Any
    w_xy::Any
    w_xz::Any
    w_zy::Any
    h_x::Any
    h_y::Any
    h_z::Any
    rev::Any
    cap::Any
    e::Any
    rev_vendita::Any
    rev_acquisto::Any
    h_xx::Any
    h_xy::Any
    h_xz::Any
    h_yy::Any
    h_zz::Any
    h_yz::Any
    aux_deg::Any
end

struct BuildStageProblem_4
    M::Any
    soc::Any
    soc_quad::Any
    charge::Any 
    discharge::Any
    deg::Any
    x::Any
    y::Any
    z::Any
    u::Any
    w_xx::Any
    w_yy::Any
    w_zz::Any
    w_xy::Any
    w_xz::Any
    w_zy::Any
    w_uu::Any
    w_xu::Any
    w_yu::Any
    w_zu::Any
    capacity::Any
    revamping::Any
    e::Any
    rev_vendita::Any
    rev_acquisto::Any
end

struct Results_4
    objective::Any
    gain_stage::Any
    cost_rev::Any
    deg_stage::Any
    soc::Any
    charge::Any
    discharge::Any
    #bin::Any
    deg::Any
    soc_quad::Any
    x::Any
    y::Any
    z::Any
    u::Any
    w_xx::Any
    w_yy::Any
    w_zz::Any
    w_uu::Any
    w_xy::Any
    w_xz::Any
    w_zy::Any
    w_xu::Any
    w_yu::Any
    w_zu::Any
    rev::Any
    cap::Any
    e::Any
    rev_vendita::Any
    rev_acquisto::Any
end
