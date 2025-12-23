# STAGE MAXIMIZATION PROBLEM FORMULATION

function BuildStageProblem_3(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam)       

    @unpack (MIPGap, MIPFocus, Method, Cuts, Heuristics) = SolverParameters;
  
    @unpack (NYears, NMonths, NStages, Big, NHoursStep, conv, bin) = InputParameters;     #NSteps,NHoursStage
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull, fix) = Battery ;         

    k_deg = 1/(2*Nfull)
    Small = 1
    disc = 7
    Beta = (max_SOC-min_SOC)/disc

    M = Model(Gurobi.Optimizer)
    set_optimizer_attribute(M, "MIPGap", 0.01)
    set_optimizer_attribute(M, "MIPFocus", 2)
    set_optimizer_attribute(M, "NorelHeurTime", 15)
    set_optimizer_attribute(M, "Cuts", 2)

    # DEFINE VARIABLES

    @variable(M, min_SOC/max_SOC <= soc[iStep=1:NSteps+1] <= 1, base_name = "Energy")                # MWh   energy_Capacity NSteps
    @variable(M, (min_SOC/max_SOC)^2 <= soc_quad[iStep=1:NSteps+1] <= 1, base_name = "Square energy")

    @variable(M, 0 <= e_charge[iStep=1:NSteps] <= max_SOH/min_SOH, base_name= "Charge")      #max_disc   0<=discharge<=1
    @variable(M, 0 <= e_discharge[iStep=1:NSteps] <= max_SOH/min_SOH, base_name= "Discharge")
    @variable(M, bin_op[iStep=1:NSteps], Bin, base_name = "Binary_operation")
    
    @variable(M, 0 <= deg[iStep=1:NSteps] <= Small*max_SOH/min_SOH, base_name = "Degradation")

    @variable(M, 0 <= revamping[iStage=1:NStages] <= (max_SOH-min_SOH)/max_SOC, base_name = "Revamping")
    @variable(M, min_SOH/min_SOH <= capacity[iStep=1:NSteps+1] <= max_SOH/min_SOH, base_name = "Energy_Capacity")        #energy_Capacity     [iStage=1:NStages]
    @variable(M, e[iStage=1:NStages], Bin, base_name ="Binary Revamp")max_SOH/min_SOH

    @variable(M, 0<= rev_vendita[iStage=1:NStages] <= max_SOH/min_SOH, base_name = "Vendita rev")
    @variable(M, -max_SOH/min_SOH <= rev_acquisto[iStage=1:NStages] <= 0, base_name = "Acquisto rev")

    #VARIABLES FOR DISCRETIZATION of Stored Energy

    @variable(M, x[iStep=1:NSteps+1], Bin, base_name = "Binary_1")
    @variable(M, y[iStep=1:NSteps+1], Bin, base_name = "Binary_2")
    @variable(M, z[iStep=1:NSteps+1], Bin, base_name = "Binary_3")

    @variable(M, w_xx[iStep=1:NSteps+1] , Bin, base_name = "xx")
    @variable(M, w_yy[iStep=1:NSteps+1], Bin, base_name = "yy")
    @variable(M, w_zz[iStep=1:NSteps+1], Bin, base_name = "zz")
    @variable(M, w_xy[iStep=1:NSteps+1], Bin, base_name = "xy")
    @variable(M, w_xz[iStep=1:NSteps+1], Bin, base_name = "xz")
    @variable(M, w_zy[iStep=1:NSteps+1], Bin, base_name = "yz")
    
    @variable(M, 0 <= h_x[iStep=1:NSteps+1] <= max_SOH/min_SOH, base_name = "Aux_1")
    @variable(M, 0 <= h_y[iStep=1:NSteps+1] <= max_SOH/min_SOH, base_name = "Aux_2")
    @variable(M, 0 <= h_z[iStep=1:NSteps+1] <= max_SOH/min_SOH, base_name = "Aux_3")

    @variable(M, 0 <= h_xx[iStep=1:NSteps+1] <= max_SOH/min_SOH, base_name = "Aux_xx")
    @variable(M, 0 <= h_xy[iStep=1:NSteps+1] <= max_SOH/min_SOH, base_name = "Aux_xy")
    @variable(M, 0 <= h_xz[iStep=1:NSteps+1] <= max_SOH/min_SOH, base_name = "Aux_xz")
    @variable(M, 0 <= h_yy[iStep=1:NSteps+1] <= max_SOH/min_SOH, base_name = "Aux_yy")
    @variable(M, 0 <= h_zz[iStep=1:NSteps+1] <= max_SOH/min_SOH, base_name = "Aux_zz")
    @variable(M, 0 <= h_yz[iStep=1:NSteps+1] <= max_SOH/min_SOH, base_name = "Aux_yz")


    # DEFINE OBJECTIVE function - length(Battery_price) = NStages+1=21
    @objective(
      M,
      MathOptInterface.MAX_SENSE, 
      sum(Power_prices[iStep]*(e_discharge[iStep]*Eff_discharge-e_charge[iStep]/Eff_charge) for iStep=1:NSteps) -               # 
      Battery_price_purchase[1]*(revamping[1]) 
      #-sum(Battery_price_purchase[iStage]*(revamping[iStage]) for iStage=1:NStages) +
      - sum(Battery_price_purchase[iStage]*(capacity[Steps_stages[iStage]+2] + rev_acquisto[iStage]) for iStage=2:NStages) 
      + sum(Battery_price_sale[iStage]*(capacity[Steps_stages[iStage]+1] - rev_vendita[iStage]) for iStage=2:NStages) +
      Battery_price_sale[NStages+1]*(capacity[end]- min_SOH/min_SOH)  
      -sum(fix*e[iStage] for iStage=1:NStages) 
      + 2300    
      )
         
    # DEFINE CONSTRAINTS
    @constraint(M,energy[iStep=1:NSteps], e_charge[iStep]-e_discharge[iStep] == Beta*((h_x[iStep+1]-h_x[iStep])+2*(h_y[iStep+1]-h_y[iStep])+4*(h_z[iStep+1]-h_z[iStep]))+min_SOC*(capacity[iStep+1]-capacity[iStep]))

    @constraint(M, en_bal[iStep=1:NSteps+1], min_SOC + Beta*(x[iStep]+2*y[iStep]+4*z[iStep]) == soc[iStep]) 
    @constraint(M, en_square[iStep=1:NSteps+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*Beta*(x[iStep]+2*y[iStep]+4*z[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep])*Beta^2)
    
    # AUXILIARY FOR FORTET EXACT INEQUALITIES
    @constraint(M, xx_1[iStep=1:NSteps+1], w_xx[iStep] <= x[iStep])
    @constraint(M, xx_2[iStep=1:NSteps+1], w_xx[iStep] >= 2*x[iStep]-1)

    @constraint(M, xy_1[iStep=1:NSteps+1], w_xy[iStep] <= x[iStep])
    @constraint(M, xy_2[iStep=1:NSteps+1], w_xy[iStep] <= y[iStep])
    @constraint(M, xy_3[iStep=1:NSteps+1], w_xy[iStep] >= x[iStep]+y[iStep]-1)

    @constraint(M, xz_1[iStep=1:NSteps+1], w_xz[iStep] <= x[iStep])
    @constraint(M, xz_2[iStep=1:NSteps+1], w_xz[iStep] <= z[iStep])
    @constraint(M, xz_3[iStep=1:NSteps+1], w_xz[iStep] >= x[iStep]+z[iStep]-1)

    @constraint(M, yy_1[iStep=1:NSteps+1], w_yy[iStep] <= y[iStep])
    @constraint(M, yy_2[iStep=1:NSteps+1], w_yy[iStep] >= 2*y[iStep]-1)

    @constraint(M, zz_1[iStep=1:NSteps+1], w_zz[iStep] <= z[iStep])
    @constraint(M, zz_2[iStep=1:NSteps+1], w_zz[iStep] >= 2*z[iStep]-1)

    @constraint(M, zy_1[iStep=1:NSteps+1], w_zy[iStep] <= z[iStep])
    @constraint(M, zy_2[iStep=1:NSteps+1], w_zy[iStep] <= y[iStep])
    @constraint(M, zy_3[iStep=1:NSteps+1], w_zy[iStep] >= z[iStep]+y[iStep]-1)
    
    # AUXILIARY CONSTRAINTS FOR ENERGY BALANCE CONSTRAINTS
    @constraint(M, h_x_1[iStep=1:NSteps+1], h_x[iStep]>= min_SOH/min_SOH*x[iStep])
    @constraint(M, h_x_2[iStep=1:NSteps+1], h_x[iStep]<= max_SOH/min_SOH*x[iStep])
    @constraint(M, h_x_3[iStep=1:NSteps+1], h_x[iStep]>= capacity[iStep]-max_SOH/min_SOH*(1-x[iStep]))
    @constraint(M, h_x_4[iStep=1:NSteps+1], h_x[iStep]<= capacity[iStep]-min_SOH/min_SOH*(1-x[iStep]))

    @constraint(M, h_y_1[iStep=1:NSteps+1], h_y[iStep]>= min_SOH/min_SOH*y[iStep])
    @constraint(M, h_y_2[iStep=1:NSteps+1], h_y[iStep]<= max_SOH/min_SOH*y[iStep])
    @constraint(M, h_y_3[iStep=1:NSteps+1], h_y[iStep]>= capacity[iStep]-max_SOH/min_SOH*(1-y[iStep]))
    @constraint(M, h_y_4[iStep=1:NSteps+1], h_y[iStep]<= capacity[iStep]-min_SOH/min_SOH*(1-y[iStep]))

    @constraint(M, h_z_1[iStep=1:NSteps+1], h_z[iStep]>= min_SOH/min_SOH*z[iStep])
    @constraint(M, h_z_2[iStep=1:NSteps+1], h_z[iStep]<= max_SOH/min_SOH*z[iStep])
    @constraint(M, h_z_3[iStep=1:NSteps+1], h_z[iStep]>= capacity[iStep]-max_SOH/min_SOH*(1-z[iStep]))
    @constraint(M, h_z_4[iStep=1:NSteps+1], h_z[iStep]<= capacity[iStep]-min_SOH/min_SOH*(1-z[iStep]))

    #ADDITIONAL CONSTRAINTS FOR DEGRADATION (BILINEAR TERMS)
    @constraint(M, h_xx_1[iStep=1:NSteps+1], h_xx[iStep]>= min_SOH/min_SOH*w_xx[iStep])
    @constraint(M, h_xx_2[iStep=1:NSteps+1], h_xx[iStep]<= max_SOH/min_SOH*w_xx[iStep])
    @constraint(M, h_xx_3[iStep=1:NSteps+1], h_xx[iStep]>= capacity[iStep]-max_SOH/min_SOH*(1-w_xx[iStep]))
    @constraint(M, h_xx_4[iStep=1:NSteps+1], h_xx[iStep]<= capacity[iStep]-min_SOH/min_SOH*(1-w_xx[iStep]))

    @constraint(M, h_xy_1[iStep=1:NSteps+1], h_xy[iStep]>= min_SOH/min_SOH*w_xy[iStep])
    @constraint(M, h_xy_2[iStep=1:NSteps+1], h_xy[iStep]<= max_SOH/min_SOH*w_xy[iStep])
    @constraint(M, h_xy_3[iStep=1:NSteps+1], h_xy[iStep]>= capacity[iStep]-max_SOH/min_SOH*(1-w_xy[iStep]))
    @constraint(M, h_xy_4[iStep=1:NSteps+1], h_xy[iStep]<= capacity[iStep]-min_SOH/min_SOH*(1-w_xy[iStep]))

    @constraint(M, h_xz_1[iStep=1:NSteps+1], h_xz[iStep]>= min_SOH/min_SOH*w_xz[iStep])
    @constraint(M, h_xz_2[iStep=1:NSteps+1], h_xz[iStep]<= max_SOH/min_SOH*w_xz[iStep])
    @constraint(M, h_xz_3[iStep=1:NSteps+1], h_xz[iStep]>= capacity[iStep]-max_SOH/min_SOH*(1-w_xz[iStep]))
    @constraint(M, h_xz_4[iStep=1:NSteps+1], h_xz[iStep]<= capacity[iStep]-min_SOH/min_SOH*(1-w_xz[iStep]))

    @constraint(M, h_yy_1[iStep=1:NSteps+1], h_yy[iStep]>= min_SOH/min_SOH*w_yy[iStep])
    @constraint(M, h_yy_2[iStep=1:NSteps+1], h_yy[iStep]<= max_SOH/min_SOH*w_yy[iStep])
    @constraint(M, h_yy_3[iStep=1:NSteps+1], h_yy[iStep]>= capacity[iStep]-max_SOH/min_SOH*(1-w_yy[iStep]))
    @constraint(M, h_yy_4[iStep=1:NSteps+1], h_yy[iStep]<= capacity[iStep]-min_SOH/min_SOH*(1-w_yy[iStep]))

    @constraint(M, h_zz_1[iStep=1:NSteps+1], h_zz[iStep]>= min_SOH/min_SOH*w_zz[iStep])
    @constraint(M, h_zz_2[iStep=1:NSteps+1], h_zz[iStep]<= max_SOH/min_SOH*w_zz[iStep])
    @constraint(M, h_zz_3[iStep=1:NSteps+1], h_zz[iStep]>= capacity[iStep]-max_SOH/min_SOH*(1-w_zz[iStep]))
    @constraint(M, h_zz_4[iStep=1:NSteps+1], h_zz[iStep]<= capacity[iStep]-min_SOH/min_SOH*(1-w_zz[iStep]))

    @constraint(M, h_yz_1[iStep=1:NSteps+1], h_yz[iStep]>= min_SOH/min_SOH*w_zy[iStep])
    @constraint(M, h_yz_2[iStep=1:NSteps+1], h_yz[iStep]<= max_SOH/min_SOH*w_zy[iStep])
    @constraint(M, h_yz_3[iStep=1:NSteps+1], h_yz[iStep]>= capacity[iStep]-max_SOH/min_SOH*(1-w_zy[iStep]))
    @constraint(M, h_yz_4[iStep=1:NSteps+1], h_yz[iStep]<= capacity[iStep]-min_SOH/min_SOH*(1-w_zy[iStep]))

    #binary variable for operation
    #@constraint(M, charging[iStep=1:NSteps], e_charge[iStep] <= max_P*NHoursStep*(1-bin_op[iStep]))
    #@constraint(M, discharging[iStep=1:NSteps], e_discharge[iStep] <= max_P*NHoursStep*bin_op[iStep])

    # CONSTRAINTS ON DEGRADATION
    @constraint(M, deg_1[iStep=1:NSteps], deg[iStep] >= (2*min_SOC*Beta-2*Beta)*(h_x[iStep+1]-h_x[iStep]+2*(h_y[iStep+1]-h_y[iStep])+4*(h_z[iStep+1]-h_z[iStep]))+Beta^2*(h_xx[iStep+1]-h_xx[iStep]+4*(h_xy[iStep+1]-h_xy[iStep])+8*(h_xz[iStep+1]-h_xz[iStep])+4*(h_yy[iStep+1]-h_yy[iStep])+16*(h_zz[iStep+1]-h_zz[iStep])+16*(h_yz[iStep+1]-h_yz[iStep])))
    @constraint(M, deg_2[iStep=1:NSteps], deg[iStep] >= (2*min_SOC*Beta-2*Beta)*(h_x[iStep]-h_x[iStep+1]+2*(h_y[iStep]-h_y[iStep+1])+4*(h_z[iStep]-h_z[iStep+1]))+Beta^2*(h_xx[iStep]-h_xx[iStep+1]+4*(h_xy[iStep]-h_xy[iStep+1])+8*(h_xz[iStep]-h_xz[iStep+1])+4*(h_yy[iStep]-h_yy[iStep+1])+16*(h_zz[iStep]-h_zz[iStep+1])+16*(h_yz[iStep]-h_yz[iStep+1])))

    #CONSTRAINT ON REVAMPING
    @constraint(M, energy_capacity[iStage=1:NStages], capacity[Steps_stages[iStage]+2] == capacity[Steps_stages[iStage]+1]+revamping[iStage]-deg[Steps_stages[iStage]+1]*k_deg) #
   
    @constraint(M, initial_e[iStep=1], capacity[iStep] == min_SOH)

    @constraint(M,en_cap1[iStage in 1:NStages, iStep in ((Steps_stages[iStage]+2):Steps_stages[iStage+1])], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k_deg)

    @constraint(M, stop_charge[iStage in 2:NStages, iStep in (Steps_stages[iStage]:(Steps_stages[iStage]+Steps_stop[iStage-1]))], e_charge[iStep] <= (1-e[iStage])*max_P) 
    @constraint(M, stop_discharge[iStage in 2:NStages, iStep in (Steps_stages[iStage]:(Steps_stages[iStage]+Steps_stop[iStage-1]))], e_discharge[iStep] <= (1-e[iStage])*max_P) 

    @constraint(M, rev_3[iStage=1:NStages], capacity[Steps_stages[iStage]+2]>= capacity[Steps_stages[iStage]+1])
   #@constraint(M, rev[iStage=1:NStages], revamping[iStage] <= (max_SOH-min_SOH))
    @constraint(M, rev[iStage=1:NStages], revamping[iStage] <= (max_SOH-min_SOH)*e[iStage])
    @constraint(M, first_e[iStage=1], e[iStage] == 1)

    #CONSTRAINT SU VARIABILI AUSILIARIE PER ACQUISTO/VENDITA
    @constraint(M, vendita[iStage=1], rev_vendita[iStage] == 0)
    @constraint(M, vendita_1[iStage=2:NStages], rev_vendita[iStage] >= 0)
    @constraint(M, vendita_2[iStage=2:NStages], rev_vendita[iStage] >= capacity[Steps_stages[iStage]+1]- e[iStage]*Big)

    @constraint(M, acquisto_1[iStage=2:NStages], rev_acquisto[iStage] >= -(1-e[iStage])*Big)
    @constraint(M, acquisto_2[iStage=2:NStages], rev_acquisto[iStage] >= -capacity[Steps_stages[iStage]+2] )
    @constraint(M, acquisto_3[iStage=1], rev_acquisto[iStage] == 0)

    return BuildStageProblem_3(
        M,
        soc,
        soc_quad,
        e_charge,
        e_discharge,
        bin_op,
        deg,
        x,
        y,
        z,
        h_x,
        h_y,
        h_z,
        w_xx,
        w_yy,
        w_zz,
        w_xy,
        w_xz,
        w_zy,
        capacity,
        revamping,
        e,
        rev_vendita,
        rev_acquisto,
        h_xx,
        h_xy,
        h_xz,
        h_yy,
        h_zz,
        h_yz,
      )
end
