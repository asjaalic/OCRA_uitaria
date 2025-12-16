# EXCEL SAVINGS
#using DataFrames
#using XLSX

function data_saving(InputParameters::InputParam,ResultsOpt::Results_3)

   @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;       #NSteps,NHoursStage  
   @unpack (e_charge, e_discharge, rev, cap, soc, soc_quad, deg, deg_stage, gain_stage, cost_rev, x, y, z, h_x, h_y, h_z, w_xx, w_yy, w_zz, w_xy, w_xz, w_zy, bin_op) = ResultsOpt;
   #@unpack (e_charge, e_discharge, rev, cap, soc, soc_quad, deg, aux_deg, deg_stage, gain_stage, cost_rev, e, rev_vendita, rev_acquisto, x, y, z, h_x, h_y, h_z, w_xx, w_yy, w_zz, w_xy, w_xz, w_zy) = ResultsOpt;
   @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull,fix,cost ) = Battery ; 

    hour=string(now())
    a=replace(hour,':'=> '-')

    nameF= "Opzione 2, min SoC $min_SOC, max SoH $max_SOH mid costs"
    nameFile="Summary " 

    folder = "$nameF"
    mkdir(folder)
    cd(folder)
    main=pwd()

    general = DataFrame()
    #general_1 = DataFrame()
    battery_costs= DataFrame()

    capacity=zeros(NStages)
    capacity_finale = zeros(NStages)

    for iStage=1:NStages
        capacity[iStage]=cap[(Steps_stages[iStage]+2)]
    end

    for iStage=2:NStages
        capacity_finale[iStage-1] = cap[Steps_stages[iStage]+1]
    end
        capacity_finale[end] = cap[end]
    
    general[!,"Stage"] = 1:1:NStages
    general[!,"Initial Capacity MWh"] = capacity[:]
    general[!,"Final Capacity MWh"] = capacity_finale[:]
    general[!,"Revamping MWh"] = rev[:]
    #general[!, "Binary revamp"] = e[:]
    general[!,"Degradation"] = deg_stage[:]
    #general[!,"Net_Revenues"] = revenues_per_stage[:]
    general[!,"Gain charge/discharge"] = gain_stage[:]
    general[!,"Cost revamping"] = cost_rev[:]
   # general[!, "Aux acquisto"] = rev_acquisto[:]
   # general[!, "Aux vendita"] = rev_vendita[:]

    battery_costs[!,"Purchase costs €/MWh"] = Battery_price_purchase[1:NStages+1]    #Battery_price_purchase
    battery_costs[!,"Sale costs €/MWh"] = Battery_price_sale[1:NStages+1] 

    XLSX.writetable("$nameFile.xlsx", overwrite=true,                                       #$nameFile
    results_stages = (collect(DataFrames.eachcol(general)),DataFrames.names(general)),
    #results_1 = (collect(DataFrames.eachcol(general_1)),DataFrames.names(general_1)),
    costs = (collect(DataFrames.eachcol(battery_costs)),DataFrames.names(battery_costs)),
    )

    for iStage=1:NStages
        steps = DataFrame()

        steps[!,"Step"] = (Steps_stages[iStage]+1):(Steps_stages[iStage+1])
        steps[!, "Energy_prices €/MWh"] = Power_prices[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Energy capacity MWh"] = cap[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "SOC MWh"] = soc[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Charge MW"] = e_charge[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Discharge MW"] = e_discharge[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "SOC_quad MWh"] = soc_quad[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Deg MWh"] = deg[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "Bin operation "] = bin_op[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]

        steps[!, "x"] = x[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "y"] = y[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "z"] = z[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "h_x"] = h_x[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "h_y"] = h_y[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]
        steps[!, "h_z"] = h_z[(Steps_stages[iStage]+1):(Steps_stages[iStage+1])]

        XLSX.writetable("$iStage stage $a.xlsx", overwrite=true,                                       #$nameFile
        results_steps = (collect(DataFrames.eachcol(steps)),DataFrames.names(steps)),
        )

    end

    cd(main)             # ritorno nella cartella di salvataggio dati


    return println("Saved data in xlsx")
end


