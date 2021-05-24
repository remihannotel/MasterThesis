using Ipopt 
using JuMP
using Interact 
using Plots 
using Blink 

# Parameters 
β = 0.99 # discount factor 
θ = 0.3 # share of capital in dirty sector 
γ = 0.4 # share of clean production 
α = 0.3 # share of capital in clean sector 
σ = 2  # coefficient of risk aversion  
L = 1  # labor force set to 1
A = 1  # total factor productivity 

######### BAU scenario ########

m = Model(with_optimizer(Ipopt.Optimizer))
@variable(m, 0.001 <= x <= 1)
@variable(m, 0.001 <= y <= 1)
@NLobjective(m, Max, (((1-γ*α-θ*(1-γ))*y^(γ*α)*x^(θ*(1-γ)) -x -y)^(1-σ) -1 )/(1-σ) + β*( (  (1+θ*(1-γ)*y^(γ*α)*x^(θ*(1-γ) -1))*x + (1 + γ*α*y^(γ*α -1)*x^(θ*(1-γ)))*y)^(1-σ) -1) /(1-σ))  
JuMP.optimize!(m)  
Kd = value(x)  #Kd = 0.097
Kc = value(y)  #Kc = 0.065 
@show termination_status(m)
has_values(m)

# Kd is dirty capital, Kc is clean capital 
# Now with the value of dirty capital, I can plug it into the energy sector function to get emissions, and then the change in temperature
# I create a bench of function for the environmental part.  

    # A function that outputs emissions related to dirty sector (expressed in 1000 GtC (Gigatone of Carbon) due to the length of time period)
function energy_emissions(K)
	return K^θ*L^(1-θ) #
end 
	
	# A function that outputs atmospheric concentration (in GtC)
function  atmos_concentration(E)  
	ϕ = 0.5 
	ξ = 0.3
	S1 = ϕ*ξ*1000*E + 438.175  
	S2 = (1-ϕ)*ξ*1000*E + 0.99*438.175
	St = S1 + S2 
	return St    
end 
	# A function that outputs change in temperature since pre-industrial period 
function varpreindustemps(S)
	λ = 3 
	Sbar = 581
	T = λ*log(S/Sbar)/log(2)
	return T 
end 
	
	# A function that outputs damages as a function of change in temperature relatively to pre-industrial temperature level (1850-1900)
	# κ1 and κ2 are parameters that are calibrated so that there is a tipping point around +6°C, from which damages increases dramatically.
function damage(T)   
	κ1 = 1/20.46
	κ2 = 1/6.081     
	damages = 1 - 1/(1+(κ1*T)^2 + (κ2*T)^6.754) 
	return damages 
end  

	# Results 
energy_bau = energy_emissions(Kd)  # ≈ 0.497
concentration_bau = atmos_concentration(energy_bau) # ≈ 1021.1
change_temp_preindus = varpreindustemps(concentration_bau) # ≈ 2.44
damage_BAU = damage(change_temp_preindus) # 1.61%
	# Plotting the damage function 
p1 = plot(0:0.1:7,x -> damage(x), label = "damage", title = "damage as a function of temperature increase", yaxis = "share of the production that is destroyed" , xaxis = "change in temperature from pre-industrial era")
	# Plotting the BAU value of damage in the function 
scatter!(p1,[2.44],[0.0161], label = "BAU damage when γ = 0.4")
scatter!(p1, [2.51], [0.0173], label = "BAU damage when γ = 0.2")
Ystar = A*Kc^(α*γ)*L^(1-γ*α-θ*(1-γ))*Kd^(θ*(1-γ)) #0.474
wagestar = (1-γ*α-θ*(1-γ))*Ystar/L #0.331  
rdstar = (θ*(1-γ))*Ystar/Kd 
rcstar = (α*γ)*Ystar/Kc 
Rdstar = 1 + rdstar #1.876 
Rcstar = 1 + rcstar #1.876 
c_youngstar = wagestar - Kd - Kc #0.169
c_oldstar = Rcstar*Kd + Rdstar*Kc #0.304 
Cstar = c_youngstar + c_oldstar #0.474 
utility_BAU = (c_youngstar^(1-σ)-1)/(1-σ) + β*(c_oldstar^(1-σ)-1)/(1-σ)   # ≈ - 7.172
damagestar = damage_BAU # 0.0161 

# After damages equilibrium : 
Y_after_damage = (1- damagestar)*Ystar #0.466
w_after_damage = (1-γ*α-θ*(1-γ))*Y_after_damage/L #0.326  
rd_after_damage = (θ*(1-γ))*Y_after_damage/Kd 
rc_after_damage = (α*γ)*Y_after_damage/Kc 
Rd_after_damage = 1 + rd_after_damage #1.862 
Rc_after_damage = 1 + rc_after_damage #1.862
c_y_after_damage = w_after_damage - Kd - Kc #0.164
c_o_after_damage = Rc_after_damage*Kd + Rd_after_damage*Kc #0.302  
Cstar_after_damage = c_y_after_damage + c_o_after_damage #0.466
utility_afterdamage = (c_y_after_damage^(1-σ)-1)/(1-σ) + β*(c_o_after_damage^(1-σ)-1)/(1-σ) # ≈ -7.389

#######  Taxation/Subsidy schemes   #######
####### The goal is to see the burden of climate change between firms and households via taxation and subsidy schemes.

## Burden on firms, with an increase of the price of dirty capital


function policyfirm(increase_p_input = 0)
	m1 = Model(with_optimizer(Ipopt.Optimizer))
	@variable(m1, 0.001 <= x <= 1)
	@variable(m1, 0.001 <= y <= 1)    
	@NLobjective(m1, Max, (((1-γ*α-θ*(1-γ))*y^(γ*α)*x^(θ*(1-γ)) -x -y)^(1-σ) -1 )/(1-σ) + β*(  ((1+ θ*(1-γ)*y^(γ*α)*x^(θ*(1-γ)-1)*(1 - increase_p_input))*x + (1 + γ*α*y^(γ*α -1)*x^(θ*(1-γ)))*y)^(1-σ) -1) /(1-σ) )
	@NLconstraint(m1,  1 + θ*(1-γ)*y^(γ*α)*x^(θ*(1-γ)-1)*(1- increase_p_input) >= 1)
	JuMP.optimize!(m1)
	Kd1 = value(x)
	Kc1 = value(y) 
	return Kd1
end

function policyfirm2(increase_p_input = 0)
		m1 = Model(with_optimizer(Ipopt.Optimizer))
		@variable(m1, 0.001 <= x <= 1)  
		@variable(m1, 0.001 <= y <= 1)  
		@NLobjective(m1, Max, (((1-γ*α-θ*(1-γ))*y^(γ*α)*x^(θ*(1-γ)) -x -y)^(1-σ) -1 )/(1-σ) + β*( (  (1+θ*(1-γ)*y^(γ*α)*x^(θ*(1-γ)-1)*(1- increase_p_input))*x + (1 + γ*α*y^(γ*α -1)*x^(θ*(1-γ)))*y)^(1-σ) -1) /(1-σ) )
		@NLconstraint(m1,  1+ θ*(1-γ)*y^(γ*α)*x^(θ*(1-γ)-1)*(1- increase_p_input) >= 1)
		JuMP.optimize!(m1)
		Kd1 = value(x)
		Kc1 = value(y) 
		return Kc1 
end 

p = plot(0:0.1:2,x -> policyfirm(x), label = "dirty capital", title = "capital when fiscal policy on firms", yaxis = "capital", xaxis = "price increase in dirty capital for firms")
plot!(p, 0:0.1:2, x -> policyfirm2(x), label = "clean capital")
plot!(p,0:0.1:2, x -> 0.08, seriestype = "hline", label = "threshold +2.4°C")

# find price increase for which Kd = 0.08
upper = policyfirm(0.7) # Kd = 0.081
lower = policyfirm(0.726) # Kd = 0.08 
price_increase = 0.726
Kd1 = lower # 0.080
Kc1 = policyfirm2(price_increase) # 0.084

# Results 
energy1 = energy_emissions(Kd1)  # ≈ 0.470
concentration1 = atmos_concentration(energy1) # ≈ 1012.85
change_temp_preindus1 = varpreindustemps(concentration1) # 2.4
damage1 = damage(change_temp_preindus1) # 1.55 % 

Ystar1 = A*Kc1^(α*γ)*L^(1-γ*α-θ*(1-γ))*Kd1^(θ*(1-γ)) #0.472
wagestar1 = (1-γ*α-θ*(1-γ))*Ystar1/L #0.331 
rdstar1 = (θ*(1-γ))*Ystar1/Kd1 
rcstar1 = (α*γ)*Ystar1/Kc1 
Rdstar1 = 1 + rdstar1 - price_increase #1.330
Rcstar1 = 1 + rcstar1 # 1.671
c_youngstar1 = wagestar1 - Kd1 - Kc1  # 0.166
c_oldstar1 = Rcstar1*Kd1 + Rdstar1*Kc1 # 0.247
Cstar1 = c_youngstar1 + c_oldstar1  #0.413
utility1 = (c_youngstar1^(1-σ)-1)/(1-σ) + β*(c_oldstar1^(1-σ)-1)/(1-σ)   # -8.056
government1 = Kd1*((rdstar1 + price_increase)/rdstar1 -1)*rdstar1 #0.0584 
damagestar1 = damage1 # 0.0155

Y_after_damage1 = (1- damagestar1)*Ystar1 #0.465
w_after_damage1 = (1-γ*α-θ*(1-γ))*Y_after_damage1/L # 0.325 
rd_after_damage1 = (θ*(1-γ))*Y_after_damage1/Kd1 
rc_after_damage1 = (α*γ)*Y_after_damage1/Kc1 
Rd_after_damage1 = 1 + rd_after_damage1 - price_increase # 1.314
Rc_after_damage1 = 1 + rc_after_damage1 # 1.661
c_y_after_damage1 = w_after_damage1 - Kd1 - Kc1 # 0.161
c_o_after_damage1 = Rc_after_damage1*Kd1 + Rd_after_damage1*Kc1 # 0.245 
Cstar_after_damage1 = c_y_after_damage1 + c_o_after_damage1 # 0.405
government1 = Kd1*((rd_after_damage1 + price_increase)/rd_after_damage1 -1)*rd_after_damage1 #0.0584 
utility_afterdamage1 = (c_y_after_damage1^(1-σ)-1)/(1-σ) + β*(c_o_after_damage1^(1-σ)-1)/(1-σ) # - 8.285

#damage curve with new point 
pdamage = plot(0:0.1:7,x -> damage(x), label = "damages", title = "damages as a function of temperature increase", yaxis = "share of the production that is destroyed", xaxis = "change in temperature from pre-industrial era")
scatter!(pdamage,[2.44],[0.0161], label = "BAU damage")
scatter!(pdamage, [2.4], [0.015], label = "damage when there is an increase in the price of the dirty input")

#Decrease of the damage but it is the same trend as in the 1st case. 
## Burden on households

function policy_households(capitalcost = 0)
	m1 = Model(with_optimizer(Ipopt.Optimizer))    
	@variable(m1, 0.001 <= x <= 1)                
	@variable(m1, 0.001 <= y <= 1)
	@NLobjective(m1, Max, (((1-γ*α-θ*(1-γ))*y^(γ*α)*x^(θ*(1-γ)) -x*(1+capitalcost) -y)^(1-σ) -1 )/(1-σ) + β*( ( (1+ θ*(1-γ)*y^(γ*α)*x^(θ*(1-γ)-1))*x + (1 + γ*α*y^(γ*α -1)*x^(θ*(1-γ)))*y   )^(1-σ) -1 ) /(1-σ) )
	JuMP.optimize!(m1)
	Kd1 = value(x)
	Kc1 = value(y) 
	return Kd1 
end 

function policy_households2(capitalcost = 0)
		m1 = Model(with_optimizer(Ipopt.Optimizer))
		@variable(m1, 0.001 <= x <= 1)              
		@variable(m1, 0.001 <= y <= 1)
		@NLobjective(m1, Max, (((1-γ*α-θ*(1-γ))*y^(γ*α)*x^(θ*(1-γ)) -x*(1+capitalcost) -y)^(1-σ) -1 )/(1-σ) + β*( ( (1+θ*(1-γ)*y^(γ*α)*x^(θ*(1-γ)-1))*x + (1 + γ*α*y^(γ*α -1)*x^(θ*(1-γ)))*y   )^(1-σ) -1 ) /(1-σ) )
		JuMP.optimize!(m1)
		Kd1 = value(x)
		Kc1 = value(y) 
		return Kc1 
end 

pl = plot(0:0.1:1, x -> policy_households(x), label = "dirty capital", title = "capital when fiscal policy on households", yaxis = "capital", xaxis = "increase in the price of dirty asset for households")
plot!(pl, 0:0.1:1, x -> policy_households2(x), label = "clean capital")
plot!(pl, 0:0.1:1, x -> 0.08, seriestype = "hline", label = "threshold + 2.4°C")

upper1 = policy_households(0.17)
lower1 = policy_households(0.15) #0.08
Kd2 = lower1 
Kc2 = policy_households2(0.15) # Kc = 0.065
price_increase2 = 0.15
# Results 
energy2 = energy_emissions(Kd2)  # ≈ 0.469
concentration2 = atmos_concentration(energy2) # ≈ 1012.6
change_temp_preindus2 = varpreindustemps(concentration2) # 2.4
damage2 = damage(change_temp_preindus2) # 1.5 % 

Ystar2 = A*Kc2^(α*γ)*L^(1-γ*α-θ*(1-γ))*Kd2^(θ*(1-γ)) #0.458
wagestar2 = (1-γ*α-θ*(1-γ))*Ystar2/L # 0.320
rdstar2 = (θ*(1-γ))*Ystar2/Kd2 
rcstar2 = (α*γ)*Ystar2/Kc2
Rdstar2 = 1 + rdstar2 #2.028
Rcstar2 = 1 + rcstar2 # 1.841
c_youngstar2 = wagestar2 - Kd2*(1+price_increase2) - Kc2 # 0.163 
c_oldstar2 = Rcstar2*Kd2 + Rdstar2*Kc2 # 0.280 
Cstar2 = c_youngstar2 + c_oldstar2 # 0.443
utility2 = (c_youngstar2^(1-σ)-1)/(1-σ) + β*(c_oldstar2^(1-σ)-1)/(1-σ)   # -7.687
government2 = Kd2*price_increase2 #0.012   
damagestar2 = damage2 # 1.5%

Y_after_damage2 = (1- damagestar2)*Ystar2 # 0.450
w_after_damage2 = (1-γ*α-θ*(1-γ))*(Y_after_damage2/L) # 0.315 
rd_after_damage2 = (θ*(1-γ))*Y_after_damage2/Kd2 # 1.013
rc_after_damage2 = (α*γ)*Y_after_damage2/Kc2 # 0.828
Rd_after_damage2 = 1 + rd_after_damage2 # 2.013
Rc_after_damage2 = 1 + rc_after_damage2 # 1.828 
c_y_after_damage2 = w_after_damage2 - Kd2*(1+price_increase2) - Kc2 # 0.158
c_o_after_damage2 = Rc_after_damage2*Kd2 + Rd_after_damage2*Kc2 # 0.278
Cstar_after_damage2 = c_y_after_damage2 + c_o_after_damage2 # 0.436
government_after_damage2 = Kd2*price_increase2 #0.012 
utility_afterdamage2 = (c_y_after_damage2^(1-σ)-1)/(1-σ) + β*(c_o_after_damage2^(1-σ)-1)/(1-σ) # -7.906 


## Price increase in dirty firm's input, with price increase for dirty asset for households 

function policy_firm_households(price_increase_firm = 0 ; capitalcost = 0)
	m1 = Model(with_optimizer(Ipopt.Optimizer))
	@variable(m1, 0.001 <= x <= 1)
	@variable(m1, 0.001 <= y <= 1)
	@NLobjective(m1, Max, (((1-γ*α-θ*(1-γ))*y^(γ*α)*x^(θ*(1-γ)) -x*(1+capitalcost) -y)^(1-σ) -1 )/(1-σ) + β*(  ((1+ θ*(1-γ)*y^(γ*α)*x^(θ*(1-γ)-1) - price_increase_firm)*x + (1 + γ*α*y^(γ*α -1)*x^(θ*(1-γ)))*y)^(1-σ) -1) /(1-σ) )
	@NLconstraint(m1,  1+ (θ*(1-γ)*y^(γ*α)*x^(θ*(1-γ)-1) - price_increase_firm) >= 1)
	JuMP.optimize!(m1)
	Kd1 = value(x)
	Kc1 = value(y) 
	return Kd1 
end 

function policy_firm_households2(price_increase_firm = 0 ; capitalcost = 0)
	m1 = Model(with_optimizer(Ipopt.Optimizer))
	@variable(m1, 0.001 <= x <= 1)
	@variable(m1, 0.001 <= y <= 1)
	@NLobjective(m1, Max, (((1-γ*α-θ*(1-γ))*y^(γ*α)*x^(θ*(1-γ)) -x*(1+capitalcost) -y )^(1-σ) -1 )/(1-σ) + β*(  ((1+ θ*(1-γ)*y^(γ*α)*x^(θ*(1-γ)-1) - price_increase_firm)*x + (1 + γ*α*y^(γ*α -1)*x^(θ*(1-γ)))*y)^(1-σ) -1) /(1-σ) )
	@NLconstraint(m1,  1+ (θ*(1-γ)*y^(γ*α)*x^(θ*(1-γ)-1) - price_increase_firm) >= 1)
	JuMP.optimize!(m1)
	Kd1 = value(x)
	Kc1 = value(y) 
	return Kc1
end 

ui_bis = @manipulate for cost in slider(0:0.01:1, value = 0, label = "increase in the price of dirty asset for households")
		p3 = plot(0:0.1:2,x -> policy_firm_households(x; capitalcost = cost), label = "dirty capital", title = "capital when policy parameters on households & firms", yaxis = "capital", xaxis = "price increase in dirty capital for firms")
		plot!(p3,0:0.1:2, x -> policy_firm_households2(x; capitalcost = cost), label = "clean capital")
		plot!(p3,0:0.1:2, x -> 0.08, seriestype = "hline", label = "threshold +2.4°C")
end 

w_both = Window() 
body!(w_both,ui_bis)


upper2 = policy_firm_households(0.6; capitalcost = 0.04)
lower2 = policy_firm_households(0.595; capitalcost = 0.04) #0.08
Kd3 = lower2 #0.08
Kc3 = policy_firm_households2(0.595; capitalcost = 0.04) # Kc = 0.08
price_increase3 = 0.595
# Results 
energy3 = energy_emissions(Kd3)  # ≈ 0.470
concentration3 = atmos_concentration(energy3) # ≈ 1012.7
change_temp_preindus3 = varpreindustemps(concentration3) # 2.4
damage3 = damage(change_temp_preindus3) # 1.5 % 

Ystar3 = A*Kc3^(α*γ)*L^(1-γ*α-θ*(1-γ))*Kd3^(θ*(1-γ)) # 0.469
wagestar3 = (1-γ*α-θ*(1-γ))*Ystar3/L # 0.328
rdstar3 = (θ*(1-γ))*Ystar3/Kd3
rcstar3 = (α*γ)*Ystar3/Kc3
Rdstar3 = 1 + rdstar3 - price_increase3 # 1.458
Rcstar3 = 1 + rcstar3 # 1.703
c_youngstar3 = wagestar3 - Kd3*(1+0.04) - Kc3  # 0.165
c_oldstar3 = Rcstar3*Kd3 + Rdstar3*Kc3 # 0.253
Cstar3 = c_youngstar3 + c_oldstar3 # 0.418
utility3 = (c_youngstar3^(1-σ)-1)/(1-σ) + β*(c_oldstar3^(1-σ)-1)/(1-σ)  # -7.986
government3 = Kd3*0.04 + Kd3*rdstar3*((rdstar3 + price_increase3)/rdstar3 -1) #0.051
damagestar3 = damage3 # 0.015

Y_after_damage3 = (1- damagestar3)*Ystar3 # 0.462
w_after_damage3 = (1-γ*α-θ*(1-γ))*(Y_after_damage3/L) # 0.323
rd_after_damage3 = (θ*(1-γ))*Y_after_damage3/Kd3 
rc_after_damage3 = (α*γ)*Y_after_damage3/Kc3
Rd_after_damage3 = 1 + rd_after_damage3 - price_increase3 # 1.442
Rc_after_damage3 = 1 + rc_after_damage3 # 1.691
c_y_after_damage3 = w_after_damage3 - Kd3*(1+0.04) - Kc3 # 0.160
c_o_after_damage3 = Rc_after_damage3*Kd3 + Rd_after_damage3*Kc3 # 0.251
Cstar_after_damage3 = c_y_after_damage3 + c_o_after_damage3 # 0.411
government3_after_damage3 = Kd3*0.04 + Kd3*rd_after_damage3*((rd_after_damage3 + price_increase3)/rd_after_damage3 -1) # 0.051
utility_afterdamage3 = (c_y_after_damage3^(1-σ)-1)/(1-σ) + β*(c_o_after_damage3^(1-σ)-1)/(1-σ) #-8.213


################################################################


# Welfare analysis 
# I compute the consumption-equivalent factor 
function welfare(u, u_bau)
    (u/u_bau)^(1/(1-σ)) - 1 
end 

c_e_factor_policyfirms = welfare(utility1, utility_BAU) # - 10.98%
c_e_factor_policyfirms_withdamage = welfare(utility_afterdamage1 ,utility_afterdamage) # ≈ -10.82% 
# increase by 10.82 % of the consumption under this policy to reach the same utility that arises in BAU scenario. 

c_e_factor_policyhouseholds = welfare(utility2, utility_BAU) # -6.70%
c_e_factor_policyhouseholds_withdamage = welfare(utility_afterdamage2,utility_afterdamage) # - 6.54%
# increase by 6.55 % of the consumption under this policy to reach the same utility that arises in BAU scenario. 

c_e_factor_policyboth = welfare(utility3,utility_BAU) #-10.20%
c_e_factor_policyboth_withdamage = welfare(utility_afterdamage3,utility_afterdamage) # - 10.04%
# increase by 9.96 % of the consumption under this policy to reach the same utility that arises in BAU scenario. 


