
### Set parameters ### 
N_y = 1     # young population 
N_o = 1     # old population 
θ = 0.3     # share of capital in production       
β = 0.99     # discount factor
δ = 0     # depreciation of capital 
α = 0.3   #share of clean capital in capital production
γ = 0.2
τ_l = 0   # taxation on labor 
τ_k = 0  # taxation on capital 
L = 1       # labor normalized to 1 
t_y = 0     # no transfers for young people  
t_o = 0     # no transfers for old people 

f = K -> -K + (N_y*(β/(1+β))*(1-θ)L^(-θ)*(γ/(1-γ))^(α*γ))^(1/(1-θ))

# I define a lower and an upper bound for K 
# to find the root of the above equation

# First : building the grid 
Kmin = 0.0001 # lower bound for K
Kmax = 0.3 # upper bound for K
n = 200 # number of grid points  

Kgrid = Kmin:(Kmax-Kmin)/(n-1):Kmax # equispaced grid for K
f_K = zeros(length(Kgrid)) # output vector with 0s 

# Now for each element i in the grid filled with zero, I replace it by f(K). 
for (i, K) in enumerate(Kgrid) 
    f_K[i] = f(K)
end 

# We plot f(K) and an horizontal line for f(K) = 0, to see the intersection. 
using Plots
plot(Kgrid, f_K,label = "Capital")  # Plotting f(K)
plot!([0], seriestype="hline", label="f(K) = 0")  # Plotting f(K) = 0

# We see on the graph that f(K) = 0 when K ≈ 0.15 
plot!([0.15], seriestype="vline", label="f(0.15)")

# Solving the model 
using Roots

Kstar = find_zero(f, (Kmin, Kmax)) # = 0.1967605106610738

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
