using Optim

"""
    Tf, T_m_w, T_i_w, N_ion, P_i, P_SPDC, F_0 = Memory_Match(n, L, effs_MM, F_targ, M, filt_params; single_click = true, r = 0, turn_off = true, P_i_input = nothing)

This function calculates the various figures-of-merit for the memory matching entanglement scheme.
The outputs are (time outputs are in seconds)

- Tf: the total average entanglement distribution time (T_ent in the main text)
- T_m_w: the average worst-case memory waiting time
- T_i_w: the average worst-case ion waiting time
- N_ion: number of times ion is driven while another is in stored state in same cavity (for SC only)
- P_i:  average number of ion photons emitted (assuming drive is not turned off after click)
- P_SPDC: BB SPDC photon probability (for a single mode)
- F_0: entanglement fidelity as P_SPDC goes to zero

# Arguments
- n: swap depth (n = 0 without central repeater; n = 1 with)
- L: total length in meters
- effs_MM = (η_m, η_d, η_i_io, η_FC, η_AFC, η_herald), where
- η_m: memory efficiency (storage + readout)
- η_d: detector efficiency (assumed same for all photon frequencies. Note this is easy to change in the code)
- η_i_io: ion efficiency (from input to output single photon efficiency)
- η_FC: frequency conversion efficiency (assumed equal for all frequencies. Note this is easy to change in the code.)
- η_AFC: additional efficiency of the AFC readout-matching sequence (bare readout loss is in η_m)
- η_herald: SPDC heralding efficiency
- F_targ: target fidelity
- M: backbone multiplexing number
- filt_params = (Λ_0, 𝜈, ɛ_m), where
- Λ_0: low-pass filter width
- 𝜈: frequency of oscillations in memory matching retrieval
- ɛ_m: parameter relating to cropped echo technique
- single_click: true for single click scheme (false for double click)
- r: dephasing parameter
- turn_off: turn off ion drives after clicks for ion node
- P_i_input: if given, use this ion photon probability directly; if `nothing` (default), P_i is chosen to minimize Tf.
"""
function Memory_Match(n, L, effs_MM, F_targ, M, filt_params; single_click = true, r = 0, turn_off = true, P_i_input = nothing)

    n in (0,1) || throw(ArgumentError("n must be 0 or 1")) #We only consider up to one repeater link in the middle of the chain
    
    L0 = L/2^n; #Distance in each fundamental BB link (L without central memory repeater [n=0], L/2 with central memory repeater [n=1])
    η_m, η_d, η_i_io, η_FC, η_AFC, η_herald = effs_MM #unpack efficiencies
    Λ_0, 𝜈, ɛ_m = filt_params; #unpack filter parameters
    k0= (1 + ɛ_m^2/2*Λ_0^4/(Λ_0^2+𝜈^2)^2)^(-1/2); #parameter associated with fast oscillation approximation in matching wavefunctions (see Sec. II of main text)
    η_filt = k0^(-2)/(1 + ɛ_m^2/2)
    F0_max = (1 + k0^4)/2
    F_targ <= F0_max || throw(ArgumentError("Target fidelity F_targ=$F_targ is unattainable: the maximum possible F₀ is $F0_max (reached as P_i → 0)."))
    L_att = 21.7*1e3 #attenuation length in meters
    v = 2*1e8 #speed of light in fiber
    η_i = η_i_io*η_d # efficiency of EN ion channel
    η_BB = η_d*exp(-L0/(2*L_att)); #efficiency of backbone generation
    η_mc = η_herald*η_m*η_d; #efficiency of memory-memory swaps
    η_mcs = η_AFC*η_filt*η_FC*η_mc #efficiency of EN memory channel
    
    x_par_pur = 2; #we assume the two links needed for purification in SC are generated sequentially.
    #x_par_pur = 3/2; #uncomment this if the two links needed for purification are generated in parallel.
    x_par_DC = 2; #we assume the two DC BB chains (upper and lower) are generated sequentially (note the SPDC Matching equations assume this implicitly)
    #x_par_DC = 3/2; #uncomment if the two DC BB chains (upper and lower) are generated in parallel
    γ_f = turn_off ? 1/2 : 1; #implementing the fast oscillation approximation
    γ_f_DC = turn_off ? 2/3 : 1; #implementing the fast oscillation approximation

    #Some helper functions here, see below these for main optimizer code block
    function Ω_SC(P_i) #Generates Ω for SC
        P_g = 1-P_i;
        if turn_off == false
            Ω =  abs(r-2) < 1e-8 ? k0^2*P_g^3*(log(1-P_i))^2/P_i^2/((1-γ_f*P_i)) : k0^2/(P_i^2*(1-r/2)^2)*(P_g^((1+r)/2) - P_g^(3/2))^2/(1-γ_f*P_i)
        else
            Ω = 4*k0^2/(9*P_i^2)*(1-P_g^(3/2))^2/(1-P_i*γ_f)
        end
        return Ω
    end

    function Rate_SC(P_i) #Calculates the rate and other relevant quantities for SC for given P_i
        F_0 = (1 + Ω_SC(P_i)^2)/2
        if n ==0 
            P_BB = η_BB/(1-η_BB)*(F_0-F_targ)/((4*F_0-1)*(1 - η_mcs + η_mcs/η_i*(1-η_i)*γ_f)) #found from solving perturbative expression for F in main text
            p_swap = η_i*η_mcs*P_i*(1 + P_i*(η_i/η_mcs - η_i - η_i*γ_f))
            P_pur = 1/p_swap^2*η_i^2*η_mcs^2/2*P_i^2*(1-γ_f*P_i)^2
            T_SL = L0/((1-(1-P_BB)^M)*v*p_swap)
            T_m_w = L/v
            T_i_w = T_SL
            N_ion = 1/p_swap
        else
            P_BB = η_BB/2/(1-η_BB)*(F_0-F_targ)/((4*F_0-1)*(2 - η_mcs -η_mc + γ_f*η_mcs/η_i*(1-η_i)))  #found from solving perturbative expression for F in main text
            p_swap1 = 1/2*(η_mcs + η_i*P_i*(2 - η_mcs*(1+γ_f)))
            p_swap2 = η_i*η_mc*P_i/2/p_swap1^2*(η_mcs + η_i*P_i*(2 - η_mc - η_mcs*(1+γ_f))) 
            P_pur = η_mcs^2*η_i^2*η_mc^2*P_i^2/8*(1-γ_f*P_i)^2/p_swap1^4/p_swap2^2
            T_SL = 3/2*L0/((1-(1-P_BB)^M)*v*p_swap1*p_swap2)
            T_m_w = L0/((1-(1-P_BB)^M)*v*p_swap1)
            T_i_w = T_SL + L0/((1-(1-P_BB)^M)*v*p_swap1)
            N_ion = 1/(p_swap1*p_swap2)
        end
        T_ent = x_par_pur*T_SL/P_pur
        P_SPDC = P_BB/(2*η_BB)
        return 1/T_ent, F_0, T_m_w, T_i_w, N_ion, P_SPDC
    end
            
    function Rate_DC(P_i) #Calculates the rate and other relevant quantities for DC for given P_i
        P_g = 1- P_i;
        Ω = abs(r-2) < 1e-8 ? k0^2/P_i^2*(P_g^2*log(P_g) + 1/2*(1-P_g^2)) : 2*k0^2/P_i^2/(1-r/2)*( 1/(1+r/2)*(1-P_g^(1+r/2)) - 1/2*(1-P_g^2))
        F_0 = (1 + Ω^2)/2
        if n ==0 
            p_swap = η_mcs^2*η_i^2/8*P_i^2
            P_BB = η_BB/(1-η_BB)*(F_0-F_targ)/((4*F_0-1)*(1 - η_mcs + η_mcs/η_i*(1-η_i)*γ_f_DC)); 
            T_ent = x_par_DC*L0/((1-(1-P_BB)^M)*v*p_swap);
            T_i_w = L/v
            T_m_w = L0/((1-(1-P_BB)^M)*v)
        else
            temp= 2- η_mcs - η_mc + η_mcs/η_i*(1-η_i)*γ_f_DC
            P_BB = η_BB/2/(1-η_BB)*(F_0-F_targ)/((4*F_0-1)*temp); #found from solving perturbative expression for F in main text
            p_swap1 = η_i*η_mcs*P_i/2*(1 - η_mcs/2 + η_mcs/(2*η_i)*γ_f_DC*(1 - η_i))
            p_swap2 = η_mcs^2*η_mc^2*η_i^2*P_i^2/(32*p_swap1^2)
            T_ent = x_par_DC*L0*3/2/((1-(1-P_BB)^M)*v*p_swap1*p_swap2)
            T_i_w = x_par_DC*L0/((1-(1-P_BB)^M)*v*p_swap1)
            T_m_w = x_par_DC*L0/((1-(1-P_BB)^M)*v*p_swap1) + L0/((1-(1-P_BB)^M)*v)
        end
        P_SPDC = P_BB/(2*η_BB)
        N_ion = nothing #not defined for DC
        return 1/T_ent, F_0, T_m_w, T_i_w, N_ion, P_SPDC
    end

    #Main optimization block
    rate_fun = single_click ? Rate_SC : Rate_DC
    optimize_P_i = P_i_input === nothing
    if !optimize_P_i
        0 < P_i_input < 1 || throw(ArgumentError("P_i_input = $P_i_input must lie in (0, 1)."))
    end

    P_i = if optimize_P_i
        res = optimize(x -> -rate_fun(1 / (1 + exp(x)))[1], -15, 15) #optimize over #x, where P_i = 1/(1+e^x)
        1 / (1 + exp(Optim.minimizer(res)))
    else
        P_i_input
    end

    Rate, F_0, T_m_w, T_i_w, N_ion, P_SPDC = rate_fun(P_i)
    Tf = 1 / Rate

    #Some error messages if something goes wrong
    if P_SPDC <= 0   # F_0 ≤ F_targ at this P_i
        msg = optimize_P_i ?
        "No P_i in (0, 1) reaches F_targ = $F_targ (F_targ is too close to F0_max = $F0_max)." :
        "P_i = $P_i gives F_0 = $F_0 ≤ F_targ = $F_targ; choose a smaller P_i."
        throw(ArgumentError(msg))
    end
    P_SPDC < 1/(2*η_BB) || throw(DomainError(P_SPDC,
    "Backbone success probability P_BB = 2η_BB·P_SPDC ≥ 1; the model is invalid here."))
    
    return Tf, T_m_w, T_i_w, N_ion, P_i, P_SPDC, F_0
end

"""
    Tf, T_m_w, T_i_w, N_ion, P_i, P_SPDC, P_b = SPDC_Match(n, L, effs_SPDC, F_targ, M, N_bin, T_EN; single_click = true, r = 0, N_cap = Inf, turn_off = true, turn_off_mem= true, P_inputs = (nothing, nothing), optimize_T = true)

This function calculates the various figures-of-merit for the SPDC matching entanglement scheme. The outputs are (time outputs are in seconds)

- Tf: the total average entanglement time (T_ent in the paper)
- T_m_w: the average worst-case memory waiting time
- T_i_w: the average worst-case ion waiting time
- N_ion: number of times ion is driven while another is in stored state in same cavity
- P_i: average number of ion photons emitted (assuming drive is not turned off after click)
- P_SPDC: probability of getting single photon out of SPDC per attempt per mode
- P_b: probability of SPDC photon in edge node per bin

# Arguments
- n: swap depth (n = 0 without central repeater; n = 1 with)
- L: total length in meters
- effs_SPDC = (η_m, η_d, η_i_io, η_FC, η_herald), where
- η_m: memory efficiency (storage + readout)
- η_d: detector efficiency
- η_i_io: ion efficiency (from input to output single photon efficiency)
- η_FC: frequency conversion efficiency (assumed equal for all)
- η_herald: SPDC heralding efficiency
- F_targ: target fidelity
- M: backbone multiplexing number
- N_bin: number of bins for SPDC in EN
- T_EN: effective duration between successive ion pulse attempts in EN (in seconds)
- single_click: true for single click scheme (false for double click)
- r: dephasing parameter
- N_cap: maximum N_ion is allowed to be, set to Inf by default
- turn_off: turn off ion drives after click(s) for EN
- turn_off_mem: turn off SPDC drive after click(s) for EN
- P_inputs: optional input for (P_i, P_b). These should be between 0 < P_b < 0.5, and 0 < P_i < 0.5/1 for SC/DC
- optimize_T: if false and P_inputs is given, no optimization is employed and P_inputs is used directly. If true, we optimize by minimizing T_ent as a function of P_i, P_b. If P_inputs is given, we use those as initial guesses. For simplicity we have not allowed for partial optimization (i.e., fixing one of the photon inputs and optimizing over the others), but the code can be straightforwardly modified to do this.
"""
function SPDC_Match(n, L, effs_SPDC, F_targ, M, N_bin,T_EN; single_click = true, r = 0, N_cap = Inf, turn_off = true, turn_off_mem= true, P_inputs = (nothing, nothing), optimize_T = true)
    #Note: see end of this function for the main routine (with or without optimization)

    #some errors for bad inputs
    n in (0,1) || throw(ArgumentError("n must be 0 or 1")) #We only consider up to one repeater link in the middle of the chain
    if all(x -> !isnothing(x), P_inputs)
        (0 < P_inputs[1] < (single_click ? 0.5 : 1.0) && 0 < P_inputs[2] < 0.5) || throw(ArgumentError("P_inputs outside allowed bounds"))
    end
    if any(isnothing, P_inputs) && optimize_T == false
        @warn("Need full input for photon probabilities if not optimizing. Reverting to optimization.",maxlog=1)
        optimize_T = true
    end
    if any(isnothing, P_inputs) && !all(isnothing, P_inputs) && optimize_T == true 
        @warn("Ignoring incomplete photon probability inputs.",maxlog=1)
    end

    η_m, η_d, η_i_io, η_FC, η_herald  = effs_SPDC #unpack efficiencies
    L_att = 21.7*1e3 #attenuation length in meters
    v = 2*1e8 #speed of light in fiber
    η_i = η_i_io*η_d # efficiency at cavity freq detection of ion photon
    η_b = η_FC*η_d #efficiency at cavity freq
    L0 = L/2^n; #Distance in each fundamental BB link (L without central memory repeater [n=0], L/2 with central memory repeater [n=1])
    η_BB = η_d*exp(-L0/(2*L_att)); #efficiency of backbone generation
    η_mc = η_herald*η_m*η_d; #efficiency of memory swaps
    η_mc < 1 || throw(ArgumentError("η_herald*η_m*η_d must be < 1 (consider introducing slight inefficiencies)"))
    α_fm = turn_off_mem ? 1/2 : 1
    γ_f = turn_off ? 1/2 : 1
    γ_f_DC = turn_off ? 2/3 : 1
    κ_0 = η_b*(N_bin-1)*α_fm
    bfac = n == 0 ? 1.0 : 0.5 #these following terms end up being useful later
    Dfac = 1 - η_b + η_b/η_i*(1-η_i)*N_bin*γ_f
    DfacDC = 1 - η_b + η_b/η_i*(1-η_i)*N_bin*γ_f_DC
    Aden = (n == 0 ? 12 : 36)*(1-η_mc)
    x_par_pur = 2; #we assume the two links needed for purification in SC are generated sequentially.
    #x_par_pur = 3/2; #uncomment this if the two links needed for purification are generated in parallel.
    #Note: this protocol assumes that the upper and lower BB links in DC are generated in series. 

    
    #Some helper functions here, see below for main optimization block
    
    #This function calculates the various rates and figures of merit for the SC protocols:
    #Outputs are: Rate, P_SPDC, T_i_w, T_m_w, N_ion
    #Inputs are: P_i, and P_b. Given in terms of the transformation P(u) = P_max/(1+e^{-u})
    function Rate_sc(u_in)
        P_i = 0.5/(1+exp(-u_in[1])); P_b = 0.5/(1+exp(-u_in[2])); β_1 = sqrt(P_b);
        #We find c2 by using the analytical perturbative expression for the average target fidelity 
        P_SPDC = ( (1 - F_targ - r*P_i*(γ_f - 1/2))/Aden - bfac*P_b*Dfac) / (1-η_BB) #covers n = 0 and n =1 cases
        c2 = (1-η_BB)*P_SPDC; P_BB = 2*P_SPDC*η_BB;
        #The rest of this function is just using the various expressions from the paper
        R_BB =  (1-(1-P_BB)^M)*v/L0; 
        p_en = η_i*P_i +  η_b*β_1^2*N_bin;
        R_EN = p_en/T_EN
        B1, B2, B3, B4, _ = gen_Bs(β_1,η_b,η_mc);
        Trχ_1 = η_mc/4*(1 + c2*(3-4*η_mc))
        Trχ_2 = 1 - η_mc/2 + c2*η_mc*(η_mc - 3/2)
        G0, G1, G2 = gen_Gs(c2,η_mc)
        C1 = B1*G1 + B2*G0; C2 = B1*G2 + B2*G1; C3 = B3*G1 + B4*G0; C4 = B3*G2 + B4*G1;
        H1 = n == 0 ? B1 : C1; H2 = n == 0 ? B2 : C2; H3 = n == 0 ? B3 : C3; H4 = n == 0 ? B4 : C4; #as described in the text, with a repeater we take B -> C
        L00 = η_b^2/4*N_bin^2*β_1^4*(1 - γ_f*P_i)^2*(C4*H3 + C3*H4)
        L11 = η_i^2/4*P_i^2*(C1*H2+C2*H1) + η_b*η_i*(1-η_i)*N_bin*β_1^2*γ_f*P_i^2*(C1*H4 + C2*H3)/2
        L01 = η_b*η_i/4*N_bin*β_1^2*P_i*(1- γ_f*P_i)*(H1*C4+H2*C3)+η_b^2/4*N_bin^2*β_1^4*(1-η_i)*γ_f*P_i*(C3*H4 + H3*C4)
        p_swap1 = (2*η_i*P_i*(B1*Trχ_1 + B2*Trχ_2)*(1-κ_0*β_1^2) + 2*η_b*β_1^2*N_bin*(1- η_i*γ_f*P_i - κ_0*β_1^2)*(B3*Trχ_1 + B4*Trχ_2))/p_en
        p_swap2 = (L00 + 2*L01 + L11)*2^(4+n)/(p_en^2*p_swap1^(n+1))
        P_pur = (L01^2 + L00*L11)*2^(9+2*n)/(p_en^2*p_swap1^(n+1)*p_swap2)^2
        if n == 0
            T_SL = 1/((p_swap1*p_swap2)*(2*R_EN+R_BB))*(1 + R_BB/R_EN*(1+2*p_swap1)/2 + 2*R_EN/(R_EN+R_BB)*(1 + R_EN/R_BB + p_swap1*R_BB/R_EN));
            T_i_w = T_SL + 1/(R_EN + R_BB) + 1/min(R_EN,R_BB)
            T_m_w = max(3/(2*R_EN), 1/(R_EN + R_BB) + 1/min(R_EN,R_BB))
        else
            T_SL = 3/((2*p_swap1*p_swap2)*(R_EN + R_BB))*(1 + R_BB/R_EN + R_EN/R_BB)
            T_i_w = T_SL + 1/R_BB + 1/(p_swap1*(R_EN + R_BB))*(1 + R_BB/R_EN + R_EN/R_BB)
            T_m_w = 1/min(R_EN,R_BB) + 1/(p_swap1*(R_EN + R_BB))*(1 + R_BB/R_EN + R_EN/R_BB)      
        end
        N_ion = (p_en*p_swap1*p_swap2)^(-1)
        T_ent = x_par_pur*T_SL/P_pur
        return 1/T_ent, P_SPDC, T_i_w, T_m_w, N_ion
    end

    #This function is useful for calculating P_SPDC in terms of F_targ
    function dc_params(u_in)
        P_i = 1/(1 + exp(-u_in[1])); P_b = 0.5/(1 + exp(-u_in[2])); P_g = 1 - P_i
        y_0 = η_b*(N_bin-1)*P_b/P_i
        if turn_off_mem #Calculate the I_0 parameter
            if abs(r-2) < 1e-8 #avoid removable singularity
                I0 = (1-2*y_0)/P_i^2*(1/2*(1 -P_g^2) + P_g^2*log(1-P_i))  + 2*y_0/P_i^2*(5/18*(1-P_g^3) + 1/3*P_g^3*log(1-P_i) + 1/2*(P_g^3 - P_g^2))
            elseif abs(r-4) < 1e-8 #limiting form
                I0 = (1-2*y_0)/P_i^2*(1 + 2*P_g^3 - 3*P_g^2)/3 + 2*y_0/P_i^2*(-1/4*(1-P_g^4) + 4/9*(1-P_g^3) + P_g^3/3*log(1-P_i))
            else
                I0 = (1-2*y_0)/P_i^2*(1 -2*P_g^(1+r/2) + P_g^2 - r/2*(1-P_g^2))/(1-r^2/4) + 2*y_0/P_i^2*( (1-P_g^(2+r/2))/((1-r/2)*(2+r/2)) - r/3*(1-P_g^3)/(1-r^2/4) - (P_g^(1+r/2) - P_g^3)/((1+r/2)*(2-r/2)));
            end
        else
            I0 =  abs(r-2) < 1e-8  ? (1-y_0*P_i)^2/P_i^2*(1/2*(1 -P_g^2) + P_g^2*log(1-P_i)) :  (1-y_0*P_i)^2/P_i^2*(1 - 2*P_g^(1+r/2) + P_g^2 - r/2*(1-P_g^2))/(1-r^2/4)
        end
        F_0 = (1 + I0^2)/2
        AdenDC = (n == 0 ? 4 : 12)*(1-η_mc)*(4*F_0-1)
        c2 = (F_0 - F_targ + 2*κ_0*P_b*(2*F_0-1))/AdenDC - bfac*P_b*DfacDC #covers n = 0 and n = 1 cases
        P_SPDC = c2/(1-η_BB)
        return P_i, P_b, P_SPDC, c2, I0
    end

    #This function calculates the various rates and figures of merit for the DC protocols:
    #Outputs are: Rate, P_SPDC, T_i_w, T_m_w
    #Inputs are: P_i, and P_b. Given in terms of the transformation P(u) = P_max/(1+e^{-u})
    function Rate_dc(u_in)
        P_i, P_b, P_SPDC, c2, I0 = dc_params(u_in)
        β_1 = sqrt(P_b);
        P_BB = 2*P_SPDC*η_BB;
        R_BB =  (1-(1-P_BB)^M)*v/L0; 
        p_en = η_i*η_b*N_bin*P_i*β_1^2;
        R_EN = p_en/T_EN;
        B1, B2, B3, B4, _ = gen_Bs(β_1,η_b,η_mc);
        Trχ_1 = η_mc/4*(1 + c2*(3-4*η_mc))
        Trχ_2 = 1 - η_mc/2 + c2*η_mc*(η_mc - 3/2)
        p_swap1 = 4*η_b*η_i*N_bin*β_1^2*P_i*(B1*Trχ_1 + B2*Trχ_2)*(B3*Trχ_1 + B4*Trχ_2)/p_en;
        p_swap11 = 1/p_en*η_i*η_b*N_bin*β_1^2*P_i*((B1+B3)*Trχ_1 + (B4+B2)*Trχ_2);
        G0, G1, G2 = gen_Gs(c2,η_mc)
        C1 = B1*G1 + B2*G0; C2 = B1*G2 + B2*G1; C3 = B3*G1 + B4*G0; C4 = B3*G2 + B4*G1; 
        H1 = n == 0 ? B1 : C1; H2 = n == 0 ? B2 : C2; H3 = n == 0 ? B3 : C3; H4 = n == 0 ? B4 : C4; 
        L00 = η_b^2/4*N_bin^2*β_1^4*(H3*C4 + H4*C3)*P_i^2*( η_i^2/16*(C1*H2 + H1*C2)*(1-κ_0*β_1^2)^4 + η_i*η_b/8*(1-η_i)*N_bin*β_1^2*γ_f_DC*(C1*H4 + C2*H3));
        L01 = η_b^2/4*N_bin^2*β_1^4*(H4*C1 + H3*C2)*P_i^2*( η_i^2/16*(C4*H1 + H2*C3)*(1-κ_0*β_1^2)^4 + η_i*η_b/8*(1-η_i)*N_bin*β_1^2*γ_f_DC*(C4*H3 + H4*C3));
        p_swap2 = 2^(9+2*n)*(L01 + L00)/(p_en^2*p_swap1^(n+1));
        if n == 0
            P_S1, P_S2, P_S3, P_P1, P_P2, P_P3 = gen_probs(R_BB,R_EN);
            T_P1 = 1/p_swap2*(1/R_EN + 1/p_swap1*(1/(2*R_EN) + 2/(2*R_EN + R_BB)));
            T_P2 = 1/p_swap2*(1/R_EN + 1/p_swap1*(p_swap11/(R_EN + R_BB) + 2/(2*R_EN + R_BB)));
            T_P3 = 1/p_swap2*(1/R_EN + 1/p_swap1*((1+p_swap11)/(R_EN + R_BB) + 1/(2*R_EN + R_BB)));
            p_s_swap1 = (C3*B4 + C4*B3 + C1*B2 + C2*B1 + 2*(C3*B2 + C4*B1));
            p_s_swap2 = 2^3/p_s_swap1*((C3*B4 + C4*B3)*(C1*B2+C2*B1) + (C3*B2 + C4*B1)*(C1*B4 + C2*B3))
            T_S1 = 1/p_s_swap2*(1/R_BB + 1/p_s_swap1*(1/R_BB + 1/(R_EN + R_BB) + 1/(2*R_EN + R_BB)));
            T_S2 = 1/p_s_swap2*(1/R_BB + 1/p_s_swap1*((1+p_swap11)/(R_EN + R_BB) + 1/(2*R_EN + R_BB)));
            T_S3 = 1/p_s_swap2*(1/R_BB + 1/p_s_swap1*(p_swap11/(R_EN+R_BB) + 2/(2*R_EN+R_BB)));
            T_ent = P_S1*T_S1 + P_S2*T_S2 + P_S3*T_S3 + P_P1*T_P1 + P_P2*T_P2 + P_P3*T_P3;
            T_i_w = 1/(R_EN + R_BB) + max(2/R_BB, 1/R_EN + 1/(R_BB + R_EN));
            T_temp1 =1/(R_EN + R_BB) + max(1/(R_EN +R_BB) + 1/R_EN, 2/R_BB);
            T_temp2 =1/(2*R_EN + R_BB) + max(1/(R_EN + R_BB) + 1/min(R_BB,R_EN), 3/(2*R_EN));
            T_m_w = max(T_temp1,T_temp2);  
        else
            T_ent = 3/2/(p_swap1*p_swap2)*(1/(R_EN + R_BB)*(1 + R_BB/R_EN + R_EN/R_BB) + p_swap11/R_BB);
            T_i_w = 2/R_BB + 2/3*T_ent*p_swap2;
            T_m_w = max(2/R_BB, 1/(R_BB + R_EN) + 1/R_EN) + 2/3*T_ent*p_swap2;
        end
        return 1/T_ent, P_SPDC, T_i_w, T_m_w
    end

    #Function used in the optimization for the SC protocols. Keeps the rate and photon numbers physical by imposing penalty if they go negative or above perturbative region.
    function objective_sc(u)
        P_i =  0.5/(1 + exp(-u[1])); P_b = 0.5/(1 + exp(-u[2]));
        P_SPDC = ( (1 - F_targ - r*P_i*(γ_f - 1/2))/Aden - bfac*P_b*Dfac) / (1-η_BB)
        if !(0 < P_SPDC < 0.5 && isfinite(P_SPDC)) #impose penalty for leaving perturb. region
            return 1e10
        end
        Rate, _, _, _, N_i = Rate_sc(u)
        (isfinite(Rate) && Rate > 0) || return 1e10
        return -Rate * (1 - 100*max(0, N_i/N_cap - 1)^2) #for physical rate, want to minimize -rate. This line also caps approximately N_ion to be less than N_cap, if desired.
    end
    
    #Function used in the optimization for the DC protocols. Keeps the rate and photon numbers physical by imposing penalty if they go negative or above perturbative region.
    function objective_dc(u)
        P_i, P_b, P_SPDC, _, _ = dc_params(u)
        if !(0 < P_SPDC < 0.5 && isfinite(P_SPDC)) #impose penalty for leaving perturb. region
            return 1e10
        end
        Rate = first(Rate_dc(u))
        return isfinite(Rate) && Rate > 0 ? -Rate : 1e10 #for physical rate, want to minimize -rate
    end
    
    #Main optimization/calculation logic here
    if optimize_T
        if single_click #SC branch
            if any(isnothing, P_inputs) #If we have no initial guesses for photon numbers
                P_i_guess   = 0.1 #These are essentially random guesses, but they seem to work. Can be tweaked.
                P_b_guess = 1e-3 #These are essentially random guesses, but they seem to work. Can be tweaked.
                initial_u = [log(P_i_guess/(0.5 - P_i_guess)), log(P_b_guess/(0.5-P_b_guess))]
            else #if photon number initial guesses are given in P_input
                initial_u = [log(P_inputs[1]/(0.5-P_inputs[1])), log(P_inputs[2]/(0.5-P_inputs[2]))]
            end
            # If the initial guess is infeasible, the penalty is flat there and the optimizer cannot move.
            # Step both P_i and P_b down (feasible region exists at small P_i, P_b) until we start inside it.
            for _ in 1:50
                objective_sc(initial_u) < 1e10 && break
                initial_u .-= 1.0
            end
            res_rough = optimize(objective_sc, initial_u, NelderMead()) #we first use a non-gradient descent optimization
            res_fine  = optimize(objective_sc, Optim.minimizer(res_rough), LBFGS()) #we then feed that output as an initial guess into a gradient descent optimizer to refine the optimization
            u_out = Optim.minimum(res_fine) <= Optim.minimum(res_rough) ? Optim.minimizer(res_fine) : Optim.minimizer(res_rough) #take the better optimization
            Rate, P_SPDC, T_i_w, T_m_w, N_out = Rate_sc(u_out)
            if N_out > 1.05*N_cap #The 5% cap is arbitrary; keeping it kind of big allows the optimizer to operate smoother
                @warn("No solution with N_ion less than 5% larger than N_cap found." *
                      " N_cap might be below the " *
                      "smallest N_ion reachable at this fidelity; maybe increase N_cap or lower F_targ.",
                    maxlog=1)
            end
            Tf = 1 / Rate; P_i = 0.5/(1 + exp(-u_out[1])); P_b = 0.5/(1 + exp(-u_out[2]))
        else #DC branch
            if any(isnothing, P_inputs) #If we have no initial guesses for photon numbers
                P_i_guess   = 0.5 #These are essentially random guesses, but they seem to work. Can be tweaked.
                P_b_guess = 1e-4 #These are essentially random guesses, but they seem to work. Can be tweaked.
                initial_u = [log(P_i_guess/(1-P_i_guess)), log(P_b_guess/(0.5-P_b_guess))]
            else
                initial_u = [log(P_inputs[1]/(1-P_inputs[1])), log(P_inputs[2]/(0.5-P_inputs[2]))]
            end
            # If the initial guess is infeasible, the penalty is flat there and the optimizer cannot move.
            # Step both P_i and P_b down (feasible region exists at small P_i, P_b) until we start inside it.
            for _ in 1:50
                objective_dc(initial_u) < 1e10 && break
                initial_u .-= 1.0
            end
            res_rough = optimize(objective_dc, initial_u, NelderMead()) #we first use a non-gradient descent optimization
            res_fine = optimize(objective_dc, Optim.minimizer(res_rough), LBFGS()) #we then feed that output as an initial guess into a gradient descent optimizer to refine the optimization
            u_out = Optim.minimum(res_fine) <= Optim.minimum(res_rough) ? Optim.minimizer(res_fine) : Optim.minimizer(res_rough) #take the better one
            P_i = 1/(1+exp(-u_out[1])); P_b = 0.5/(1+exp(-u_out[2])); 
            Rate, P_SPDC, T_i_w, T_m_w = Rate_dc(u_out); N_out = nothing;
            Tf = 1/Rate;
        end
    else #no optimization
        P_i, P_b = P_inputs
        if single_click
            u_in = [log(P_i/(0.5-P_i)), log(P_b/(0.5-P_b))]
            Rate, P_SPDC, T_i_w, T_m_w, N_out = Rate_sc(u_in)
        else
            u_in = [log(P_i/(1-P_i)), log(P_b/(0.5-P_b))]
            Rate, P_SPDC, T_i_w, T_m_w = Rate_dc(u_in); N_out = nothing;
        end
        Tf = 1/Rate
    end

    #Some error messages if something went wrong
    if P_SPDC <= 0
        msg = optimize_T ?
            "No (P_i, P_b) found that reaches F_targ = $F_targ." :
            "P_i = $P_i, P_b = $P_b give P_SPDC ≤ 0, so F_targ = $F_targ is not reached; try smaller values."
        throw(ArgumentError(msg))
    end
    P_SPDC < 1/(2*η_BB) || throw(DomainError(P_SPDC,
        "Backbone success probability P_BB = 2η_BB·P_SPDC ≥ 1; the model is invalid here."))
        (isfinite(Rate) && Rate > 0) || throw(DomainError(Rate,
        "Unphysical rate at P_i = $P_i, P_b = $P_b."))
    return Tf, T_m_w, T_i_w, N_out, P_i, P_SPDC, P_b
end


##Helper functions

function gen_Bs(β_1,η_b,η_mc)
    B1 = 1 - β_1^2 + (1- η_b)*(1- η_mc)*β_1^2;
    B2 = η_mc/2*(1 - η_b)*β_1^2;
    B3 = (1 - η_mc)*(1+ 2*β_1^2*(1-η_b)*(1-η_mc))
    B4 = η_mc/2*(1 + 4*β_1^2*(1-η_b)*(1-η_mc))
    Bc = 1 - 1/2*β_1^2  + 2*β_1^2*(1-η_b)*(1-η_mc);
    return B1, B2, B3, B4, Bc
end

function gen_Gs(c2,η_mc)
    G0 = (1-η_mc)*(1 - 3*η_mc*c2);
    G1 = η_mc/4*(1 + 3*c2*(1-2*η_mc))
    G2 = η_mc^2/4*c2
    return G0, G1, G2
end

function gen_probs(R_BB,R_EN)
    P_S1 = 2*R_EN^2/((2*R_EN + R_BB)*(R_EN + R_BB));
    P_S2 = 2*R_EN^2*R_BB/((2*R_EN + R_BB)*(R_EN + R_BB)^2);
    P_S3 = 2*R_EN^2*R_BB/((2*R_EN + R_BB)^2*(R_EN + R_BB));
    P_P1 = (R_BB/(2*R_EN + R_BB))^2;
    P_P2 = 2*R_EN*R_BB^2/((2*R_EN + R_BB)^2*(R_EN + R_BB));
    P_P3 = 2*R_EN*R_BB^2/((2*R_EN + R_BB)*(R_EN + R_BB)^2);
    return P_S1, P_S2, P_S3, P_P1, P_P2, P_P3
end
