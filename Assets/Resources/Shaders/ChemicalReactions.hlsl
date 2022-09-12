void UpdateChemicalReactions(uint3 id)
{          
    // Read the fluid source velocity : we will need it for each advection 
#ifdef FLOW
    float2 v = src[xy(id)].v;
#endif

    // =============
    // 1) GRAY-SCOTT
    // =============

    float A1 = src[xy(id)].A1;
    float dA1_dt = 0;
    float B1 = src[xy(id)].B1;
    float dB1_dt = 0;
    
    // Chemical power released
    #ifdef HEAT    
        float chemicalPower = 0;
    #endif

    // Diffusion
    float A1_mx = src[mx(id)].A1;
    float A1_px = src[px(id)].A1;
    float A1_my = src[my(id)].A1;
    float A1_py = src[py(id)].A1;
    float lap_A1 = (A1_px + A1_mx + A1_py + A1_my - 4 * A1) / (dx * dx);                
    
    float B1_mx = src[mx(id)].B1;
    float B1_px = src[px(id)].B1;
    float B1_my = src[my(id)].B1;
    float B1_py = src[py(id)].B1;
    float lap_B1 = (B1_px + B1_mx + B1_py + B1_my - 4 * B1) / (dx * dx);  

    dA1_dt += A_diffusion * lap_A1;
    dB1_dt += B_diffusion * lap_B1;

    // Advection
#ifdef FLOW        
    float advA1 = 0;
    float advB1 = 0;
    advA1 = Advection(v + float2(VX0, VY0), A1, A1_mx, A1_px, A1_my, A1_py);
    advB1 = Advection(v + float2(VX0, VY0), B1, B1_mx, B1_px, B1_my, B1_py); 
    dA1_dt += - advA1;
    dB1_dt += - advB1;
#endif        

    // GS REACTIONS
    // Effective supply and removal
    float GS_F_eff = GS_F;
    float GS_R_eff = GS_R;
    // // Adjust GS supply and depletion to generate a phase diagram
    // GS_F_eff = 0.02 + 0.04 * (1.0 * id.x / RESX);
    // GS_R_eff = 0.05 + 0.02 * (1.0 * id.y / RESY);

    // Infeed for A1 and decay for B1
    float A_FEED = 1;
    dA1_dt += GS_F_eff * (A_FEED - A1);             
    dB1_dt += - (GS_F_eff + GS_R_eff) * B1;

    // GS autocatalytic reaction     
    float k = k0;
#ifdef HEAT 
    float T = src[xy(id)].T;       
    if(T>0) k = k0 * exp(- Ea / T);     // Arrhenius
#endif
    float rate_GS1 = k * A1 * B1 * B1;
    dA1_dt += - rate_GS1;
    dB1_dt += rate_GS1;
#ifdef HEAT        
    chemicalPower += rate_GS1 * DeltaH1;    // Enthalpy
#endif
        


    // ===============
    // 2) GRAY-SCOTT 2
    // ===============

#ifdef GS2        
    float A2 = src[xy(id)].A2;
    float B2 = src[xy(id)].B2;
    float dA2_dt = 0;
    float dB2_dt = 0;

    float A2_mx = src[mx(id)].A2;
    float A2_px = src[px(id)].A2;
    float A2_my = src[my(id)].A2;
    float A2_py = src[py(id)].A2;   
    float lap_A2 = (A2_px + A2_mx + A2_py + A2_my - 4 * A2) / (dx * dx);

    float B2_mx = src[mx(id)].B2;
    float B2_px = src[px(id)].B2;
    float B2_my = src[my(id)].B2;
    float B2_py = src[py(id)].B2;    
    float lap_B2 = (B2_px + B2_mx + B2_py + B2_my - 4 * B2) / (dx * dx);

    dA2_dt += A_diffusion * lap_A2;    
    dB2_dt += B_diffusion * lap_B2;

    #ifdef FLOW
        float advA2 = 0;
        float advB2 = 0;
        advA2 = Advection(v + float2(VX0, VY0), A2, A2_mx, A2_px, A2_my, A2_py);  
        advB2 = Advection(v + float2(VX0, VY0), B2, B2_mx, B2_px, B2_my, B2_py);
        dA2_dt += - advA2;
        dB2_dt += - advB2;
    #endif

    // Infeed for A2 and decay for B2
    dA2_dt += GS_F_eff * (A_FEED - A2);             
    dB2_dt += - (GS_F_eff + GS_R_eff) * B2;

    float rate_GS2 = k * A2 * B2 * B2;
    dA2_dt += - rate_GS2;
    dB2_dt += + rate_GS2;
    #ifdef HEAT        
        chemicalPower += rate_GS2 * DeltaH2;
    #endif

    // GS2 will not be affected by other reactions : we write it now
    dest[xy(id)].A2 = max(0, A2 + dt * dA2_dt);
    dest[xy(id)].B2 = max(0, B2 + dt * dB2_dt);

#endif // GS2

// Chemical power won't change anymore : write it
#ifdef HEAT
    dest[xy(id)].chem = chemicalPower;
#endif



    // ==================
    // 3) LEARNING SYSTEM
    // ==================

#ifdef LEARNING

    float Tox = src[xy(id)].Tox;
    float dTox_dt = 0;
    float Ant = src[xy(id)].Ant;
    float dAnt_dt = 0;
    float Sti = src[xy(id)].Sti;
    float dSti_dt = 0;
    float STM_S = src[xy(id)].STM_S;
    float dSTM_S_dt = 0;
    float LTM = src[xy(id)].LTM;
    float dLTM_dt = 0;

    // TRANSPORT : Diffusion
    float Sti_mx = src[mx(id)].Sti;
    float Sti_px = src[px(id)].Sti;
    float Sti_my = src[my(id)].Sti;
    float Sti_py = src[py(id)].Sti;    
    float lap_Sti = (Sti_px + Sti_mx + Sti_py + Sti_my - 4 * Sti) / (dx * dx);
    dSti_dt += Sti_diffusion * lap_Sti;

    float Ant_mx = src[mx(id)].Ant;
    float Ant_px = src[px(id)].Ant;
    float Ant_my = src[my(id)].Ant;
    float Ant_py = src[py(id)].Ant;    
    float lap_Ant = (Ant_px + Ant_mx + Ant_py + Ant_my - 4 * Ant) / (dx * dx);
    dAnt_dt += Ant_diffusion * lap_Ant;

    float Tox_mx = src[mx(id)].Tox;
    float Tox_px = src[px(id)].Tox;
    float Tox_my = src[my(id)].Tox;
    float Tox_py = src[py(id)].Tox;    
    float lap_Tox = (Tox_px + Tox_mx + Tox_py + Tox_my - 4 * Tox) / (dx * dx);
    dTox_dt += Tox_diffusion * lap_Tox;

    float STM_S_mx = src[mx(id)].STM_S;
    float STM_S_px = src[px(id)].STM_S;
    float STM_S_my = src[my(id)].STM_S;
    float STM_S_py = src[py(id)].STM_S;    
    float lap_STM_S = (STM_S_px + STM_S_mx + STM_S_py + STM_S_my - 4 * STM_S) / (dx * dx);
    dSTM_S_dt += STM_diffusion * lap_STM_S;

    float LTM_mx = src[mx(id)].LTM;
    float LTM_px = src[px(id)].LTM;
    float LTM_my = src[my(id)].LTM;
    float LTM_py = src[py(id)].LTM;    
    float lap_LTM = (LTM_px + LTM_mx + LTM_py + LTM_my - 4 * LTM) / (dx * dx);
    dLTM_dt += LTM_diffusion * lap_LTM;  

    // TRANSPORT : Advection
    #ifdef FLOW    
        float advSti = 0;    
        advSti = Advection(v, Sti, Sti_mx, Sti_px, Sti_my, Sti_py);
        dSti_dt += - advSti;

        float advTox = 0;
        advTox = Advection(v, Tox, Tox_mx, Tox_px, Tox_my, Tox_py);
        dTox_dt += - advTox;
        
        float advAnt = 0;
        advAnt = Advection(v, Ant, Ant_mx, Ant_px, Ant_my, Ant_py);
        dAnt_dt += - advAnt;

        float advSTM_S = 0;
        advSTM_S = Advection(v, STM_S, STM_S_mx, STM_S_px, STM_S_my, STM_S_py);
        dSTM_S_dt += - advSTM_S;

        float advLTM = 0;
        advLTM = Advection(v, LTM, LTM_mx, LTM_px, LTM_my, LTM_py);
        dLTM_dt += - advLTM;
    #endif 

    
    // DECAYS
    dSti_dt += - Sti_decay * Sti;      
    dTox_dt += - Tox_decay * Tox;
    dSTM_S_dt += - STM_decay * STM_S;
    dLTM_dt += - LTM_decay * LTM;
    dAnt_dt += - Ant_decay * Ant;
    

    // CHEMICAL REACTIONS
    // ------------------
      

    // *** Toxin reacts with B1
    // ========================

    float k_ToxB_eff = k_TB;

    dB1_dt -= k_ToxB_eff * Tox * B1;
    dTox_dt -= k_ToxB_eff * Tox * B1; 


    // *** Antidote production
    // =======================

    if(memory_mechanism == 1)   // Direct network
    {
        dAnt_dt += kD * Tox * B1;
    }
    
    if(memory_mechanism == 2)   // Preemptive network
    {
        dAnt_dt += kP * Sti * B1;
    }

    if(memory_mechanism == 3)   // Associative
    {
        // STM production catalyzed by stimulus in the presence of B
        dSTM_S_dt += kM * Sti * B1;            
        
        // LTM produced from STM conversion in the presence of tox
        dLTM_dt += kL * Tox * STM_S;
        dSTM_S_dt -= kL * Tox * STM_S;

        // Antidote production catalyzed by Sti and LTM            
        dAnt_dt += kA * LTM * Sti * B1;
    }

    // Optional cap on the antidote production rate
    // float MAX_ANT_RATE = 0.02;
    // dAnt_dt = min(dAnt_dt, MAX_ANT_RATE);


    // *** Antidote reacts with Toxin
    // ==============================
    dTox_dt -= k_NT * Ant * Tox;
    dAnt_dt -= k_NT * Ant * Tox;

    // *** Antidote detrimental to B
    // =============================
    dB1_dt -= k_NB * Ant * B1 ;
    dAnt_dt -= k_NB * Ant * B1 ;

    // We are done with the learning network : we write the dest
    dest[xy(id)].Tox = max(0, Tox + dt * dTox_dt);
    dest[xy(id)].Ant = max(0, Ant + dt * dAnt_dt);
    dest[xy(id)].Sti = max(0, Sti + dt * dSti_dt);
    dest[xy(id)].STM_S = max(0, STM_S + dt * dSTM_S_dt);
    dest[xy(id)].LTM = max(0, LTM + dt * dLTM_dt);


    // EXTERNAL BOLUSES
    // Periodically delivered : Sti then Tox
    // It starts with a "quiet time" without any delivery
    // Sti bolus concentration is constant
    // Tox bolus concentration is modulated by the epsilon environment factor, that varies slowly

    if(time > Bolus_quiet_time)
    {
        // 1) Stimulus bolus
        if(fmod(time, Bolus_period) < Bolus_duration)        
        {
            if(!Bolus_AtBoundaries || (Bolus_AtBoundaries && IsAtBoundary(id))) dest[xy(id)].Sti = Bolus_concentration;
        }


        // 2) Toxin bolus
        // First Compute environment parameter epsilon, the factor between Sti and Tox bolus concentration
        float epsilon_environment = 0;
        if((time-Bolus_quiet_time) < Bolus_shift)
        {
            epsilon_environment = ((time-Bolus_quiet_time)/Bolus_shift);
        }
        else if ((time-Bolus_quiet_time) < 2 * Bolus_shift)
        {
            epsilon_environment = (2 - (time-Bolus_quiet_time)/Bolus_shift);
        }

        // Toxin bolus delivery
        if(fmod(time - Bolus_delta, Bolus_period) < Bolus_duration)        
        {
            if(!Bolus_AtBoundaries || (Bolus_AtBoundaries && IsAtBoundary(id))) dest[xy(id)].Tox = epsilon_environment * Bolus_concentration;
        }
    }

#endif



    // ========================
    // 4) UPDATE GS DESTINATION
    // ========================
    
    dest[xy(id)].A1 = max(0, A1 + dt * dA1_dt);
    dest[xy(id)].B1 = max(0, B1 + dt * dB1_dt);
    // Boundary condition : no B on the boundary
    if(!perX & (id.x == 0 || id.x == RESX - 1 )) dest[xy(id)].B1 = 0;
    if(!perY & (id.y == 0 || id.y == RESY - 1 )) dest[xy(id)].B1 = 0;
}