void InitiateRho(uint3 id)
{
	int xy = id.x * RESY + id.y;
	
	// Initial density according to hydrostatic law
	float zeta = G * RESY * dx / (Rs * T0);
    float rho_down = RHO0;
    if (zeta > 0)
    {
        rho_down = RHO0 * zeta / (1 - exp(-zeta));
    }        
	src[xy].rho = rho_down * exp(-zeta * id.y / RESY);
		
	// Initial density according to overpressure
    float RADIUS = 0.1 * RESX;
    float DELTA_RHO = 0.0;
    float d2 = (id.x - RESX / 2) * (id.x - RESX / 2) + (id.y - RESY / 2) * (id.y - RESY / 2);
    src[xy].rho += DELTA_RHO * exp(-d2 / (RADIUS * RADIUS));
}


// -------------------------------------------------------------------------------------------------------------------

void UpdateTemperature(uint3 id)
{
	//bool CENTRAL = true;
	
	float T = src[xy(id)].T;
	float Tmx = src[mx(id)].T;
	float Tpx = src[px(id)].T;
	float Tmy = src[my(id)].T;
	float Tpy = src[py(id)].T;

	float lap = (Tpx + Tmx + Tpy + Tmy - 4 * T) / (dx * dx);

	float2 v = src[xy(id)].v;
	float2 vpx = src[px(id)].v;
	float2 vmx = src[mx(id)].v;
	float2 vpy = src[py(id)].v;
	float2 vmy = src[my(id)].v;
	
	// chemical power from reaction
    float powerChem = src[xy(id)].chem;

    
	// Advection
	// Upwind scheme
	float gx, gy;
	if (v.x > 0)
		gx = (T - Tmx) / dx;
	else
		gx = (Tpx - T) / dx;

	if (v.y > 0)
		gy = (T - Tmy) / dx;
	else
		gy = (Tpy - T) / dx;

	float adv = v.x * gx + v.y * gy;
	
	// Divergence v is ignored
	float divv = 0;
    //if (CENTRAL)
    //    divv = 0.5 * (1 / dx) * (vpx[0] - vmx[0]) + 0.5 * (1 / dx) * (vpy[1] - vmy[1]);
    //else
    //    divv = (1 / dx) * (v[0] - vmx[0]) + (1 / dx) * (v[1] - vmy[1]);
	
	dest[xy(id)].T = T + (LAMBDA / (RHO0 * Cv)) * lap * dt - (adv + T * divv + powerChem) * dt;
	
	// Boundary conditions
	float Tup = T0 - 0.5 * DeltaT - AmpliT * sin(OmegaT * time);
	float Tdown = T0 + 0.5f * DeltaT - AmpliT * sin(OmegaT * time);

	if (id.y == 0) 
		dest[xy(id)].T += dt / (RHO0 * Cv * dx) * BOUNDARY_HEAT_TRANSFER_COEFFICIENT * (Tdown - T);
	if (id.y == (RESY - 1))
		dest[xy(id)].T += dt / (RHO0 * Cv * dx) * BOUNDARY_HEAT_TRANSFER_COEFFICIENT * (Tup - T);
}


// -------------------------------------------------------------------------------------------------------------------

void UpdateVelocity(uint3 id)
{
	bool CENTRAL = true;
	
	float2 v = src[xy(id)].v;
	float2 vpx = src[px(id)].v;
	float2 vmx = src[mx(id)].v;
	float2 vpy = src[py(id)].v;
	float2 vmy = src[my(id)].v;
	
	float rho = src[xy(id)].rho;
	float rhopx = src[px(id)].rho;
	float rhomx = src[mx(id)].rho;
	float rhopy = src[py(id)].rho;
	float rhomy = src[my(id)].rho;
	
	float T = src[xy(id)].T;
	float Tmx = src[mx(id)].T;
	float Tpx = src[px(id)].T;
	float Tmy = src[my(id)].T;
	float Tpy = src[py(id)].T;
	
	float P = Rs * rho * T;
	float Ppx = Rs * rhopx * Tpx;
	float Pmx = Rs * rhomx * Tmx;
	float Ppy = Rs * rhopy * Tpy;
	float Pmy = Rs * rhomy * Tmy;
	
	// 1) Advection
	float2 adv = float2(0, 0);
		
	// Upwind
	float2 gx;
	float2 gy;

	if (v.x > 0)
		gx = (v - vmx) / dx;
	else
		gx = (vpx - v) / dx;

	if (v.y > 0)
		gy = (v - vmy) / dx;
	else
		gy = (vpy - v) / dx;

	adv = v.x * gx + v.y * gy;
	
    /*
	// Central
	if (CENTRAL)
		adv = v[0] * 0.5 * (1 / dx) * (vpx - vmx) + v[1] * 0.5 * (1 / dx) * (vpy - vmy);
	else
		adv = v[0] * (1 / dx) * (v - vmx) + v[1] * (1 / dx) * (v - vmy);
	*/

	// 2) Pressure gradient
	float2 gP;
	if (CENTRAL)
		gP = float2(0.5 * (Ppx - Pmx) / dx, 0.5 * (Ppy - Pmy) / dx);
	else
		gP = float2((P - Pmx) / dx, (P - Pmy) / dx);


	// 3) Viscosity

	float2 lapV = (vpx + vmx + vpy + vmy - 4 * v) / (dx * dx);


	// 4) Gravity

	float2 gravity = float2(0, -G);

	// ADD UP ALL CONTRIBUTIONS
	float2 newv = v + dt * (-adv - (1 / rho) * gP + mu / rho * lapV + gravity);
	
	
	// CLAMP
	float VMAX = 0.5 * dx / dt;
	if (abs(newv.x) > VMAX)
		newv.x = sign(newv.x) * VMAX;
	if (abs(newv.y) > VMAX)
		newv.y = sign(newv.y) * VMAX;
	
	dest[xy(id)].v = newv;
	
	// No slip Boundary conditions
	if (!perY && id.y == 0)
		dest[xy(id)].v = float2(0, 0);
	if (!perY && id.y == RESY - 1)
		dest[xy(id)].v = float2(0, 0);
	if (!perX && id.x == 0)
        dest[xy(id)].v.x = 0;
		//dest[xy(id)].v = float2(0, 0);
	if (!perX && id.x == RESX - 1)
        dest[xy(id)].v.x = 0;
		//dest[xy(id)].v = float2(0, 0);
}

// -------------------------------------------------------------------------------------------------------------------

void UpdateRho(uint3 id)
{	
	bool CENTRAL = true;
	
	// Just a bit of diffusion to compensate the effect of central difference scheme
	float DIFF_RHO = 0.001;
	
	float2 v = src[xy(id)].v;
	float2 vpx = src[px(id)].v;
	float2 vmx = src[mx(id)].v;
	float2 vpy = src[py(id)].v;
	float2 vmy = src[my(id)].v;
	
	float rho = src[xy(id)].rho;
	float rhopx = src[px(id)].rho;
	float rhomx = src[mx(id)].rho;
	float rhopy = src[py(id)].rho;
	float rhomy = src[my(id)].rho;
		
	float2 rhov = rho * v;
	float2 rhov_px = rhopx * vpx;
	float2 rhov_mx = rhomx * vmx;
	float2 rhov_py = rhopy * vpy;
	float2 rhov_my = rhomy * vmy;

	float div_rhov;
	if (CENTRAL)
		div_rhov = (0.5f / dx) * (rhov_px[0] - rhov_mx[0]) + (0.5f / dx) * (rhov_py[1] - rhov_my[1]);
	else
		div_rhov = (1 / dx) * (rhov[0] - rhov_mx[0]) + (1 / dx) * (rhov[1] - rhov_my[1]);

	float lap = (rhomx + rhopx + rhomy + rhopy - 4 * rho) / (dx * dx);
	dest[xy(id)].rho = rho + dt * (-div_rhov + DIFF_RHO * lap);
}


// *******************************************************************************************************************



