// *** ShaderAux.hlsl *** 
//
// Useful auxiliary functions for compute shader calculation, inclusing accessors, color palette & Perlin noise



// *** Advection
// Computes the upflow advection of C in vector field v, using neighbour values Cmx, Cpx, Cmy, Cpy
float Advection(float2 v, float C, float Cmx, float Cpx, float Cmy, float Cpy)
{
    float adv = 0;
    if (v.x > 0) adv += v.x * (C - Cmx) / dx;
    else adv += v.x * (Cpx - C) / dx;
    if (v.y > 0) adv += v.y * (C - Cmy) / dx;
    else adv += v.y * (Cpy - C) / dx;
    return adv;
}


// ==============
// INDICES ACCESS
// ==============

// WARNING : Assumes RESX, RESY, perX, perY have been defined.
// and access is column-major using : id.x * RESY + id.y

bool IsAtBoundary(uint3 id)
{
   return (id.x == 0 || id.x == RESX - 1 || id.y == 0 || id.y == RESY - 1);
}


int xy(uint3 id) {
	return id.x * RESY + id.y;
}

int px(uint3 id) {

	uint x = id.x + 1;
	uint y = id.y;

	if (x == RESX) {
		if (perX) x = 0;
		else x = RESX - 1;
	}

	return x * RESY + y;
}

int mx(uint3 id) {

	int x = id.x - 1;
	int y = id.y;

	if (x == -1) {
		if (perX) x = RESX - 1;
		else x = 0;
	}

	return x * RESY + y;
}

int py(uint3 id) {

	uint x = id.x;
	uint y = id.y + 1;

	if (y == RESY) {
		if (perY) y = 0;
		else y = RESY - 1;
	}

	return x * RESY + y;
}

int my(uint3 id) {

	int x = id.x;
	int y = id.y - 1;

	if (y == -1) {
		if (perY) y = RESY - 1;
		else y = 0;
	}

	return x * RESY + y;
}



// ==============
// COLOR PALETTES
// ==============

float4 colorT5(float T, float TMIN, float TMAX) {

	float x = min(0.999, (T - TMIN) / (TMAX - TMIN));
	float y = (5 * x) % 1;

	if (x < 0.2) return float4(0, 0, y, 1);

	if (x < 0.4) return float4(0, y, 1, 1);

	if (x < 0.6) return float4(0, 1, 1 - y, 1);

	if (x < 0.8) return float4(y, 1, 0, 1);

	return float4(1, 1 - y, 0, 1);
}

float4 colorTerrainOLD(float h, float HMAX)
{
    float z = h / HMAX;
    float x = 3*fmod(z, 0.33);

    if (z < 0)
    {
        return float4(0.0, 0.0, 1.0, 1.0);
    }
    if (z < 0.33)
    {
        float4 col1 = float4(0.133, 0.5, 0.133, 1);
        float4 col2 = float4(1.0, 1.0, 0.6, 1);
        return (1 - x) * col1 + x * col2;
    }
    if (z < 0.66)
    {
        float4 col1 = float4(1.0, 1.0, 0.6, 1.0);
        float4 col2 = float4(0.5, 0.36, 0.33, 1.0);
        return (1 - x) * col1 + x * col2;
    }
    if (z < 0.99)
    {
        float4 col1 = float4(0.5, 0.36, 0.33, 1.0);
        float4 col2 = float4(1.0, 1.0, 1.0, 1.0);
        return (1 - x) * col1 + x * col2;
    }
    return float4(1.0, 1.0, 1.0,1.0);
}

float4 colorTerrain(float h, float HMAX)
{
    float z = h / HMAX;
    float x = 2*fmod(z, 0.5);

    if (z < 0)
    {
        return float4(0.0, 0.0, 1.0, 1.0);
    }
    if (z < 0.5)
    {
        float4 col1 = float4(0.133, 0.5, 0.133, 1);
        float4 col2 = float4(1.0, 1.0, 0.6, 1);
        return (1 - x) * col1 + x * col2;
    }
    if (z < 0.99)
    {
        float4 col1 = float4(1.0, 1.0, 0.6, 1.0);
        float4 col2 = float4(0.5, 0.36, 0.33, 1.0);
        return (1 - x) * col1 + x * col2;
    }
    return float4(0.5, 0.36, 0.33, 1.0);
}

float4 colorLerp(float f, float FMIN, float FMAX, float4 colmin, float4 colmax)
{
    float x = clamp((f - FMIN) / (FMAX - FMIN), 0, 1);
    return (1 - x) * colmin + x * colmax;
}

float4 color2Lerp(float f, float FMIN, float FMAX, float4 colmin, float4 colmid, float4 colmax)
{
    float x = clamp((f - FMIN) / (FMAX - FMIN), 0, 0.9999);
	float y = (2 * x) % 1;
	
	if(x<0.5) return (1 - y) * colmin + y * colmid;
	
	return (1 - y) * colmid + y * colmax;    
}


float4 PureHue(float x)
{	
	float u = 0.1666666666;
	float y = (6 * clamp(x,0,0.99999)) % 1;
	if(x<u) 	   return float4(1, 		y, 		0, 		1); 	
	else if(x<2*u) return float4(1- y, 	1, 		0, 		1);
	else if(x<3*u) return float4(0, 		1, 		y, 		1);
	else if(x<4*u) return float4(0, 		1 - y, 	1, 		1);
	else if(x<5*u) return float4(y, 		0, 		1, 		1);
	else 		   return float4(1, 		0, 		1 - y, 	1);
}


float4 colorHueSat(float2 v, float RMAX, float4 colzero)
{	
	float r = sqrt(v.x * v.x + v.y * v.y);
	float angle_norm = 0.5 + 0.5 * (atan2(v.y,v.x) / PI);	
	angle_norm = (angle_norm + 0.25) % 1;	// shift so that Red is at the North direction
	
	float4 colpure = PureHue(angle_norm);

	float x = clamp(r/RMAX, 0, 1);

	return (1-x) * colzero + x * colpure;
}


// Note : the mid-function 'return' will cause a warning : use of potentially uninitialized variable in Unity
float4 colorInferno(float v, float VMIN, float VMAX)
{
    float4 c0 = (1 / 255.0) * float4(0, 0, 4, 255);
    float4 c1 = (1 / 255.0) * float4(22, 11, 57, 255);
    float4 c2 = (1 / 255.0) * float4(66, 10, 104, 255);
    float4 c3 = (1 / 255.0) * float4(106, 23, 110, 255);
    float4 c4 = (1 / 255.0) * float4(147, 38, 103, 255);
    float4 c5 = (1 / 255.0) * float4(186, 54, 85, 255);
    float4 c6 = (1 / 255.0) * float4(221, 81, 58, 255);
    float4 c7 = (1 / 255.0) * float4(243, 118, 27, 255);
    float4 c8 = (1 / 255.0) * float4(252, 165, 10, 255);
    float4 c9 = (1 / 255.0) * float4(246, 215, 70, 255);
    float4 c10 = (1 / 255.0) * float4(246, 215, 70, 255);
    
    float x = min(0.999, (v - VMIN) / (VMAX - VMIN));
    float y = (10 * x) % 1;

    if (x < 0.1)               
        return (1 - y) * c0 + y * c1;
    if (x < 0.2)               
        return (1 - y) * c1 + y * c2;
    if (x < 0.3)               
        return (1 - y) * c2 + y * c3;
    if (x < 0.4)               
        return (1 - y) * c3 + y * c4;
    if (x < 0.5)               
        return (1 - y) * c4 + y * c5;
    if (x < 0.6)               
        return (1 - y) * c5 + y * c6;
    if (x < 0.7)               
        return (1 - y) * c6 + y * c7;
    if (x < 0.8)               
        return (1 - y) * c7 + y * c8;
    if (x < 0.9)               
        return (1 - y) * c8 + y * c9;
                   
    return (1 - y) * c9 + y * c10;
    
}


// ============
// PERLIN NOISE
// ============

// Inspired by https://www.shadertoy.com/view/XscGzl

float map(float value, float old_lo, float old_hi, float new_lo, float new_hi)
{
	float old_range = old_hi - old_lo;
    if (old_range == 0.0) {
	    return new_lo; 
	} else {
	    float new_range = new_hi - new_lo;  
	    return (((value - old_lo) * new_range) / old_range) + new_lo;
	}
}

float hash(float x)
{
	return frac(sin(x) * 43758.5453123);
}

// A pseudorandom gradient
float3 gradient(float3 cell)
{
    // To get spatial periodic noise, NCELL should be the same in dynperlin
    float NCELL = 16;
    cell.x = cell.x % NCELL;
    cell.y = cell.y % NCELL;

	float h_i = hash(cell.x);
	float h_j = hash(cell.y + pow(h_i, 3.0));
	float h_k = hash(cell.z + pow(h_j, 5.0));
    float ii = map(frac(h_i + h_j + h_k), 0.0, 1.0, -1.0, 1.0);
    float jj = map(frac(h_j + h_k), 0.0, 1.0, -1.0, 1.0);
	float kk = map(h_k, 0.0, 1.0, -1.0, 1.0);
    float l = length(float3(ii, jj, kk));
    if(l>0) return normalize(float3(ii, jj, kk));
    else return float3(1,0,0);
}


float fade(float t)
{
   	float t3 = t * t * t;
    float t4 = t3 * t;
    float t5 = t4 * t;
    return (6.0 * t5) - (15.0 * t4) + (10.0 * t3);        
}    


float PerlinNoise(float3 coord)
{
    float3 cell = floor(coord);
    float3 unit = frac(coord);
   
    float3 unit_000 = unit;
    float3 unit_100 = unit - float3(1.0, 0.0, 0.0);
    float3 unit_001 = unit - float3(0.0, 0.0, 1.0);
    float3 unit_101 = unit - float3(1.0, 0.0, 1.0);
    float3 unit_010 = unit - float3(0.0, 1.0, 0.0);
    float3 unit_110 = unit - float3(1.0, 1.0, 0.0);
    float3 unit_011 = unit - float3(0.0, 1.0, 1.0);
    float3 unit_111 = unit - float3(1.0, 1.0, 1.0);

    float3 c_000 = cell;
    float3 c_100 = cell + float3(1.0, 0.0, 0.0);
    float3 c_001 = cell + float3(0.0, 0.0, 1.0);
    float3 c_101 = cell + float3(1.0, 0.0, 1.0);
    float3 c_010 = cell + float3(0.0, 1.0, 0.0);
    float3 c_110 = cell + float3(1.0, 1.0, 0.0);
    float3 c_011 = cell + float3(0.0, 1.0, 1.0);
    float3 c_111 = cell + float3(1.0, 1.0, 1.0);

    float wx = fade(unit.x);
    float wy = fade(unit.y);
    float wz = fade(unit.z);
 
    float x000 = dot(gradient(c_000), unit_000);
	float x100 = dot(gradient(c_100), unit_100);
	float x001 = dot(gradient(c_001), unit_001);
	float x101 = dot(gradient(c_101), unit_101);
	float x010 = dot(gradient(c_010), unit_010);
	float x110 = dot(gradient(c_110), unit_110);
	float x011 = dot(gradient(c_011), unit_011);
	float x111 = dot(gradient(c_111), unit_111);
   
    // (0,0,0) - (1,0,0)
    // (0,0,1) - (1,0,1)
    // (0,1,0) - (1,1,0)
    // (0,1,1) - (1,1,1)
    float y0 = lerp(x000, x100, wx);
    float y1 = lerp(x001, x101, wx);
    float y2 = lerp(x010, x110, wx);
    float y3 = lerp(x011, x111, wx);
    
	float z0 = lerp(y0, y2, wy);
    float z1 = lerp(y1, y3, wy);
    
    return lerp(z0, z1, wz);
}	



