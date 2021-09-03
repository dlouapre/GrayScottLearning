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


