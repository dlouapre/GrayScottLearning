int xy(uint3 id) {
	return id.x * RESY + id.y;
}

int px(uint3 id) {

	int x = id.x + 1;
	int y = id.y;

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

	int x = id.x;
	int y = id.y + 1;

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



