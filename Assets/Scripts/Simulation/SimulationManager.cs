// Below are preprocessor directives to chose the desired simulation features. 
// They should be identical in SimulationEngine.compute and SimulationManager.cs

#define FLOW
#define HEAT
#define GS2
#define LEARNING


// Preproc directive should be used only when they save time :
// - for Point in the ComputeShader and equations using those variables
// - for the CellData structure and everything that uses it by accessing a given variable (like initial condition currentState[...].<something> )


// *********************************************************
// *** Instructions for adding a new variable / equation ***
// *********************************************************
// This is not straightforward :-)
// It should be done in a consistent way in both SimulationManager.cs (SM) and SimulationEngine.compute (SE)
// If needed, define a corresponding preprocessor directive in SM and SE
//
// 1) New variable definition
// Add it to CellData in SM
// Add it to Point in SE (same order as CellData in SM)
// 
// 2) Texture definition (if you want to visualize that variable as a texture)
// Define the corresponding RWTexture2D in SE
// Define the associated texture in enum TEX in ImageTexture
// Set the texture link in SetTextures() in SM
// Populate the texture in the dump kernel in SE with the relevant variable
//
// 3) Evolution
// Define the relevant evolution equations in the Update() part of the compute shader kernel in SE
// For efficiency it is better do it such that dest[] is not read, and is written only once.
// When adding an equation, if parameters are needed they need : 
// - to be added in the inspector
// - corresponding parameter should be added to cbuffer GLOB
// - Parameters should be then passed to the cbuffer GLOB in SetGlobalParameters()
// Optionnaly expose the parameter value in Parameters()
//
// Optionnaly :
// If you want interactivity on this variable, add it to ApplyMouseInput in SE, the number should be the same as the enum TEX
// *********************************************************


using System.IO;
using System.Globalization;
using UnityEngine;
using UnityEngine.UI;
using System.Collections;
using System.Collections.Generic;
using UnityEngine.SceneManagement;
using UnityEngine.Serialization;
using System.Runtime.InteropServices;   // For Marshal.SizeOf

public class SimulationManager : MonoBehaviour
{
    static SimulationManager S;

    // ==================
    // INSPECTOR SETTINGS
    // ==================

    public bool auto_play = false;

    [Header("Space / time discretization")]
    public int RESX = 256;          // should be an integer multiple of 32
    public int RESY = 128;          // should be an integer multiple of 32
    public float DX = 1;            // Cell size
    public bool PERX = false;       // Spatial Periodic boundary conditions ?
    public bool PERY = false;       // 
    public float DT = 1;            // Depends on diffusion or reaction constant. Use 1 for pure gray scott, 0.01f for fluid simulation or if there is diffusion boundary condition with high diffusion coefficient
    public int STEPS = 250;         // Number of simulation steps per frame

    [Header("Gray-Scott")]
    public float CHI_A = 0.2f;          // Diffusion coefficient
    public float CHI_B = 0.1f;
    public float GS_F = 0.03f;          // Gray scott reaction parameters F and R
    public float GS_R = 0.061f;
    public float k0 = 1f;               // Autocatalysis reaction constant
    public float initialProportionOfPixelsB1 = 0.25f;


    [Header("GS2")]
    public float initialProportionOfPixelsB2 = 0;


    [Header("GS Temperature")]
    public float Ea = 1.7f;             // Activation energy
    public float DeltaH1 = 0* 0.025f;   // Enthaply of GS chemical reactions for Gray-Scott systems 1 and 2 
    public float DeltaH2 = 0*-0.02f;    // negative means exothermic


    public enum AntidoteProduction {None, Direct, Preemptive, Associative};

    [Header("Learning systems")]    
    public AntidoteProduction antidoteProduction = AntidoteProduction.None;
    
    public bool Bolus_AtBoundaries = true;      // Bolus delivered at boundaries ? (or uniformly)
    public float Bolus_quiet_time = 2000;       // How much quiet time before toxin starts (so that GS system can grow)       
    public float Bolus_period = 1000;           // Period of boluses
    public float Bolus_duration = 10;           
    public float Bolus_concentration = 1;
    public float Bolus_delta = 100;             // Time difference between stimulus and toxin
    public float Bolus_shift = 100;             // Characteristic time of the bolus concentration evolution slope

    // Diffusion and decay rates
    public float Sti_Diffusion = 2.5f;
    public float Sti_decay = 0.01f;
    public float Tox_Diffusion = 2.5f;
    public float Tox_decay = 0.015f;
    public float Ant_Diffusion = 0.6f;
    public float Ant_decay = 0.004304f;
    public float STM_Diffusion = 0.6f;
    public float STM_decay = 0.00717f;
    public float LTM_Diffusion = 0.6f;
    public float LTM_decay = 0.000083f;

    // Reaction constants
    public float k_AntTox = 50;        // Antidote reacts with Toxin
    public float k_ToxB = 1;                                        // Toxin reaction with B
    public float k_AntB = 0.01f;                                    // Antidote reaction with B
    public float k_Reactive = 0;       // Reactive production of Antidote from Toxin
    public float k_Preemptive;         // Preemptive production of Antidote from Stimulus
    public float k_Prod_STM = 4;       // STM production from stimulus
    public float k_Prod_LTM = 1;       // LTM production from STM conversion in the presence of toxin
    public float k_Output = 3.757f;    // Antidote output in the presence of LTM and Stimulus

    
    [Header("Thermal simulation")]
    public bool simu_Temperature = false;   // Simulate temperature ?
    public float T0 = 1.5f;                 // Average temperature
    public float DeltaT = 0*0.5f;           // Top/bottom temperature difference
    public float LAMBDA = 0.1f;             // Thermal conductivity. Max is around 25 for dx=1 dt=0.01
    public float Cv = 1;                    // Heat capacity

    // pulsation and amplitude of temporal Temperature variation at boundary
    public float OmegaT = 0 * 1f / 10000f;
    public float AmpliT = 0 * 0.4f;
    public float BoundaryTransferCoefficient = 1;
    public bool initialTemperatureGradient = true;


    [Header("Fluid simulation")]
    public bool simu_v = false;             // Fluid velocity ?
    public float VX0 = 0;                   // Stationnary flow
    public float VY0 = 0;
    public float RHO0 = 1;
    public float Rs = 10;
    public float G = 0.05f;           // we want G.(RESY*DX) of the order of Rs*T0 (hydrostatic vs thermal pressure)
    public float mu = 0.01f;          // max is around 12 for dx=1 dt=0.01
    


    // ===============================================================================================================

    // =========
    // CELL DATA
    // =========
    
    // Contain all variables to be passed to the compute shader
    // WARNING : Order here should match the struct Point in the compute shader. 

    struct CellData
    {
        public float A1;
        public float B1;

#if LEARNING
        public float St;
        public float Tx;
        public float An;
        public float STMem_S;
        public float LTMem;
#endif

#if HEAT
        public float T;
        public float chem;
#endif

#if FLOW
        public float rho;
        public Vector2 v;
#endif

#if GS2
        public float A2;
        public float B2;
#endif

    }


// ===============================================================================================================

    // TEXTURES FOR VISUALISATION
    
    // This list should match :
    // - the RWTexture2D declared in the compute shader
    // - the numeric ids of texture used in the interactivity part of the compute shader.
    // - the shader.SetTexture list in SimulationManager
    public enum TEX {Temperature, Density, Velocity, GS1, GS2, Toxin, Antidote, Stimulus, ShortTermMemory_S, LongTermMemory, COUNT};

// ===============================================================================================================

    static void SetGlobalParameters()
    {
        shader.SetFloat("time", 0);

        // General parameters
        shader.SetBool("simu_T", S.simu_Temperature);
        shader.SetBool("simu_v", S.simu_v);
        shader.SetInt("memory_mechanism", (int) S.antidoteProduction);

        // Grid
        shader.SetInt("RESX", S.RESX);
        shader.SetInt("RESY", S.RESY);
        shader.SetFloat("dx", S.DX);
        shader.SetFloat("dt", S.DT);
        shader.SetBool("perX", S.PERX);
        shader.SetBool("perY", S.PERY);

        // Thermal and fluid model
        shader.SetFloat("T0", S.T0);
        shader.SetFloat("DeltaT", S.DeltaT);        
        shader.SetFloat("LAMBDA", S.LAMBDA);
        shader.SetFloat("Cv", S.Cv);
        shader.SetFloat("BOUNDARY_HEAT_TRANSFER_COEFFICIENT", S.BoundaryTransferCoefficient);
        shader.SetFloat("OmegaT", S.OmegaT);
        shader.SetFloat("AmpliT", S.AmpliT);

        shader.SetFloat("VX0", S.VX0);
        shader.SetFloat("VY0", S.VY0);
        shader.SetFloat("RHO0", S.RHO0);
        shader.SetFloat("Rs", S.Rs);
        shader.SetFloat("G", S.G);
        shader.SetFloat("mu", S.mu);


        // Gray Scott
        shader.SetFloat("A_diffusion", S.CHI_A);
        shader.SetFloat("B_diffusion", S.CHI_B);
        shader.SetFloat("GS_F", S.GS_F);
        shader.SetFloat("GS_R", S.GS_R);
        shader.SetFloat("Ea", S.Ea);
        shader.SetFloat("k0", S.k0);

        shader.SetFloat("DeltaH1", S.DeltaH1);
        shader.SetFloat("DeltaH2", S.DeltaH2);


        // Learning
        
        shader.SetBool("Bolus_AtBoundaries", S.Bolus_AtBoundaries);
        shader.SetFloat("Bolus_quiet_time", S.Bolus_quiet_time);
        shader.SetFloat("Bolus_period", S.Bolus_period);
        shader.SetFloat("Bolus_concentration", S.Bolus_concentration);
        shader.SetFloat("Bolus_duration", S.Bolus_duration);
        shader.SetFloat("Bolus_delta", S.Bolus_delta);
        shader.SetFloat("Bolus_shift", S.Bolus_shift);

        shader.SetFloat("Sti_diffusion", S.Sti_Diffusion);
        shader.SetFloat("Sti_decay", S.Sti_decay);        

        shader.SetFloat("Tox_diffusion", S.Tox_Diffusion);
        shader.SetFloat("Tox_decay", S.Tox_decay);
        shader.SetFloat("Ant_diffusion", S.Ant_Diffusion);
        shader.SetFloat("Ant_decay", S.Ant_decay);

        shader.SetFloat("STM_diffusion", S.STM_Diffusion);
        shader.SetFloat("STM_decay", S.STM_decay);
        shader.SetFloat("LTM_diffusion", S.LTM_Diffusion);
        shader.SetFloat("LTM_decay", S.LTM_decay);

        shader.SetFloat("kD", S.k_Reactive);
        shader.SetFloat("k_TB", S.k_ToxB);
        shader.SetFloat("k_NB", S.k_AntB);

        shader.SetFloat("kM", S.k_Prod_STM);
        shader.SetFloat("kL", S.k_Prod_LTM);        
        shader.SetFloat("kA", S.k_Output);
        shader.SetFloat("kP", S.k_Preemptive);
        shader.SetFloat("k_NT", S.k_AntTox);

    }

    // =============================================================================

    public static string Parameters()
    {
        string sf_glob = $"Mouse input size={S.input_size}, strength={S.input_strength} // {output_filename} : RES={S.RESX}x{S.RESY} DX={S.DX} DT={S.DT} STEPS={S.STEPS} PER_X={S.PERX} PER_Y={S.PERY}";
        string sf_gs = $"CHI_A={S.CHI_A}, CHI_B={S.CHI_B}, GS_F={S.GS_F}, GS_R={S.GS_R}, k0={S.k0}, p_init_B1={S.initialProportionOfPixelsB1} ";
        string sf_thermo = $"Ea={S.Ea}, DeltaH1={S.DeltaH1} ";
        string sf_gs2 = $"DeltaH2={S.DeltaH2}, p_init_B2={S.initialProportionOfPixelsB2} ";
        string s_ht = $"HeatTransfer={S.simu_Temperature} T0={S.T0}, DeltaT={S.DeltaT}, LAMBDA={S.LAMBDA}, Cv={S.Cv}, h_bnd={S.BoundaryTransferCoefficient}, OmegaT={S.OmegaT}, AmpliT={S.AmpliT} ";
        string s_flow = $"NavierStokes={S.simu_v} RHO0={S.RHO0}, Rs={S.Rs}, G={S.G}, mu={S.mu}, Ra ={Rayleigh():E2}, Gr ={Grashof():E2}, Pr ={Prandtl():E2}";
        string s_tox = $"Boluses quiet_time={S.Bolus_quiet_time}, period={S.Bolus_period}, time_shift={S.Bolus_delta}, concentration={S.Bolus_concentration}, Sti_decay={S.Sti_decay}, Sti_diff={S.Sti_Diffusion}, Tox_decay={S.Tox_decay}, Tox_diff={S.Tox_Diffusion}";
        string s_learn = $"Mode={S.antidoteProduction} : Ant_diff={S.Ant_Diffusion}, Ant_decay={S.Ant_decay}, STM_diff={S.STM_Diffusion}, STM_decay={S.STM_decay}, LTM_diff={S.LTM_Diffusion}, LTM_decay={S.LTM_decay}";
        string s_react = $"Reaction constants : k_Reactive ={S.k_Reactive}, k_Preemptive ={S.k_Preemptive}, k_STM ={S.k_Prod_STM}, k_LTM ={S.k_Prod_LTM}, k_Output ={S.k_Output}, k_AntTox ={S.k_AntTox}, k_ToxB ={S.k_ToxB}, k_AntB ={S.k_AntB}";

        string sf = "";
        sf += $"** GLOBALS ** {sf_glob}\n";
        sf += $"** GRAY-SCOTT ** {sf_gs}\n";
        sf += $"** HEAT ** {s_ht}, {sf_thermo}\n";
        sf += $"** GS 2 ** {sf_gs2}\n";
        sf += $"** LEARNING ** {s_tox}  /  {s_learn}  /  {s_react}\n";
        sf += $"** FLOW ** {s_flow}\n";

        return sf;
    }

    public static float Rayleigh()
    {
        float Lc = S.RESY * S.DX;
        float nu = S.mu / S.RHO0;
        float beta = 1 / S.T0;
        float alpha = S.LAMBDA / (S.RHO0 * S.Cv);
        return S.G * beta / (nu * alpha) * S.DeltaT * (Lc * Lc * Lc);
    }

    public static float Grashof()
    {
        float beta = 1 / S.T0;
        float Lc = S.RESY * S.DX;
        return S.G * beta * S.DeltaT * (Lc * Lc * Lc) * S.RHO0 * S.RHO0 / (S.mu * S.mu);
    }

    public static float Prandtl()
    {
        return S.mu * S.Cv / S.LAMBDA;
    }


// ===============================================================================================================

    // Helpers
    bool IsInBounds(int i, int j)
    {
        return (i >= 0) && (i < RESX) && (j >= 0) && (j < RESY);
    }

    public static RenderTexture GetTexture(TEX TEX_ID)
    {
        return Textures[(int) TEX_ID];
    }


    void OnValidate()
    {
        Debug.Log("Updated parameters");
        if(Time.time > 0 && shader != null)
        {
            SetGlobalParameters();
        }
    }


    void OnApplicationQuit()
    {
        CreateOutput();        
        
        buffer[0].Release();
        buffer[1].Release();
    }


    public static void SetInteractiveTexture(TEX tex)
    {
        S.InteractiveTexture = tex;
        RawImage rawImage = ImageTexture.RawImages[(int) S.InteractiveTexture];
        rawImage.rectTransform.GetWorldCorners(S.InteractiveTextureWorldCorners);            
    }

// ===============================================================================================================

    // COMPUTE SHADER CIRCUITRY : 
    // shader, kernel ID and compute buffer to store and pass data
    static ComputeShader shader;
    int NTGX, NTGY;
    static int initKernel;
    static int mainKernel;
    static int dumpKernel1;
    static int dumpKernel2;
    static ComputeBuffer[] buffer = new ComputeBuffer[2];
    static CellData[] currentState;
    // int ids for following and swapping the buffers
    static int current;
    static int next;    


    // TEXTURE PROCESSING
    static RenderTexture[] Textures;
    const int NTEX = (int) TEX.COUNT;


    // SIMULATION & INTERACTIVITY    
    public static float time;
    Vector3[] InteractiveTextureWorldCorners = new Vector3[4];
    TEX InteractiveTexture;
    int input_size = 5;
    float input_strength = 1;


    // RECORDING & STATISTICS
    public static float Bmean = 0;
    public static int Nspots = 0;
    
    static string output_filename;
    const float RECORD_INTERVAL = 10;
    static string times_list;
    static string Bmean_list;
    static string Nspots_list;

    // public static float MaxConcentration = 0;
    // static string Ant_list;
    // static string Tox_list;


// ===============================================================================================================

    private void Awake()
    {
        S = this;

        // Call the shader, find the kernels
        shader = Resources.Load<ComputeShader>("Shaders/SimulationEngine");
        initKernel = shader.FindKernel("CSInit");
        mainKernel = shader.FindKernel("CSMain");
        dumpKernel1 = shader.FindKernel("CSDump1");
        dumpKernel2 = shader.FindKernel("CSDump2");

        // Create the array of texture that will receive the data
        Textures = new RenderTexture[NTEX];
        for(int i = 0; i < Textures.Length; i++)
        {
            Textures[i] = new RenderTexture(RESX, RESY, 0);
            Textures[i].enableRandomWrite = true;
            Textures[i].filterMode = FilterMode.Point;
            Textures[i].wrapMode = TextureWrapMode.Clamp;
            Textures[i].Create();
        }

        // Number of thread groups
        NTGX = 1 + ((S.RESX - 1) / 32);
        NTGY = 1 + ((S.RESY - 1) / 32);     

        
    }

// ===============================================================================================================

    void Start()
    {

        UnityEngine.Random.InitState(1);

        // PREPARE RECORDING OUTPUT & STATISTICS
        output_filename = System.DateTime.Now.ToString("yyyy_MM_dd_HH_mm_ss") + " " + SceneManager.GetActiveScene().name + " " + antidoteProduction + ".csv";
        time = 0;
        times_list = "0";
        Bmean_list = "0";
        Nspots_list = "0";

        // Tox_list = "0";
        // Ant_list = "0";


        // INITIALISATION OF STATE & BUFFER
        currentState = new CellData[RESX * RESY];
        InitializeCurrentState();

        int STRIDE = System.Runtime.InteropServices.Marshal.SizeOf(typeof(CellData));
        buffer[0] = new ComputeBuffer(RESX * RESY, STRIDE);
        buffer[0].SetData(currentState);
        buffer[1] = new ComputeBuffer(RESX * RESY, STRIDE);
        buffer[1].SetData(currentState);

        // Initialise variables used for the swaping
        current = 0;
        next = 1;

        // Set Parameters in the global cbuffer of the compute shader
        SetGlobalParameters();
        Debug.Log(Parameters());


        // Call the Initiation kernel
        shader.SetBuffer(initKernel, "src", buffer[current]);
        shader.SetBuffer(initKernel, "dest", buffer[next]);
        shader.Dispatch(initKernel, NTGX, NTGY, 1);

        // Set Buffer for dumpKernel
        shader.SetBuffer(dumpKernel1, "src", buffer[current]);
        shader.SetBuffer(dumpKernel1, "dest", buffer[next]);
        shader.SetBuffer(dumpKernel2, "src", buffer[current]);
        shader.SetBuffer(dumpKernel2, "dest", buffer[next]);

        // Create link between compute shader texture and texture display
        SetTextures();

        shader.Dispatch(dumpKernel1, NTGX, NTGY, 1);
        shader.Dispatch(dumpKernel2, NTGX, NTGY, 1); 
    }


    // INITIALIZE HERE THE VARIABLES THAT HAVE TO BE INITIALIZED
    void InitializeCurrentState()
    {
        for (int i = 0; i < RESX; i++)
        {
            for (int j = 0; j < RESY; j++)
            {               

                // Initial conditions on the GS system
                currentState[i * RESY + j].A1 = 1;
                currentState[i * RESY + j].B1 = 0;
                float rnd1 = UnityEngine.Random.value;
                if (rnd1 < S.initialProportionOfPixelsB1)
                {
                    currentState[i * RESY + j].A1 = 0.5f;
                    currentState[i * RESY + j].B1 = 0.25f;
                }

//                 // Pearson initial conditions
//                 const float SIZE = 4;
//                 const float SIGMA = 1;
//                 if (i > RESX/2 - SIZE && i < RESX/2 + SIZE && j > RESY/2 - SIZE && j < RESY/2 + SIZE)
//                 {
//                    float z =  (i - RESX/2) * (i - RESX/2) + (j - RESY/2) * (j - RESY/2);
//                    float f = Mathf.Exp(-z/(SIGMA*SIGMA));
//                    currentState[i * RESY + j].A1 = 1-f; // 0.5f + 0.1f * UnityEngine.Random.value;
//                    currentState[i * RESY + j].B1 = f; // 0.25f + 0.1f * UnityEngine.Random.value;

                // SECOND GS SYSTEM
#if GS2
                // Bartlett initial conditions system 2
                currentState[i * RESY + j].A2 = 1;
                currentState[i * RESY + j].B2 = 0;
                float rnd2 = UnityEngine.Random.value;
                if (rnd2 < S.initialProportionOfPixelsB2)
                {
                    currentState[i * RESY + j].A2 = 0.5f;
                    currentState[i * RESY + j].B2 = 0.25f;
                }
#endif

                // Heritable trait
#if HERITABILITY
                if(simu_heritable)
                {
                    if(i > 6 * RESX / 8 && j > 6 * RESY / 8) currentState[i * RESY + j].H = 1f;
                    else currentState[i * RESY + j].H = 0f;
                    // // Checkerboard initial distribution of trait
                    // const int NSQUARES = 4;
                    // float x = (float) i / (float) RESX;
                    // float y = (float) j / (float) RESY;
                    // const float PI = 3.141592654f;
                    // int sx = (int) Mathf.Sign(Mathf.Sin(NSQUARES * PI * x));
                    // int sy = (int) Mathf.Sign(Mathf.Sin(NSQUARES * PI * y));
                    // if(sx*sy == 1) currentState[i * RESY + j].H = 1f;
                    // else currentState[i * RESY + j].H = 0f;
                }
#endif

                // THERMAL AND FLUID FLOW
#if HEAT
                currentState[i * RESY + j].T = T0 + (initialTemperatureGradient ? 0.5f * DeltaT - DeltaT * ((float)j / RESY) : 0);
#endif
#if FLOW
                currentState[i * RESY + j].rho = RHO0;
                currentState[i * RESY + j].v = Vector2.zero;
#endif

            }
        }


    }

    void SetTextures()
    {
        // It assumes the texture are properly named in the compute shader. This list should match the enum in ImageTexture
        shader.SetTexture(dumpKernel1, "Result_T", Textures[(int) TEX.Temperature]);
        shader.SetTexture(dumpKernel1, "Result_rho", Textures[(int)TEX.Density]);
        shader.SetTexture(dumpKernel1, "Result_v", Textures[(int)TEX.Velocity]);
        shader.SetTexture(dumpKernel1, "Result_GS1", Textures[(int)TEX.GS1]);
        shader.SetTexture(dumpKernel1, "Result_GS2", Textures[(int)TEX.GS2]);
        shader.SetTexture(dumpKernel2, "Result_Tox", Textures[(int)TEX.Toxin]);
        shader.SetTexture(dumpKernel2, "Result_Ant", Textures[(int)TEX.Antidote]);
        shader.SetTexture(dumpKernel2, "Result_Sti", Textures[(int)TEX.Stimulus]);
        shader.SetTexture(dumpKernel2, "Result_STM_S", Textures[(int)TEX.ShortTermMemory_S]);
        shader.SetTexture(dumpKernel2, "Result_LTM", Textures[(int)TEX.LongTermMemory]);
    }


// ================================================================================================================

    private void Update()
    {
        
        if (Input.GetKeyDown(KeyCode.Return)) S.auto_play = !S.auto_play;
        if (Input.GetKeyDown(KeyCode.Space) || S.auto_play)
        {
            MouseInteractivity();
            RunSimulation();
            ComputeStatistics();
        }
        const float CHANGE_STRENGTH_FACTOR = 1.2f;
        const float CHANGE_SIZE_FACTOR = 1.3f;
        if(Input.GetKeyDown(KeyCode.UpArrow)) input_strength *= CHANGE_STRENGTH_FACTOR;
        if(Input.GetKeyDown(KeyCode.DownArrow)) input_strength /= CHANGE_STRENGTH_FACTOR;
        if(Input.GetKeyDown(KeyCode.RightArrow)) input_size = (int) Mathf.Floor(input_size * CHANGE_SIZE_FACTOR);
        if(Input.GetKeyDown(KeyCode.LeftArrow)) input_size = (int) Mathf.Ceil(input_size / CHANGE_SIZE_FACTOR);
    }

    void MouseInteractivity()
    {
        if (Input.GetMouseButton(0) || Input.GetMouseButton(1))
        {
            // Get Mouse Position on texture
            float xmin = InteractiveTextureWorldCorners[0].x;
            float xmax = InteractiveTextureWorldCorners[2].x;
            float ymin = InteractiveTextureWorldCorners[0].y;
            float ymax = InteractiveTextureWorldCorners[2].y;
            
            float x_pos = (float)(Input.mousePosition.x);
            float y_pos= (float)(Input.mousePosition.y);

            float ix = (float)(x_pos - xmin) / (xmax - xmin);
            float iy = (float)(y_pos - ymin) / (ymax - ymin);

            int mouse_i = (int)Mathf.Floor(ix*RESX);
            int mouse_j = (int)Mathf.Floor(iy*RESY);
        
            // Tell the shader there is a mouse input, where and how
            shader.SetBool("inputMouse", true);
            shader.SetInt("inputSize", S.input_size);
            shader.SetFloat("inputStrength", S.input_strength);
            shader.SetInt("inputTexture", (int) S.InteractiveTexture);
            shader.SetInt("mouse_i", mouse_i);
            shader.SetInt("mouse_j", mouse_j);

            if (Input.GetMouseButton(0)) shader.SetBool("remove", false);
            else shader.SetBool("remove", true);
        }
        else
        {
            shader.SetBool("inputMouse", false);
        }
    }

    // ===============
    // SIMULATION LOOP
    // ===============

    void RunSimulation()
    {
        


        // ===========
        //  MAIN LOOP
        // ===========

        for (int s = 0; s < S.STEPS; s++)
        {
            shader.SetFloat("dt", S.DT);
            time += DT;

            // Set and dispatch
            shader.SetFloat("time", time);
            shader.SetBuffer(mainKernel, "src", buffer[current]);
            shader.SetBuffer(mainKernel, "dest", buffer[next]);
            shader.Dispatch(mainKernel, NTGX, NTGY, 1);

            // Swap buffers
            current = (current + 1) % 2;
            next = (next + 1) % 2;
        }        

        // Retrieve data from buffer
        buffer[current].GetData(currentState);

        // Dump results into textures
        shader.SetBuffer(dumpKernel1, "src", buffer[current]);
        shader.SetBuffer(dumpKernel1, "dest", buffer[next]);
        shader.Dispatch(dumpKernel1, NTGX, NTGY, 1);

        shader.SetBuffer(dumpKernel2, "src", buffer[current]);
        shader.SetBuffer(dumpKernel2, "dest", buffer[next]);
        shader.Dispatch(dumpKernel2, NTGX, NTGY, 1);
    }

// ================================================================================================================

    // ===================
    // STATISTICS & OUTPUT
    // ===================

    void CreateOutput()
    {
        string OUTPUT_PATH = "output";
        Directory.CreateDirectory(OUTPUT_PATH);         // created if does not exist   

        using (StreamWriter sw = File.CreateText("output" + Path.DirectorySeparatorChar.ToString() + output_filename))
        {
            sw.WriteLine(times_list);
            sw.WriteLine(Bmean_list);
            sw.WriteLine(Nspots_list);

            // sw.WriteLine(Tox_list);
            // sw.WriteLine(Ant_list);
        }
        Debug.Log("Output " + output_filename);
    }


    void ComputeStatistics()
    {
        Nspots = CountSpots();
        Bmean = ComputeBmean();

        // UpdateMaxConcentration();
        
        if (time % RECORD_INTERVAL < 1.5 * STEPS * DT)
        {
            times_list += "," + time.ToString();
            Bmean_list += "," + Bmean.ToString();
            Nspots_list += "," + Nspots.ToString();

            // Tox_list += "," + ComputeSpatialAverageTox().ToString();
            // Ant_list += "," + ComputeSpatialAverageAnt().ToString();
        }
    }


    // *** Count the number of Gray-Scott spots by using floodfill from pixels above a certain threshold for (B - A) 
    int CountSpots()
    {
        // First determine which are the pixels above the threshold
        const float THRESHOLD = -0.3f;
        int[,] to_segment = new int[RESX, RESY];
        for(int i = 0; i < RESX; i++)
        {
            for(int j = 0; j < RESY; j++)
            {
                if(currentState[i * RESY + j].B1 - currentState[i * RESY + j].A1 > THRESHOLD)
                {
                    to_segment[i, j] = 1;
                }
            }
        }


        // Then floodfill them one by one and removed the neighbours from the to_segment list
        int res = 0;
        for (int i = 0; i < RESX; i++)
        {
            for (int j = 0; j < RESY; j++)
            {
                if(to_segment[i, j] == 1)
                {
                    res++;
                    List<(int, int)> to_flood = new List<(int, int)>();
                    to_flood.Add((i, j));

                    while (to_flood.Count > 0)
                    {
                        
                        int k = to_flood[0].Item1;
                        int l = to_flood[0].Item2;

                        // Flood each neighbour if it exists, and add them to the list
                        if (IsInBounds(k - 1, l) && to_segment[k - 1, l] == 1)
                        {
                            to_segment[k - 1, l] = 0;
                            to_flood.Add((k - 1, l));
                        }
                        if (IsInBounds(k + 1, l) && to_segment[k + 1, l] == 1)
                        {
                            to_segment[k + 1, l] = 0;
                            to_flood.Add((k + 1, l));
                        }
                        if (IsInBounds(k, l - 1) && to_segment[k, l - 1] == 1)
                        {
                            to_segment[k, l - 1] = 0;
                            to_flood.Add((k, l - 1));
                        }
                        if (IsInBounds(k, l + 1) && to_segment[k, l + 1] == 1)
                        {
                            to_segment[k, l + 1] = 0;
                            to_flood.Add((k, l + 1));
                        }

                        to_flood.RemoveAt(0);
                    }                    
                }
            }
        }

        return res;
    }


    float ComputeBmean()
    {
        float res = 0;
        for (int i = 0; i < RESX; i++)
        {
            for (int j = 0; j < RESY; j++)
            {
                res += currentState[i * RESY + j].B1;
            }
        }
        return res / (RESX * RESY);        
    }


//     float ComputeSpatialAverageTox()
//     {
//         float res = 0;
// #if LEARNING
//         for (int i = 0; i < RESX; i++)
//         {
//             for (int j = 0; j < RESY; j++)
//             {                
//                 res += currentState[i * RESY + j].Tx;
//             }
//         }
//         res /= (RESX * RESY);
// #endif
//         return res;
//     }


//     float ComputeSpatialAverageAnt()
//     {
//         float res = 0;
// #if LEARNING
//         for (int i = 0; i < RESX; i++)
//         {
//             for (int j = 0; j < RESY; j++)
//             {
//                 res += currentState[i * RESY + j].An;
//             }
//         }
//         res /= (RESX * RESY);
// #endif
//         return res;
//     }


//     void UpdateMaxConcentration()
//     {
//         MaxConcentration = 0;
// #if LEARNING
//         for (int i = 0; i < RESX; i++)
//         {
//             for (int j = 0; j < RESY; j++)
//             {
//                 MaxConcentration = Mathf.Max(MaxConcentration, currentState[i * RESY + j].An);
//             }
//         }
// #endif
//     }





}
