using System.IO;
using System.Globalization;
using UnityEngine;
using UnityEngine.UI;
using System.Collections;
using System.Collections.Generic;
using UnityEngine.SceneManagement;

public class SimulationManager : MonoBehaviour
{
    static SimulationManager S;

    // ===============================================================================================================

    public enum AntidoteProduction {None, Direct, Preemptive, Associative};

    // INSPECTOR SETTINGS

    public bool auto_play = false;
    public float DT_FAST_FACTOR = 10;

    [Header("Simulation modules")]
    public bool simu_Temperature = true;
    public bool simu_v = false;
    public bool simu_GS_transport = true;
    public bool simu_GS_reactions = true;
    public AntidoteProduction antidoteProduction = AntidoteProduction.None;
    public bool compute_Bmean = false;

    [Header("Interactivity")]
    public ImageTexture.TEX InteractiveTexture;    
    public int InputSizeAdd = 4;
    public int InputSizeRemove = 20;

    [Header("Space / time discretization")]
    // Resolution & scale: should be an integer multiple of 32
    public int RESX = 256;
    public int RESY = 128;
    public float DX = 1;
    public float DT = 1;    //0.01f for fluid simulation;
    public int STEPS = 250;
    public bool PERX = false;
    public bool PERY = false;

    [Header("Thermal simulation")]
    public float T0 = 1.5f;
    public float DeltaT = 0*0.5f;
    public float LAMBDA = 0.1f;      // max is around 25 for dx=1 dt=0.01
    public float Cv = 1;
    // pulsation and amplitude of temporal Temperature variation at boundary
    public float OmegaT = 0 * 1f / 10000f;
    public float AmpliT = 0 * 0.4f;
    public float BoundaryTransferCoefficient = 1;
    public bool initialTemperatureGradient = true;


    [Header("Fluid simulation")]
    public float RHO0 = 1;
    public float Rs = 10;
    public float G = 0.05f;           // we want G.(RESY*DX) of the order of Rs*T0 (hydrostatic vs thermal pressure)
    public float mu = 0.01f;          // max is around 12 for dx=1 dt=0.01
    

    [Header("Gray-Scott systems")]
    public float CHI_A = 0.2f;
    public float CHI_B = 0.1f;
    public float GS_F = 0.03f;
    public float GS_R = 0.061f;
    public float Ea = 1.7f;    
    public float DeltaH1 = 0* 0.025f;    
    public float DeltaH2 = 0*-0.02f;   // negative means exothermic
    public float initialProportionOfPixelsB1 = 0.25f;
    public float initialProportionOfPixelsB2 = 0;


    [Header("Learning systems")]

    public bool Bolus_AtBoundaries = true;
    public bool Bolus_AllBoundaries = true;
    public float Bolus_quiet_time = 8500;
    public float Bolus_period = 1000;
    public float Bolus_duration = 10;
    public float Bolus_concentration = 1;
    public float Bolus_delta = 100;
    public float Bolus_shift = 100;

    public float Bolus_advection_x = 1;

    public float Sti_Diffusion = 2.5f;
    public float Pre_Diffusion = 2.5f;


    public float Sti_decay = 0.01f;
    public float Tox_Diffusion = 2.5f;
    public float Tox_decay = 0.015f;

    public float Ant_Diffusion = 0.6f;
    public float Ant_decay = 0.004304f;

    public float STM_Diffusion = 0.6f;
    public float STM_decay = 0.00717f;

    public float LTM_Diffusion = 0.6f;
    public float LTM_decay = 0.000083f;
  
    public float kR = 0;
    public float kT = 1;
    public float kS = 4;
    public float kL = 1;
    public float kD;
    public float kO = 3.757f;
    public float kA = 50;
    public float k_ToxB = 1;
    public float k_AntB = 0.01f;

    

    // ===============================================================================================================

    static bool running = false;
    public static float time;
    int NTGX, NTGY;
    public static float Bmean = 0;
    public static int Nspots = 0;
    public static float MaxConcentration = 0;


    // Evolution recording
    static string output_filename;
    const float RECORD_INTERVAL = 10;
    static string times_list;
    static string Bmean_list;
    static string Nspots_list;
    static string Ant_list;
    static string Tox_list;


    // Shader, kernel ID and compute buffer to store and pass data
    static ComputeShader shader;
    static int initKernel;
    static int mainKernel;
    static int dumpKernel1;
    static int dumpKernel2;
    static ComputeBuffer[] buffer = new ComputeBuffer[2];

    struct CellData
    {
        public float T;
        public float rho;
        public Vector2 v;
        public float x1;
        public float y1;
        public float chem;
        //public float x2;
        //public float y2;
        public float Tx;
        public float An;
        public float St;
        //public float Ct;
        public float STMem_S;
        //public float STMem_C;
        public float LTMem;
    }
    const int STRIDE = 12 * 4;

    static CellData[] currentState;

    // int ids for following and swapping the buffers
    static int current;
    static int next;    

    static RenderTexture[] Textures;
    const int NTEX = (int) ImageTexture.TEX.COUNT;

    Vector3[] InteractiveTextureWorldCorners = new Vector3[4];

    // ---------------------------------------------------------------------------------------------------------------

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

    // ---------------------------------------------------------------------------------------------------------------

    void Start()
    {

        UnityEngine.Random.InitState(1);

        output_filename = System.DateTime.Now.ToString("yyyy_MM_dd_HH_mm_ss") + " " + SceneManager.GetActiveScene().name + " " + antidoteProduction + ".csv";

        time = 0;
        times_list = "0";
        Bmean_list = "0";
        Nspots_list = "0";
        Tox_list = "0";
        Ant_list = "0";

        // Initial values
        // --------------

        currentState = new CellData[RESX * RESY];

        for (int i = 0; i < RESX; i++)
        {
            for (int j = 0; j < RESY; j++)
            {
                // Thermal & fluid
                currentState[i * RESY + j].rho = RHO0;
                currentState[i * RESY + j].T = T0 + (initialTemperatureGradient ? 0.5f * DeltaT - DeltaT * ((float) j/RESY) : 0);
                currentState[i * RESY + j].v = Vector2.zero;

                //// Pearson initial conditions
                //if (i > 0.45f * RESX && i < 0.55f * RESX && j > 0.45f * RESY && j < 0.55f * RESY)
                //{
                //    currentState[i * RESY + j].x1 = 0.5f + 0.1f * UnityEngine.Random.value;
                //    currentState[i * RESY + j].y1 = 0.25f + 0.1f * UnityEngine.Random.value;
                //}
                //else
                //{
                //    currentState[i * RESY + j].x1 = 1;
                //    currentState[i * RESY + j].y1 = 0;
                //}

                // Bartlett initial conditions
                currentState[i * RESY + j].x1 = 1;
                currentState[i * RESY + j].y1 = 0;
                float rnd1 = UnityEngine.Random.value;
                if (rnd1 < S.initialProportionOfPixelsB1)
                {
                    currentState[i * RESY + j].x1 = 0.5f;
                    currentState[i * RESY + j].y1 = 0.25f;
                }
                //// Bartlett initial conditions system 2
                //currentState[i * RESY + j].x2 = 1;
                //currentState[i * RESY + j].y2 = 0;
                //float rnd2 = UnityEngine.Random.value;
                //if (rnd2 < S.initialProportionOfPixelsB2)
                //{
                //    currentState[i * RESY + j].x2 = 0.5f;
                //    currentState[i * RESY + j].y2 = 0.25f;
                //}
            }
        }


        buffer[0] = new ComputeBuffer(RESX * RESY, STRIDE);
        buffer[0].SetData(currentState);
        buffer[1] = new ComputeBuffer(RESX * RESY, STRIDE);
        buffer[1].SetData(currentState);

        // Initialise variables used for the swaping
        current = 0;
        next = 1;

        SetGlobalParameters();


        // Call the Initiation kernel
        shader.SetBuffer(initKernel, "src", buffer[current]);
        shader.SetBuffer(initKernel, "dest", buffer[next]);
        shader.Dispatch(initKernel, NTGX, NTGY, 1);

        // Set Buffer for dumpKernel
        shader.SetBuffer(dumpKernel1, "src", buffer[current]);
        shader.SetBuffer(dumpKernel1, "dest", buffer[next]);
        shader.SetBuffer(dumpKernel2, "src", buffer[current]);
        shader.SetBuffer(dumpKernel2, "dest", buffer[next]);

        // It assumes the texture are properly named in the compute shader
        shader.SetTexture(dumpKernel1, "Result_T", Textures[(int) ImageTexture.TEX.Temperature]);
        shader.SetTexture(dumpKernel1, "Result_rho", Textures[(int)ImageTexture.TEX.Density]);
        shader.SetTexture(dumpKernel1, "Result_v", Textures[(int)ImageTexture.TEX.Velocity]);
        shader.SetTexture(dumpKernel1, "Result_GS1", Textures[(int)ImageTexture.TEX.GS1]);
        shader.SetTexture(dumpKernel1, "Result_GS2", Textures[(int)ImageTexture.TEX.GS2]);
        shader.SetTexture(dumpKernel1, "Result_Tox", Textures[(int)ImageTexture.TEX.Toxin]);
        shader.SetTexture(dumpKernel1, "Result_Ant", Textures[(int)ImageTexture.TEX.Antidote]);

        shader.SetTexture(dumpKernel2, "Result_Sti", Textures[(int)ImageTexture.TEX.Stimulus]);
        shader.SetTexture(dumpKernel2, "Result_Pre", Textures[(int)ImageTexture.TEX.Control]);
        shader.SetTexture(dumpKernel2, "Result_STM_S", Textures[(int)ImageTexture.TEX.ShortTermMemory_S]);
        shader.SetTexture(dumpKernel2, "Result_STM_P", Textures[(int)ImageTexture.TEX.ShortTermMemory_C]);
        shader.SetTexture(dumpKernel2, "Result_LTM", Textures[(int)ImageTexture.TEX.LongTermMemory]);

        shader.Dispatch(dumpKernel1, NTGX, NTGY, 1);
        shader.Dispatch(dumpKernel2, NTGX, NTGY, 1);

        Debug.Log(Parameters());

        GetInteractiveTextureCorners();
    }

    private void Update()
    {
        if (Input.GetKeyDown(KeyCode.Return)) running = !running;
        if (Input.GetKeyDown(KeyCode.Space) || running)
        {
            RunSimulation();
            if (compute_Bmean)
            {
                UpdateMaxConcentration();
                UpdateBmean();
                Nspots = CountSpots();
                if (time % RECORD_INTERVAL < 1.5 * STEPS * DT)
                {
                    times_list += "," + time.ToString();
                    Bmean_list += "," + Bmean.ToString();
                    Nspots_list += "," + Nspots.ToString();
                    Tox_list += "," + ComputeSpatialAverageTox().ToString();
                    Ant_list += "," + ComputeSpatialAverageAnt().ToString();

                    CheckStop();
                }
            }
        }        
    }

    void CheckStop()
    {
        if (auto_play & ((time > 3 * Bolus_shift + Bolus_quiet_time) | (time > Bolus_quiet_time && Bmean < 1e-9)))
        {
            CreateOutput();
            if (antidoteProduction == AntidoteProduction.Associative)
            {
                Application.Quit();
                running = false;
            }
            else
            {
                antidoteProduction += 1;
            }
            Start();
        }
    }


    void RunSimulation()
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
        
            shader.SetBool("inputMouse", true);
            
            shader.SetInt("mouse_i", mouse_i);
            shader.SetInt("mouse_j", mouse_j);

            if (Input.GetMouseButton(0))
            {
                shader.SetBool("remove", false);
            }

            if (Input.GetMouseButton(1))
            {
                shader.SetBool("remove", true);
            }

        }
        else
        {
            shader.SetBool("inputMouse", false);
        }
        for (int s = 0; s < S.STEPS; s++)
        {
            
            
            if (time < Bolus_quiet_time)
            {
                shader.SetFloat("dt", DT_FAST_FACTOR * S.DT);
                time += DT_FAST_FACTOR * DT;
            }
            else
            {
                shader.SetFloat("dt", S.DT);
                time += DT;
            }
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



    void UpdateBmean()
    {
        Bmean = 0;
        for (int i = 0; i < RESX; i++)
        {
            for (int j = 0; j < RESY; j++)
            {
                Bmean += currentState[i * RESY + j].y1;
            }
        }
        Bmean /= (RESX * RESY);        
    }

    float ComputeSpatialAverageTox()
    {
        float res = 0;
        for (int i = 0; i < RESX; i++)
        {
            for (int j = 0; j < RESY; j++)
            {                
                res += currentState[i * RESY + j].Tx;
            }
        }
        res /= (RESX * RESY);
        return res;
    }

    float ComputeSpatialAverageAnt()
    {
        float res = 0;
        for (int i = 0; i < RESX; i++)
        {
            for (int j = 0; j < RESY; j++)
            {
                res += currentState[i * RESY + j].An;
            }
        }
        res /= (RESX * RESY);
        return res;
    }

    void UpdateMaxConcentration()
    {
        MaxConcentration = 0;
        for (int i = 0; i < RESX; i++)
        {
            for (int j = 0; j < RESY; j++)
            {
                MaxConcentration = Mathf.Max(MaxConcentration, currentState[i * RESY + j].An);
            }
        }        
    }


    int CountSpots()
    {
        // First determine which are the pixels above the threshold
        const float THRESHOLD = -0.3f;
        int[,] to_segment = new int[RESX, RESY];
        for(int i = 0; i < RESX; i++)
        {
            for(int j = 0; j < RESY; j++)
            {
                if(currentState[i * RESY + j].y1 - currentState[i * RESY + j].x1 > THRESHOLD)
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


    bool IsInBounds(int i, int j)
    {
        return (i >= 0) && (i < RESX) && (j >= 0) && (j < RESY);
    }


    public static RenderTexture GetTexture(ImageTexture.TEX TEX_ID)
    {
        return Textures[(int) TEX_ID];
    }

    public static string Parameters()
    {
        string s;
        string sf = "FILE " + output_filename + "\n";
        s = "HeatTransfer={0}, NavierStokes={1}, ChemicalTransport={2}, ChemicalReactions={3}, RES = {4} x {5}, DX={6}, DT={7}, STEPS={8}, PER_X={9}, PER_Y={10}\n";
        sf += string.Format(CultureInfo.InvariantCulture, s, S.simu_Temperature, S.simu_v, S.simu_GS_transport, S.simu_GS_reactions, S.RESX, S.RESY, S.DX, S.DT, S.STEPS, S.PERX, S.PERY);

        s = "T0={0}, DeltaT={1}, LAMBDA={2}, Cv={3}, H={4}, OmegaT={5}, AmpliT={6}\n";
        sf += string.Format(CultureInfo.InvariantCulture, s, S.T0, S.DeltaT, S.LAMBDA, S.Cv, S.BoundaryTransferCoefficient, S.OmegaT, S.AmpliT);

        s = "RHO0={0}, Rs={1}, G={2}, mu={3}, Ra ={4:E2}, Gr ={5:E2}, Pr ={6:E2}\n";
        sf += string.Format(CultureInfo.InvariantCulture, s, S.RHO0, S.Rs, S.G, S.mu, Rayleigh(), Grashof(), Prandtl());

        s = "CHI_A={0}, CHI_B={1}, GS_F={2}, GS_R={3}, Ea={4}, k0={5}, DeltaH1={6}, DeltaH2={7}, p_init_B1={8}, p_init_B2={9},\n";
        sf += string.Format(CultureInfo.InvariantCulture, s, S.CHI_A, S.CHI_B, S.GS_F, S.GS_R, S.Ea, Mathf.Exp(S.Ea / S.T0), S.DeltaH1, S.DeltaH2, S.initialProportionOfPixelsB1, S.initialProportionOfPixelsB2);

        s = "Boluses quiet_time={0}, period={1}, time_shift={2}, concentration={3}, Sti decay={4}, Advection.x = {5}, Sti_diffusion={6}, P_diffusion={7}, Tox_diffusion={8}, Tox_decay={9}\n";
        sf += string.Format(CultureInfo.InvariantCulture, s, S.Bolus_quiet_time, S.Bolus_period, S.Bolus_delta, S.Bolus_concentration, S.Sti_decay, S.Bolus_advection_x, S.Sti_Diffusion, S.Pre_Diffusion, S.Tox_Diffusion, S.Tox_decay);

        s = "Memory={6} : STM_diff={0}, LTM_diff={1}, Ant_diff={2}, STM_decay={3}, LTM_decay={4}, Ant_decay={5}\n";
        sf += string.Format(CultureInfo.InvariantCulture, s, S.STM_Diffusion, S.LTM_Diffusion, S.Ant_Diffusion, S.STM_decay, S.LTM_decay, S.Ant_decay,S.antidoteProduction);

        s = "Reaction constants : kR={6}, kT={7}, kS={0}, kL={1}, kD={8}, kO={2}, kA={3}, k_ToxB={4}, k_AntB={5}";
        sf += string.Format(CultureInfo.InvariantCulture, s, S.kS, S.kL, S.kO, S.kA, S.k_ToxB, S.k_AntB, S.kR, S.kT, S.kD);
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

    static void SetGlobalParameters()
    {
        shader.SetFloat("time", 0);

        // General parameters
        shader.SetBool("simu_T", S.simu_Temperature);
        shader.SetBool("simu_v", S.simu_v);
        shader.SetBool("simu_GS_transport", S.simu_GS_transport);
        shader.SetBool("simu_GS_reactions", S.simu_GS_reactions);
        shader.SetInt("memory_mechanism", (int) S.antidoteProduction);
        //shader.SetBool("toxin_production", S.ToxinProduction);        

        // Interactivity
        
        shader.SetInt("inputSizeAdd", S.InputSizeAdd);
        shader.SetInt("inputSizeRemove", S.InputSizeRemove);

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
        shader.SetFloat("H", S.BoundaryTransferCoefficient);
        shader.SetFloat("OmegaT", S.OmegaT);
        shader.SetFloat("AmpliT", S.AmpliT);

        shader.SetFloat("RHO0", S.RHO0);
        shader.SetFloat("Rs", S.Rs);
        shader.SetFloat("G", S.G);
        shader.SetFloat("mu", S.mu);


        // Gray Scott
        shader.SetFloat("CHI_A", S.CHI_A);
        shader.SetFloat("CHI_B", S.CHI_B);
        shader.SetFloat("GS_F", S.GS_F);
        shader.SetFloat("GS_R", S.GS_R);
        shader.SetFloat("Ea", S.Ea);
        float k0 = Mathf.Exp(S.Ea / S.T0);
        shader.SetFloat("k0", k0);

        shader.SetFloat("DeltaH1", S.DeltaH1);
        shader.SetFloat("DeltaH2", S.DeltaH2);


        // Learning
        
        shader.SetBool("Bolus_AtBoundaries", S.Bolus_AtBoundaries);
        shader.SetBool("Bolus_AllBoundaries", S.Bolus_AllBoundaries);
        shader.SetFloat("Bolus_quiet_time", S.Bolus_quiet_time);
        shader.SetFloat("Bolus_period", S.Bolus_period);
        shader.SetFloat("Bolus_concentration", S.Bolus_concentration);
        shader.SetFloat("Bolus_duration", S.Bolus_duration);
        shader.SetFloat("Bolus_delta", S.Bolus_delta);
        shader.SetFloat("Bolus_shift", S.Bolus_shift);

        shader.SetFloat("CHI_Sti", S.Sti_Diffusion);
        //shader.SetFloat("CHI_Pre", S.Pre_Diffusion);
        shader.SetFloat("Sti_decay", S.Sti_decay);
        // no decay for Pre since it decays into Tox

        shader.SetFloat("Bolus_advection_x", S.Bolus_advection_x);

        shader.SetFloat("CHI_Tox", S.Tox_Diffusion);
        shader.SetFloat("Tox_decay", S.Tox_decay);
        shader.SetFloat("CHI_Ant", S.Ant_Diffusion);
        shader.SetFloat("Ant_decay", S.Ant_decay);

        shader.SetFloat("CHI_STM", S.STM_Diffusion);
        shader.SetFloat("STM_decay", S.STM_decay);
        shader.SetFloat("CHI_LTM", S.LTM_Diffusion);
        shader.SetFloat("LTM_decay", S.LTM_decay);

        shader.SetFloat("kR", S.kR);
        shader.SetFloat("kT", S.kT);
        shader.SetFloat("k_ToxB", S.k_ToxB);
        shader.SetFloat("k_AntB", S.k_AntB);

        shader.SetFloat("kS", S.kS);
        shader.SetFloat("kL", S.kL);        
        shader.SetFloat("kO", S.kO);
        shader.SetFloat("kD", S.kD);
        shader.SetFloat("kA", S.kA);

    }

    void OnValidate()
    {
        Debug.Log("Updated parameters");
        if(Time.time > 0)
        {
            SetGlobalParameters();
            GetInteractiveTextureCorners();
        }

    }

    void GetInteractiveTextureCorners()
    {
        // Update the interactive texture map        
        shader.SetInt("inputTexture", (int)InteractiveTexture);
        RawImage rawImage = ImageTexture.RawImageFromTex(InteractiveTexture);
        rawImage.rectTransform.GetWorldCorners(InteractiveTextureWorldCorners);
    }


    void OnApplicationQuit()
    {
        if (!auto_play) CreateOutput();        
    }

    void CreateOutput()
    {
        using (StreamWriter sw = File.CreateText("Outputs/" + output_filename))
        {
            sw.WriteLine(times_list);
            sw.WriteLine(Bmean_list);
            sw.WriteLine(Nspots_list);
            sw.WriteLine(Tox_list);
            sw.WriteLine(Ant_list);
        }
        Debug.Log("Output " + output_filename);
    }


}
