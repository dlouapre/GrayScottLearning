# Gray-Scott Learning

This is the source code for our paper

Provenance of Life: Chemical Autonomous Agents Surviving through Associative Learning
by Stuart Bartlett and David Louapre

Physical Review E https://journals.aps.org/pre/abstract/10.1103/PhysRevE.106.034401

DOI https://doi.org/10.1103/PhysRevE.106.034401


Video abstract of the paper : https://www.youtube.com/watch?v=xZlXIVKBbcc 

Video demo of the code : https://youtu.be/vQzWdrHA4Cw

david.louapre@gmail.com

## General

This code runs on Unity3D, and has been created using version 2020.3.22f1, although it should be able to run on slightly different versions of Unity (some project conversion might happen when loaded for the first time). Some basic knowledge of the Unity editor is required (hierarchy, editor, inspector, running the game...)

The code heavily uses GPU for simulation, through the use of compute shaders, so it might not run well on old machines. You can always check how fast the code is running by looking at the FPS in the "Stats" tab of the Unity "Game" panel.


## What is simulated
Simulated variables and systems include :

- Chemical compounds A and B of the Gray-Scott reaction diffusion system
- Optionnaly a second GS system that exist in parallel of the first one
- Optionnaly temperature and heat transfer through diffusion and advection
- Optionnaly density and velocity of a fluid flow using a compressible perfect fluid approximation 
- Chemical compounds for the associative learning network described in the paper : Stimulus, Toxin, Antidote, Short term memory, Long term memory

## Running existing scenarios
The simplest way to use the code is to look at existing scenarios. For this, in the "Project" panel of Unity go to Assets>Scenes and double click on a scene name to load it. It should then appear in the "Hierarchy" panel of Unity and the corresponding display should load in the "Game" panel. A quick description of the scenario is written at the top of the Game view.

(see https://learn.unity.com/tutorial/exploring-the-editor-layout for Unity editor layout description)

In the Hierarchy panel, clicking on "SimulationManager" will load it in the "Inspector" panel. There you can visualise and change all the simulation parameters.

To start the simulation, click on the Play button, the scene will start in the Game panel. You will be able to visualize some of the variables of the simulation. 

The color code for each variable is defined in the compute shader file SimulationEngine in Assets>Resources>Shaders) in the CSDump kernels.

- Gray Scott order parameter B - A : [-1,0]
- Temperature [1,2]
- Density between [0.5,1.5]
- Flow velocity is color coded : hue for direction and value for the norm of the velocity [0,1]
- Chemical compounds Tox, Ant, Sti, M, L : concentration [0,1]

The lower text box shows most of the parameter values.

When the scene is running, you can start/pause the simulation by pressing `Return`. When the simulation is paused, you can run it step by step by pressing `Space` key.

When the simulation is running, you can interact with any of the displayed variables by using the mouse. Left click will add, right click will remove (temperature, concentration, density...)

You can adjust the size and the strength of the mouse interaction using arrows on the keyboard, values can be read in the lower text box (mouse input size & strength)

When you quit the simulation, an output file is created in the "output" folder. Its name reflects date, time and scenario. First line contains time stamps, second line the average B concentration, third line the number of Gray-Scott spots.

Parameters can be tweaked in the inspector panel, either before starting the simulation, or while it is running (but remember that any change made while the simulation is running will be reversed when exiting the simulation)

You can always duplicate an existing scenario for your own experiments : in the Project panel, click on a scene name in Assets>Scene and press Cmd+D / Ctrl-D to duplicate

### Parameters
This is a list of the parameters exposed in the SimulationManager in Inspector

- Auto_play : whether the simulation will start right away, or if `Return` should be pressed to launch it (or `Space` to run it step-by-step)
- RESX, RESY : resolution in pixels of the simulation. WARNING : Should be a multiple of 32. Typical value 256 x 128
- DX : cell size (1)
- PERX, PERY : whether periodic boundary conditions are applied in each dimension
- DT : time step (0.1 for Gray-Scott, 0.01 if compressible fluid simulation is running)
- STEPS : how many time steps of the simulation per frame (10 to 100), you can adjust based on the FPS displayed in the Stats tab of the Game panel when the simulation is running
- CHI_A,B : diffusion constant of the GS components, 0.2 and 0.1
- GS_F : the GS reaction parameter F : 0.03
- GS_R : the GS reaction parameter F : 0.061
- K0 : base reaction constant for the autocatalytic reaction in GS : 1
- Initial proportion of pixels : for random initialization of B component : 0 (to start empty) or 0.3
- Ea : activation energy for the GS reaction. Either 0 for temperature insensitivity, or typically 1.7 (in which case use k0=3.1, see homeostasis scenario & paper)
- Delta H : enthalpy for GS reactions. Negative is exothermic (-0.02) Positive is endothermic (0.025)
- Antidote production : mode of antidote production/learning, see paper.
- Bolus at boundaries : select for diffusion delivery of Tox and Sti, otherwise uniform delivery will be used
- Bolus quiet time : initial time without boluses, so that the GS system can develop (2000)
- Bolus period : time between two boluses (100)
- Bolus duration (1)
- Bolus concentration (1)
- Bolus delta : time between Stimulus and Toxin (50)
- Bolus shift : characteristic time for the evolution of the environment parameter. After quiet time, it represents the time to grow to maximal toxin delivery.
- _diffusion : diffusion constant of the different chemical compounds, typically 0.1 to 1, except for boundary diffusion delivery of the stimulus and toxin, in which case higher values have to be used (in that case you should decrease DT since the maximally allowed diffusion value for stability is (0.25 - DX - DX / DT)
- _decay : decay rate
Reaction constants
- k_AntTox : N + T => 0
- k_ToxB : T + B => 0
- k_AntB : N + B => 0   (antidote mild impact on B, should be smaller than k_ToxB and k_AntTox)
- k_Reactive : antidote production indirect reaction network
- k_Preemptive : antidote production inpremeptive reaction
- k_STM : short term memory production
- k_LTM : long term memory production
- k_Output : antidote production in associative network
- simu_temperature : whether to simulate heat transfer
- T0 : base temperature
- Delta_T : static difference between top and bottom boundary temperature (for Rayleigh Benard)
- Lambda : thermal conductivity
- Cv : heat capacity
- Omega_T : pulsation for oscillating boundary conditions
- Ampli_T : amplitude for oscillating boundary conditions
- Boundary Transfer Coefficient : transfer coefficient at boundary (h)
- Initial temperature gradient : whether to start with initial temperature gradient (governed by DeltaT)
- Simu_v : whether to simulate fluid flow velocity
- VX0, VY0 : stationnary flow velocity
- RHO0 : base density (1)
- Rs : base specific constant (so that P = Rs - rho - T) (10)
- G : Gravity (0.05)
- Mu : viscosity (0.1), should be lower than 0.25*DX*DX*DT (see diffusion constant)

### Existing scenes

#### Minimal GS
A simple Gray-Scott reaction diffusion system. There is no B at the beginning so it has to be added with a left click of the mouse.

#### Pure GS
Same but with a large resolution and random initial conditions for B. Press <Return> to launch the simulation.

#### GS thermostatic
A Gray-Scott system with a temperature dependent reaction rate, inside a temperature gradient. We can see lamellar structure at the bottom, solitons in the middle and nothing at the top.

#### GS stationnary flow
A GS system in a stationnary fluid flow with periodic boundary conditions. Velocity is high so spots have a hard time surviving. They can be created by holding the left mouse button for 1 second and release it. The soliton will flow.

#### GS homeostasis
Demonstration of thermal homeostasy similar to Bartlett & Bullock "A precarious existence: Thermal homeostasis of simple dissipative structures."
Two independent GS systems are used. The first one is exothermic, the second one endothermic.
Boundary temperature slowly decreases.
(Heat transfer only, no fluid flow)

#### Rayleigh-Benard
GS within a compressible fluid flow with Rayleigh-Benard convection. Temperature boundary conditions induce forced convection forming convection cells. GS spots can survive inside the convection cells.

#### Uniform learning
Associative learning network described in the paper, with uniform delivery of stimulus & toxin. Following the paper, learning mode can be changed in the parameters of the SimulationManager. (none, direct, preemptive, associative)

#### Diffusion learning
Associative learning network described in the paper, with boundary diffusion delivery of stimulus & toxin. Following the paper, learning mode can be changed in the parameters of the SimulationManager. (none, direct, preemptive, associative)

## Creating new scenario

To create a new scenario, you can always duplicate an existing one and build from here. To change the variable that is displayed in each window, in the Hierarchy panel go to Canvas : every DisplayPanel corresponds to a variable. When opening it in the Inspector, the dropdown menu besides TEX_ID indicates the variable that will be displayed.

You can create and modify more DisplayPanel by duplicating an existing one or draging a "DisplayPanel" from the Project>Assets>Prefabs into the Canvas. Position and scale of the panel can be modified in the Rect Transform section in the Inspector.

## Modifying the simulation
### Code organisation
If you want to modify the details of the simulation, tweak the code, add some variables or equations, etc. you have to be aware that this is not straightfoward. The reason is for performance the code is run on the GPU, which makes the code fairly complicated to handle and debug. 

If you feel comfortable doing it, here are some guidelines.

The first essential file is the SimulationManager.cs script. It is located in Assets>Scripts>Simulation. That files handles the values of the different parameters, pass them to the compute shader, calls it and gather the results to be displayed.

The compute shader code itself is in Assets>Resources>Shaders. It is written in HLSL. The main file is SimulationEngine, and it calls the other shader files like ChemicalReactions, ThermalFluid and ShaderAux.

For performance reasons, you will notice some #define preprocessor directives in both the SimulationManager (SM) and SimulationEngine (SE). They can be commented out to remove some sections of the simulation when they are not used. Be sure to do exactly the same in both SM and SE files.

### Adding variables and equations
If you want to add a new type of variables, some new equations, etc., here is the procedure.

It should be done in a consistent way in both SimulationManager.cs (SM) and SimulationEngine.compute (SE)
If needed, define a corresponding preprocessor directive in SM and SE

#### New variable definition

- Add it to CellData in SM
- Add it to Point in SE (same order as CellData in SM)

#### Texture definition (if you want to visualize that variable as a texture)

- Define the corresponding RWTexture2D in SE
- Define the associated texture in enum TEX in ImageTexture
- Set the texture link in SetTextures() in SM
- Populate the texture in the dump kernel in SE with the relevant variable

#### Evolution

Define the relevant evolution equations in the Update() part of the compute shader kernel in SE
For efficiency it is better do it such that dest[] is not read, and is written only once.

When adding an equation, if parameters are needed they need : 

- to be added in the inspector
- corresponding parameter should be added to cbuffer GLOB
- Parameters should be then passed to the cbuffer GLOB in SetGlobalParameters()

Optionnaly expose the parameter value in Parameters()

Optionnaly : If you want interactivity on this variable, add it to ApplyMouseInput in SE, the number should be the same as the enum TEX

