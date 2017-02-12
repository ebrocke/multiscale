# multiscale
[1] Brocke, E., Bhalla, U. S., Djurfeldt, M., Hellgren Kotaleski, J., &amp; Hanke, M. (2016). Efficient Integration of Coupled Electrical-Chemical Systems in Multiscale Neuronal Simulations. Frontiers in Computational Neuroscience, 10, 97. http://doi.org/10.3389/fncom.2016.00097
[2]

################################################################################
				Hierarchy of folders
################################################################################
.
├── createfigures	:	contains scripts for reading and plotting the results. See 'data' folder.
│   ├── createfigure.m
│   └── read_multirate_datapoints.m
├── data	:	contains *.mat object data files to generate figures in [1,2]. Each subfolder contains plotme.m script to generate the figure.
│   ├── multirate : multirate datafiles [2]
│   │   ├── Figure1
│   │   ├── Figure2
│   │   └── Figure3
│   └── refsol	:	contains *.mat data files with the reference solution
│       ├── ode15s_v2
│       └── ode15s_v3
├── model_v2	:	electrical-chemical  model with the fast communication signal, calcium is solved on the biochemical side
│   ├── common
│   ├── hh	:	electrical multicompartmental model with HH currents
│   ├── kkit	:	biochemical MAPK model
│   ├── modelconst.mat
│   ├── save_consts.m
│   ├── sim_combine2j2.m	:	can be used to obtains the solution with ode15s
│   └── yTypicalSolution.txt
├── model_v3	:	electrical-chemical  model with the slow communication signal, calcium is solved on the electrical side.
│   ├── common
│   ├── hh	:	electrical multicompartmental model with HH currents
│   ├── kkit	:	biochemical MAPK model
│   ├── modelconst.mat
│   ├── save_consts.m
│   ├── sim_combine2j2.m	:	can be used to obtain the solution with ode15s
│   └── yTypicalSolution.txt
├── README.md
└── solvers	:	contains implementation of the adaptive and fixed stepsize solvers
    ├── methods	:	numerical integration methods
    ├── sim_combine2j2.m	  :	      implements MAIN function of the coupled integration

################################################################################
				To run
################################################################################

sim_combine2j2.m has to be called to start the simulation. 
Parameters of interest can be assigned inside sim_combine2j2.m:
* MODEL	      : name of the model to simulate. Default: MODEL='model_v3'.
* iterMethod  : can be one of {'Jac','GSCELLFirst','GSERKFirst','GSSlowFirst','GSFastFirst'}.  Default: iterMethod = 'GSSlowFirst'.
* multirate   : <true, false>. Default: multirate = true.
* slv_param   : tolerance value (or the number of steps for a fixed step size controller). Default: slv_param = 1e-5.
* T_SIM	      : simulation time. Default: T_SIM = 2 [s].
* solvErkNM, solvCellNM : reference to the numerical method. Default: solvErkNM = solvCellNM = @BDF2_DEF.
* odeslv     : reference to the solver. Default: odeslv = @adaptive_solver.