## Configure FABM

There are various configurations for the use of FABM

### YAML File

Set as `FABMConfigFile` in `FABMConfig`, in `.yaml` format. 

`instances` lists all modules that are simulated in the code, with the following sub-configurations:

- `long_name` is the long name of the module used in the output
- `model` is the biogeochemical model used by the module
- `parameters` sets parameters for the module
- `initialization` sets the scalar initial values of the variables in the module
- `coupling` sets the coupling to variables from other modules

### Initial Conditions

To add depth-varying initial conditions for a FABM variable, add an initial conditions file of the same format as the Simstrat initial conditions to the folder specified as `FABMInitialPath` in `FABMConfig`.

### Inflow

To add inflow for a FABM variable, add an inflow file of the same format as the Simstrat inflows to the folder specified as `FABMInflowPath` in `FABMConfig`.

### Manipulations

The file `manipulations.nml` in `FABMConfigPath` in `FABMConfig` allows for manipulations of FABM variables during the simulation. Set the following options:

- `var_name` is the variable short name, must be among the FABM state variables (a list can be found in `_variables.dat` in `Path` in `Output`)
- `start_time` is the time in days at which the manipulation starts, must be smaller than `End d` in `Simulation`
- `end_time` is the time in days at which the manipulation stops, must be larger than `Start d` in `Simulation`
- `action_type` sets the type of action of the manipulation, there are two options:
    - `1` for a multiplication only executed once at `start_time` (`end_time` is irrelevant for this option)
    - `2` for an addition executed at every timestep between `start_time` and `end_time`
- `action_val` is the value of the action. For an action of type `2` it is multiplied by the timestep and divided by the depth range
- `threshold` is only relevant for an action of type `2`, where the action is not executed if all values in the depth range are above `threshold`
- `start_depth` is the depth below which the action is executed, cannot be lower than the lowest depth 
- `end_depth` is the depth above which the action is executed, cannot be greater than the highest depth 

> **N.B.** Depth input is converted to the nearest grid point for the simulation.


For example for the variable `npzd phytoplankton`, to double its density at day 12005 between 20 and 60 meters depth, add the following manipulation:

~~~bash
&manipulation
   var_name = 'npzd_phy',
   start_time = 12005,
   action_type = 1,
   action_val = 2.0,
   start_depth = -20,
   end_depth = -60
/
~~~

Whereas to add a density of 10<sup>-5</sup> mmol/m<sup>3</sup> during days 12005 to 12010 between 80 and 120 meters depth, whenever the density falls below 10<sup>-4</sup> mmol/m<sup>3</sup>, add the following manipulation:

~~~bash
&manipulation
   var_name = 'npzd_phy',
   start_time = 12005,
   end_time = 12010,
   action_type = 2,
   action_val = 0.00001,
   threshold = 0.0001,
   start_depth = -80,
   end_depth = -120
/
~~~

### Output Diagnostic Variables

To output diagnostic variables, first set `OutputDiagnosticVars` in `FABMConfig` to `true`. Then go to `FABMConfigPath` in `FABMConfig`.

After having started the simulation once, you will find a list of all diagnostic variables under `list_diagnostic_interior.dat` and `list_diagnostic_horizontal.dat` for interior and horizontal (bottom and surface) diagnostic variables, respectively. In the file `output_diagnostics.dat`, you can add the diagnostic variables to output with the following options:

- Add the short name of a FABM diagnostic variable to output (e.g. `attenuation_coefficient_of_photosynthetic_radiative_flux`)
- Add `select_all_*` to output all FABM diagnostic variables that start with `*` (e.g. `select_all_npzd`) or `select_all_` to output all FABM diagnostic variables
- Add `select_output_*` to output all FABM diagnostic variables that start with `*` (e.g. `select_output_npzd`) and have `Output` set `Yes` by the biogeochemical model or `select_output_` to output all FABM diagnostic variables that have Output set true

### Output Repaired Variables

FABM repairs variables below minimum or above maximum by restricting them to their extrema, if `RepairStates` is set to `true` in `FABMConfig`. To output repaired variables, first set `OutputRepairedVars` in `FABMConfig` to `true`. Then go to `FABMConfigPath` in `FABMConfig`.

After having run the simulation once, you will find a list of all repaired variables and the boundary they have crossed under `list_repaired.dat`. Upon re-running the simulation, next to the variable output there is an additional output for the variable showing:

- The value of the variable had it not been repaired, if the boundary is crossed by the variable
- The boundary if the variable stays within bounds

### Further Configurations

The following are additional configurations that can be found in `FABMConfig`:

- `BottomEverywhere` selects whether FABM Bottom State variables should be calculated at every layer (`True`) or only at the lowermost layer (`False`)
- `BioshadeFeedback` selects whether light extinction should be calculated by the biogeochemical models (`True`) or from input (`False`)
- `BackgroundExtinction` sets background extinction due to pure water extinction added to the light extinction calculated by the biogeochemical models if `BioshadeFeedback` is set to `True`

> **N.B.** If `BottomEverywhere` is set to `True`, also all horizontal diagnostic variables have output at every layer, even if for some of them (e.g. surface diagnostic variables) this is not physically meaningful. This does not affect the simulation.