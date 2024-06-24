# Examples

## Example Input file
This example simulates a rectangular VO<sub>2</sub> device, supplied with a direct voltage through a series resistor. The input file is written in xml form. Below is an example of the input file.
```xml
<?xml version="1.0"?>
<input>
 <external>
  <temperature unit='K'>300.1</temperature>
  <voltage unit='V'>5.6</voltage>
  <resistor unit='Ohm'>2.5e5</resistor>
  <capacitor unit='nF'>0.0</capacitor>
  <heatdiss unit='W/m^2K'>3e6</heatdiss>
  <Lx unit='nm' mesh='100'>100.0</Lx>
  <Ly unit='nm' mesh='40'>40.0</Ly>
  <Lz unit='nm'>20.0</Lz>
 </external>
 <time>
  <endtime unit='ns'>2e3</endtime>
  <savemethod>auto</savemethod>
  <saveperiod>20</saveperiod>
 </time>
 <initialization>
  <temperature unit='K'>300.0</temperature>
  <SOP>1.119</SOP>
  <EOP>-1.293</EOP>
  <Tcvariance method='nucleus1'>
   <Tcshift unit='K'>-5.0</Tcshift>
   <radius unit='nm'>3.0</radius>
  </Tcvariance>
 </initialization>
 <solverparameters>
  <Newtonabsolutetolerance>1e-5</Newtonabsolutetolerance>
  <Newtonrelativetolerance>1e-3</Newtonrelativetolerance>
  <Newtonmaxiteration>15</Newtonmaxiteration>
  <timesteptolerance>1e-2</timesteptolerance>
  <directsolver>pastix</directsolver>
  <loglevel>INFO</loglevel>
 </solverparameters>
</input>
```

Almost all the parameters are self-explanatory. The units are fixed and serve only to remind the user of the corresponding parameter's unit. The `external` section defines external parameters: 

Name          | Explanation
------------- | -------------------
`temperature` | Ambient temperature
`voltage`     | Direct voltage applied
`resistor`    | Resistance of the series resistor
`capacitor`   | The capacitance representing the parasitic capacitance or the external capacitor parallelly connected with the series resistor
`heatdiss`    | Heat transfer coefficient from the device to the environment
`Lx`          | Width of the device
`Ly`          | Length of the device, along the electric field
`Lz`          | Thickness of the device

The `time` section defines the parameters related to the time:

Name         | Explanation
------------ | -----------
`endtime`    | Simulation end time; beginning time is zero
`savemethod` | Method of saving time-dependent solutions, used with `saveperiod` parameter. `auto` means to save every `saveperiod`-step solution; `fixed` means to save every `saveperiod`-nanosecond solution
`saveperiod` | See `savemethod`

The `initialization` section defines parameters for initialization:

Name          | Explanation
------------- | -----------
`temperature` | Initial temperature of the device
`SOP`         | Initial value for the structural order parameter
`EOP`         | Initial value for the electronic order parameter
`Tcvariance`  | `method`: how to set up a nucleus of the high-temperature phase. `nucleus1` means to set up a semicircle with a radius of `radius` and a transition temperature shift of `Tcshift`, located at the midpoint of the $y = 0$ edge

The `solverparameters` section defines the parameters for both the Newton's iteration solver (for nonlinear differential equations) and the linear solver (for solving linear equations in each iteration of the Newton's method):

Name                      | Explanation
------------------------- | -----------
`Newtonabsolutetolerance` | Absolute tolerance for the Newton iteration
`Newtonrelativetolerance` | Relative tolerance for the Newton iteration. Meeting either the absolute tolerance or the relative tolerance is considered converged
`Newtonmaxiteration`      | Limit of the iteration times for the Newton's method
`timesteptolerance`       | Relative tolerance for the adaptive time stepping error
`directsolver`            | Which direct solver to use for solving the linear problem
`loglevel`                | Log level; see [FEniCS manual](https://fenics.readthedocs.io/projects/dolfin/en/2017.2.0/apis/api_log.html "FEniCS log level")