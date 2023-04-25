## Loop Select
0. Guidance Mode
1. Angle Mode
2. Angular Rate Mode

### ADEPTSim.m
- Line 295, if `ctrl.loopselect == 1` then the aoa and sideslip angle commands change with time

### RunController.m
- Line 84, if `loopselect == 1` then overwrite `cmd_a` with `ctrl.cmd_a`, otherwise, `cmd_a` is calculated by the respective control law (NDI or INDI)
- Line 106, if `loopselect == 2` then overwrite `cmd_r` with `ctrl.cmd_r`, otherwise `cmd_r` is calculated by by multiplying the proportional gain for the angle controller, `kp_a` with the angle error `err_a` and then modifying that result based on equations in Josh's paper.