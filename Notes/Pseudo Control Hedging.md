Pseudo Control Hedging (PCH) is a technique used in #control-allocation to handle situations where control saturation occurs. Control saturation happens when the actuators' limits are reached, and they can no longer provide the necessary control effort to follow the desired command. This can lead to degraded performance or even instability in the system.

#controls #control

The primary idea behind PCH is to adjust the control commands in a way that ensures the most critical control objectives are met while minimizing the deviation from the desired command. This is achieved by prioritizing control objectives and redistributing the control effort within the available actuator limits.

Here's a brief overview of how Pseudo Control Hedging works:

1.  Define the control objectives: First, you need to identify and prioritize the control objectives (e.g., maintaining altitude, controlling attitude, etc.) that your system needs to achieve.
    
2.  Compute desired control commands: Based on the control objectives and the control law (NDI or INDI), compute the desired control commands for the actuators.
    
3.  Check for saturation: Determine whether any of the actuators are saturated, i.e., reaching their limits. If there is no saturation, the desired control commands can be directly applied.
    
4.  Apply PCH: If saturation occurs, apply the PCH technique to modify the control commands. This is done by redistributing the control effort among the actuators in a way that ensures the most critical control objectives are met while minimizing the deviation from the desired commands.
    
5.  Implement adjusted control commands: Apply the adjusted control commands to the actuators to maintain control of the system.
    

PCH allows the control system to continue functioning effectively even when some actuators reach their limits by prioritizing the most important control objectives and redistributing the control effort accordingly. This makes it a valuable technique for handling control saturation in various control systems, including re-entry vehicles.

Remember that the implementation of PCH might require some tuning and customization based on the specific characteristics of your control system and actuators. Be prepared to experiment with different configurations to find the best balance between performance and control effort distribution.

Control law from NDI calculates the desired actuator command by inverting the classical state space equation $\dot{x} = a(x) + b(x)\delta_{comm}$ where $\delta_{comm}$ is the commanded actuator deflection. Inversion of the the equation to solve for the actuator command yields $\delta_{comm} = b^{-1}(x)[\nu - a(x)]$ where $\nu$ is the virtual control input. The virtual control input can be estimated by $\hat{\nu} = a(x) + b(x)\delta_{act}$ where $\delta_{act}$ is the actual control input with actuator position and rate limits taken into account. The PCH signal $\nu_h$ can be calculated by subtracting the estimated from the commanded virtual input: $\nu_{h} = \nu - \hat{\nu}$ 

## Pseudo Control Hedging and its Application for Safe Flight Envelope Protection

- 
