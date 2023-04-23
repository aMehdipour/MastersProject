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