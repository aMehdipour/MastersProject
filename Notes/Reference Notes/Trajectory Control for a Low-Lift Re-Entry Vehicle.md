
# Notes

- Vehicle must be able to fly within its own safety limits
- Must be able to respond to off-nominal entry conditions
- In the past, trajectory control was mainly based on path prediction, in which the controller computes a new trajectory in advance based on the actual flight conditions, and path following, in which the controller tries to reach a predetermined trajectory.
- Traditionally, reference trajectories are selected from the acceleration - velocity plane, since range can be expressed as the integral of drag over velocity. Thus, measuring velocity, a vehicle is guided to the desired final state by scheduling drag or vertical acceleration
	- This means that a certain drag and acceleration command is pulled from a pre-determined schedule/trajectory
- Gemini and Apollo track a stored L/D reference
	- Schedule is corrected based on estimated downrage error, and acceleration and velocity feedback are used to eliminate tracking errors
- Self-reliant predictor methods have been suggested for #low-lift entry capsules
	- #low-lift guidance in refs 2 and 13 #TODO fill out refs
	- Based on current flight conditions, these systems use numerical fast-time integration of the future trajectory to determine the control input during re-entry.
- The numerical fast-time techniques can be combined with optimization algorithms in refs 14-15 #TODO fill out refs
- Method in this paper used pre-determined trajectories for alt, vel, and $\gamma$ 
- The paper investigates perturbation-guidance for #low-lift RVs
- Is reliant on GPS for pos and vel updates
- Uses an H$\infty$ quadratic cost function
	- Linear quadratic regulator optimization ref 25 #TODO fill out refs
- 


## From ChatGpt

- #Perturbation-Guidance 
	Perturbation guidance is a control technique used in guidance and control systems for aerospace vehicles. It involves continuously making small adjustments to the vehicle's trajectory in response to changes in the vehicle's environment, such as variations in atmospheric conditions or external disturbances.
	
	In perturbation guidance, the vehicle's trajectory is computed based on a nominal or pre-planned trajectory, which serves as the initial guidance path. The nominal trajectory is perturbed or adjusted continuously during flight, in response to the changing environment or to ensure that the vehicle follows the desired path accurately.
	
	The adjustments made to the trajectory are usually small and incremental, and they are calculated based on feedback from sensors on the vehicle that measure its current position, velocity, and other relevant parameters. The perturbations are typically designed to be small enough to ensure that the vehicle remains within its safe operating range and to prevent excessive fuel consumption.
	
	Perturbation guidance is commonly used in space missions, where the vehicles have to navigate through complex and unpredictable environments. It is also used in aircraft autopilot systems, where it can help to ensure that the aircraft remains on course and avoids collisions with other aircraft or obstacles.
