
# Notes

- Paper propses a modified #predictor-corrector algo that is able to meet constraints such as entry heating and maximum decelerative load
- Relies on the use of bank angle control to control the vertical component of aerodynamic lift
- Augments #predictor-corrector with objective-oriented altitude rate
- Separates reentry dynamics based on time scale
- #equations-of-motion for a gliding vehicle over a spherical, rotating Earth
	 1. $\dot r = V \sin \gamma$
	 2. $\dot \theta = \frac{V \cos \gamma \sin \psi}{r \cos \phi}$
	 3. $\dot \phi = \frac{V \cos \gamma \cos \psi}{r}$
	 4. $\dot V = -D-(\frac{\sin \gamma}{r^{2}}) \Omega^{2} r \cos\phi(\sin\gamma  \cos\phi  -  \cos\gamma  \sin\phi   \cos\psi)$
	 5.$\dot \gamma = \frac{1}{V}[L \cos\sigma + (V^{2}- \frac{1}{r})(\frac{\cos\gamma}{r})+ 2 \Omega V \cos\phi \sin\psi + \Omega^{2} r \cos\phi(\cos \gamma \cos \phi + \sin\gamma \cos \psi \sin \phi)]$
	 6.  $\dot  \psi = \frac{1}{V} [\frac{L\sin \sigma}{\cos \gamma} + \frac{V^{2}{r}\cos}\gamma \sin \psi \tan \phi - 2\Omega V (\tan \gamma \cos \psi \cos \phi - \sin \phi) + \frac{\Omega^{2}r}{\cos \gamma} \sin \psi \sin \phi \cos \phi]$
	 7. Where $r$ is the radial distance from the center of Earth to the CG of the vehicle, $\theta$ and $\phi$ are the longitude and latitude, V is the *Earth relative* velocity, $\gamma$ is the flight path angle of the Earth-relative velocity vector, and $\psi$ is the heading angle of the same velocity vector, measured clockwise in the local horizontal plane from the north.
	8. Length is normalized by Earth's equatorial radius $R_{eq} = 6,378,135 m$, time is normalized by $t_{scale} = \sqrt{\frac{R_{eq}}{g_0}}$ where $g_{0}= 9.81 m/s^2$ , velocity is normalized by $V_{scale} = \sqrt{g_oR_{eq}}$ 
	9. The differentiations in eq's 1-6 are with respect to dimensionless time $\tau = \frac{t}{t_{scale}}$ 
	10. $L$ and $D$ are the nondimensionalized lift and drag acceleration in $g_0$
	11. $\sigma$ is the bank angle, which is the roll angle of the vehicle about the relative velocity vector, positive to the right. **Note** The bank angle is not the same as the body roll angle when the angle of attack is not zero
	12. $\Omega$ is the rotation rate of the Earth
	13. Since time is not important, an energy-like variable, $e$ will be used as the independent variable in the algo where $e = \frac{1}{r} - \frac{V^2}{2}$ which makes it the negative of the specific mechanical energy used in orbital mechanics #TODO Look up specific mechanical energy in orbital mechanics
	14. If you ignore Earth rate, $\frac{de}{d\tau} = DV >0$ which makes $e$ a monotonically increasing variable.
	15. State vector $\mathbf{x} = \begin{pmatrix} r & \theta & \phi & \gamma & \psi \end{pmatrix}^T$  
	16. Velocity is detemined by $V = \sqrt{2(\frac{1}{r - e})}$
	17. Angle of attack $\alpha$ enters in $L$ and $D$ through the dependence on $\alpha$ by $C_L$ and $C_D$ 
	18. As in most entry guidance developments, the angle-of-attack profile is considered to be fixed as a given function of the Mach number or relative velocity. This profile is determined by considerations in ranging capability (in the downrange and/or crossrange direction), thermal protection, and flight control. However, the guidance algorithm to be developed is independent of any particular angle-of-attack profile.
	19. Eq's 1-6 (excluding 4) can be rewritten with $e$ as the independent variable like $\mathbf{x}' = \frac{d\mathbf{x}}{de} = \mathbf{f}(\mathbf{x}, \sigma)$ where $\mathbf{x}(e_{0}) = \mathbf{x}_0$ 
	20. End constraints
		1. Specified distance from target location $s(\tau_{f}) = s_{f}^{\star}$  where s denotes a great circle range to the landing site, normalized by distance $R_{eq}$ which puts $s$ in radiance and a function of long and lat. 
		2. Specified altitude $r(\tau_{f}) = r_{f}^{\star}$
		3. Specified velocity $V(\tau_{f}) = V_{f}^{\star}$
		4. The altitude and velocity constraints can be combined together in terms of energy, giving $e_{f} = \frac{1}{r_{f}^{\star}} - \frac{V_{f}^{\star 2}}{2}$ giving a final single constrain in terms of energy of $s(e_{f}) = s_{f}^{\star}$ 
	21. Path constraints 
		1. Heating rate $\dot Q = k_{Q} \sqrt{\rho} V^{3.15} \leq \dot Q_{max}$ units: ($\frac{W}{m^2}$) at a stagnation poiont on the surface of a vehicle with a curvature radius of 0.3048m and $k_{Q} = 9.4369\times10^{-5}\times(\sqrt{g_{0}R_{eq}})^{3.15}$ V is non-dimensional as defined before, $\rho$ has the units of $\frac{kg}{m^3}$
		2. Load factor $a = \sqrt{L^{2}+D^{2}} \leq a_{max}$ units of $a_{max}$ are in g's
		3. Dynamic pressure $\bar q = \frac{g_{0} R_{eq} \rho V^{2}}{2} \leq \bar q_{max}$ 
	22. The goal is to find the bank-angle command at each instant, based on the current state such that the trajectory of the system will satisfy the boundary conditions on $\mathbf{x}(e_{0})$ and $s(e_{f}) = s_{f}^{\star}$ as well as the path constraints on heating, load factor, and dynamic pressure outlined in bullet 21
-  Baseline predictor-corrector guidance algorithm determines a complete profile of the bank-angle magnitude and the corresponding feasible logitudinal trajectory from the current point to the end of the entry phase.
-  Let $s$ again denote the range to go (in radians) on the surface of a spherical Earth along the great circle connecting the current location of the vehicle and the site of the final destination.  When the offset between the azimuth of this great circle and the heading angle is ignored, the differential equation that governs $s$ is
	1. $\dot s = \frac{\mathrm{d}s}{\mathrm{d}\tau} = -V \cos \frac{\gamma}{r}$ 
	2. Let $\mathbf{y} = \begin{pmatrix} r & \theta & \phi & \gamma & \psi \end{pmatrix}^T$ denote the #equations-of-motion equations 1-3, 5-6, and the $\dot s$ equation above by $\frac{\mathrm{d}y}{\mathrm{d}e} = \hat f(\mathbf{y}, \sigma, e)$  and let $V = \sqrt{2(\frac{1}{r - e})}$ be used wherever $V$ is needed
	3. It is essential to include $s$ in $\frac{\mathrm{d}y}{\mathrm{d}e} = \hat f(\mathbf{y}, \sigma, e)$  instead of representing the range to go by the great-circle arc length between a point on the trajectory and the landing site, computed by spherical trigonometric functions. The result from the spherical trigonometric functions cannot distinguish an undershoot case from an overshoot case, which will make the convergence of the the subsequent algorithm very difficult.
- No assumption on the smallness of the magnitude of $\gamma$ is needed, terefore this is a good approach for low-lift RVs like capsules which will see a large negative $\gamma$ towards the end of the trajectory
- 
	
	