module Movement

using LinearAlgebra: norm, cross
using DifferentialEquations
using FLOWUnsteady

#using ResumableFunctions
#using SimJulia

using Dierckx

export CMassProp, Reference
export AbstractAeroModel, AbstractPropulsionModel, AbstractInertialModel,
    AbstractAtmosphereModel
export UniformGravitationalField,
    ConstantAtmosphere
export CO, COUNTER, COCOUNTER
export equations


# Aliases
uns = FLOWUnsteady
vlm = FLOWUnsteady.vlm


# ------ General Structs -------

"""
    State(x, y, z, phi, theta, psi, u, v, w, p, q, r)

State of the aircraft: positions in inertial frame, euler angles,
velocities in body frame, angular velocities in body frame.
"""
struct State{TF}
    x::TF  # position (inertial frame)
    y::TF
    z::TF
    phi::TF  # orientation, euler angles
    theta::TF
    psi::TF
    u::TF  # velocity (body frame)
    v::TF
    w::TF
    p::TF  # angular velocity (body frame)
    q::TF
    r::TF
end

"""
    MassProp(m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz)

Mass and moments of inertia in the body frame.
Ixx = int(y^2 + z^2, dm)
Ixz = int(xz, dm)

Most aircraft are symmetric about y and so there is a convenience method
to specify only the nonzero components.
MassProp(m, Ixx, Iyy, Izz, Ixz)
"""
struct MassProp{TF}
    m::TF
    Ixx::TF
    Iyy::TF
    Izz::TF
    Ixz::TF
    Ixy::TF
    Iyz::TF
end

# most aircraft are symmetric in y
MassProp(m, Ixx, Iyy, Izz, Ixz) = MassProp(m, Ixx, Iyy, Izz, Ixz, zero(Ixx), zero(Ixx))

"""
function necessary to evaluate the mass properties of the structure: in input
    are given the mass in kg and number of rotors, the radius of center mass in m,
    the distance between the center mass and the rotors. For output is given the
    total mass, the moment of inertia respect to the x,y and z axis
"""
function massprop(m,Mmotor,l,rotors,Rc)
    mtot=m+enumerate(rotors)*Mmotor
    Ixx= (2/5)*m*Rc^2+2*l^2*Mmotor
    Iyy= Ixx
    Izz= (2/5)*m*Rc^2+4*l^2*Mmotor
    return mtot,Ixx,Iyy,Izz
end

"""
    Reference(S, b, c)

The reference area, span, and chord used in the aerodynamic computations.
"""
struct Reference{TF}
    S::TF  # area
end

"""
function used to elaborate the referenc surface, that is equal to the sum of the
    each rotor surface
"""
function reference(Rtip,rotors)
    Sref=Rtip^2*pi*enumerate(rotors)
    return Sref
end

"""
    inertialtobody(state)

Construct a rotation matrix from inertial frame to body frame

The assumed order of rotation is
1) psi radians about the z axis,
2) theta radians about the y axis,
3) phi radians about the x axis.

This is an orthogonal transformation so its inverse is its transpose.
"""
#function inertialtobody(state)

#    R = Array{eltype(state.phi)}(undef, 3, 3)

#    cphi, ctht, cpsi = cos.([state.phi, state.theta, state.psi])
#    sphi, stht, spsi = sin.([state.phi, state.theta, state.psi])

#    R[1, 1] = ctht*cpsi
#    R[1, 2] = ctht*spsi
#    R[1, 3] = -stht

#    R[2, 1] = sphi*stht*cpsi - cphi*spsi
 #   R[2, 2] = sphi*stht*spsi + cphi*cpsi
 #   R[2, 3] = sphi*ctht

#    R[3, 1] = cphi*stht*cpsi + sphi*spsi
#    R[3, 2] = cphi*stht*spsi - sphi*cpsi
#    R[3, 3] = cphi*ctht

#    return R

#end

"""
    windtobody(alpha, beta)

Rotation matrix from wind frame to body frame.
- alpha: angle of attack
- beta: sideslip angle
"""
#function windtobody(alpha, beta)

#    ca, cb = cos.([alpha, beta])
#    sa, sb = sin.([alpha, beta])

#    Rwb = [ca*cb  -ca*sb  -sa;
#           sb      cb     0.0;
 #          sa*cb  -sa*sb  ca]

#    return Rwb
#end





"""
    ConstantAtmosphere(Wi, Wb, rho, asound, g)
    TV defines vectors, TF defines a constant number

Constant atmospheric properties.
"""
struct ConstantAtmosphere{TF, TV}
    Wi::TV
    Wb::TV
    rho::TF
    asound::TF
    g::TF
end






"""

Computates the Thrust and Torque of each rotor in time
"""

function run (tinit, tfinal)

   for (t in seq(tinit, tfinal, by=tstep))
  	 for rotor in rotors
#	 	 prfrmnc, gammas = calc_distributedloads(self, Vinf, RPM, rho; t=t,
 #                                       include_comps=include_comps,
#                                        return_performance=return_performance,
#                                        Vref=Vref, sound_spd=sound_spd,
#                                        Uinds=Uinds,
#                                        _lookuptable=_lookuptable, _Vinds=_Vinds,
#                                        tiploss_correction=tiploss_correction,
#                                        AR_to_360extrap=AR_to_360extrap)                         
#     	  	 CT[rotor][t] = prfrmc[1][2]
#     	   	 CQ[rotor][t] = prfrmc[1][3]
#		T[rotor][t] = CT[rotor]*(rho*n^2*D^4)
#		Q[rotor][t] = CQ[rotor]*(rho*n^2*D^5) 
		T[rotor][t], Q[rotor][t] = calc_thrust_torque(Rotor)
     	   	end

   end
   return T, Q
end

"""

Computates the spline of the Thrust and Torque of each rotor in time
"""
function spline(t, T, torque)
i=0
   for (t in seq(tinit, tfinal, by=tstep))
   x[i]=t
   i++
   end
 for rotor in rotors
y[rotor] = T[rotor]
z[rotor] = Q[rotor]
splT[rotor] = Spline1D(x, y[rotor])
splQ[rotor] = Spline1D(x, z[rotor])
end
return splT, splQ
end

"""
equations used in order to obtain the second derivatives of the state of the vehicle,
the forces vector contains the force produced by each rotor and the weight of the
whole structure
"""
function equations(atm, MassProp, State, l, rotors, splQ, splT)


    ct, cp, cpsi = cos.([state.theta, state.phi, state.psi])
    st, sp, spsi = sin.([state.theta, state.phi, state.psi])

    theta2=((-splT[1]-splT[2]+splT[3]+splT[4])*l)/
    .Ixx
    phi2=((-splT[1]+splT[2]+splT[3]-splT[4])*l)/massprop.Iyy
    psi2=(splQ[1]-splQ[2]+splQ[3]-splQ[4])/massprop.Izz
    x2=sum(splT[i]*(sp*spsi+cp*cpsi*st) for i in 1:enumerate(rotors))/massprop.m
    y2=sum(splT[i]*(spsi*st*cp-cpsi*sp) for i in 1:enumerate(rotors))/massprop.m
    z2=sum(splT[i]*(cp*ct)-atm.g for i in 1:enumerate(rotors))/massprop.m
    acc=[theta2, phi2, psi2, x2, y2, z2]
    return acc
end

#function integrationv(acc,state)
 #   Pkg.add("SymPy")
#    using SymPy
#    theta1=integrate(acc[1], t)+state.q
#    phi1=integrate(acc[2],t)+state.p
#    psi1=integrate(acc[3],t)+state.r
#    x1=integrate(acc[4],t)+state.u
 #   y1=integrate(acc[5],t)+state.v
 #   z1=integrate(acc[6],t)+state.w
 #   vel=[theta1, phi1, psi1, x1, y1, z1]
#    return vel
#end

#function integrationp(vel,state)
 #   theta=integrate(vel[1],t)+state.theta
 #   phi=integrate(vel[2],t)+state.phi
#    psi=integrate(vel[3],t)+state.psi
#    x=integrate(vel[4],t)+state.x
#    y=integrate(vel[5],t)+state.y
#    z= integrate(vel[6],t)+state.z
#    pos=[theta,phi,psi,x,y,z]
#    return pos
#end

#function simulate
#	sim = Simulation()
#	for t in 1:100
#  	run(sim, t)
#  	forces[i][t], torque[i][t] = calcforce (tspan, rotors)
#	end
#end





function integrate(massprop,State, splQ, splT, l,n_motor, atm, tspan, s0, acc, Reference)
p = massprop,state, splQ, splT,l, enumerate(rotors), atm, Reference
prob = DifferentialEquations.ODEProblem(equations, s0, tspan, p)
sol = DifferentialEquations.solve(prob)

return sol

end

function integrate2(massprop,state,splQ, splT, l,n_motor, atm, tspan, s0, acc, Reference)
p = massprop,state, splQ, splT,l, enumerate(rotors), atm, Reference
prob = DifferentialEquations.SecondOrderODEProblem(equations, u0, s0, tspan, p)
sol = DifferentialEquations.solve(prob)
nothing # hide

return sol
end

