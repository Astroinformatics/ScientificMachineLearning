#=
    Axiliary functions for orbital mechanics
=#
using DelimitedFiles

function soln2orbit(soln, model_params=nothing)
    #=
        Performs change of variables:
        (χ(t),ϕ(t)) ↦ (x(t),y(t))
    =#
    if size(soln,1) == 2
        χ = soln[1,:]
        ϕ = soln[2,:]
        if length(model_params) == 3
            p, M, e  = model_params
        else
            error("model_params must have length 3 when size(soln,2) = 2")
        end
    elseif size(soln,1) == 4
        χ = soln[1,:]
        ϕ = soln[2,:]
        p = soln[3,:]
        e = soln[4,:]
    else
        error("size(soln,2) must be either 2 or 4")
    end

    r = p ./ (1 .+ e .* cos.(χ))
    x = r .* cos.(ϕ)
    y = r .* sin.(ϕ)

    orbit = vcat(x',y')
    return orbit
end

function orbit2tensor(orbit, component, mass=1.0)
    #=
        Construct trace-free moment tensor Ι(t) for orbit from BH orbit (x(t),y(t))

        component defines the Cartesion indices in x,y. For example,
        I_{22} is the yy component of the moment tensor.
    =#
    x = orbit[1,:]
    y = orbit[2,:]

    Ixx = x.^2
    Iyy = y.^2
    Ixy = x.*y
    trace = Ixx .+ Iyy

    if component[1] == 1 && component[2] == 1
        tmp = Ixx .- (1.0 ./ 3.0).*trace
    elseif component[1] == 2 && component[2] == 2
        tmp = Iyy .- (1.0 ./ 3.0).*trace
    else
        tmp = Ixy
    end

    return mass .* tmp
end

function d_dt(v::AbstractVector, dt)
    # uses second-order one-sided difference stencils at the endpoints; see https://doi.org/10.1090/S0025-5718-1988-0935077-0
    a = -3/2*v[1] + 2*v[2] - 1/2*v[3]
    b = (v[3:end] .- v[1:end-2])/2
    c = 3/2*v[end] - 2*v[end-1] + 1/2*v[end-2]
    return [a;b;c]/dt
end

function d2_dt2(v::AbstractVector, dt)
    # uses second-order one-sided difference stencils at the endpoints; see https://doi.org/10.1090/S0025-5718-1988-0935077-0
    a = 2*v[1] - 5*v[2] + 4*v[3] - v[4]
    b = v[1:end-2] .- 2*v[2:end-1] .+ v[3:end]
    c = 2*v[end] - 5*v[end-1] + 4*v[end-2] - v[end-3]
    return [a;b;c]/(dt^2)
end

function h_22_quadrupole_components(dt, orbit, component, mass=1.0)
    #=
        x(t) and y(t) inputs are the trajectory of the orbiting BH.

       WARNING: assuming x and y are on a uniform grid of spacing dt
        x_index and y_index are 1,2,3 for x, y, and z indicies.
    =#

    mtensor = orbit2tensor(orbit, component, mass)
    mtensor_ddot = d2_dt2(mtensor,dt)

    # return mtensor
    return 2*mtensor_ddot
end

function h_22_quadrupole(dt, orbit, mass=1.0)
    h11 = h_22_quadrupole_components(dt, orbit, (1,1), mass)
    h22 = h_22_quadrupole_components(dt, orbit, (2,2), mass)
    h12 = h_22_quadrupole_components(dt, orbit, (1,2), mass)
    return h11, h12, h22
end

function h_22_strain_one_body(dt, orbit)

    h11, h12, h22 = h_22_quadrupole(dt, orbit)

    h₊ = h11 - h22
    hₓ = 2.0*h12

    scaling_const = sqrt(pi/5)
    return scaling_const*h₊, -scaling_const*hₓ
end

function h_22_quadrupole_two_body(dt, orbit1, mass1, orbit2, mass2)
    h11_1, h12_1, h22_1 = h_22_quadrupole(dt, orbit1, mass1)
    h11_2, h12_2, h22_2 = h_22_quadrupole(dt, orbit2, mass2)
    h11 = h11_1 + h11_2
    h12 = h12_1 + h12_2
    h22 = h22_1 + h22_2
    return h11, h12, h22
end

function h_22_strain_two_body(dt, orbit1, mass1, orbit2, mass2)
    # compute (2,2) mode strain from orbits of BH 1 of mass1 and BH2 of mass 2

    @assert abs(mass1 + mass2 - 1.0) < 1e-12 "Masses do not sum to unity"

    h11, h12, h22 = h_22_quadrupole_two_body(dt, orbit1, mass1, orbit2, mass2)

    h₊ = h11 - h22
    hₓ = 2.0*h12

    scaling_const = sqrt(pi/5)
    return scaling_const*h₊, -scaling_const*hₓ
end

function one2two(path, m1, m2)
    #=
        We need a very crude 2-body path

        Assume the 1-body motion is a newtonian 2-body position vector r = r1 - r2
        and use Newtonian formulas to get r1, r2
        (e.g. Theoretical Mechanics of Particles and Continua 4.3)
    =#

    M = m1 + m2
    r1 = m2/M .* path
    r2 = -m1/M .* path

    return r1, r2
end

function compute_waveform(dt, soln, mass_ratio, model_params=nothing)

    @assert mass_ratio <= 1.0 "mass_ratio must be <= 1"
    @assert mass_ratio >= 0.0 "mass_ratio must be non-negative"

    orbit = soln2orbit(soln, model_params)
    if mass_ratio > 0
        mass1 = mass_ratio/(1.0+mass_ratio)
        mass2 = 1.0/(1.0+mass_ratio)

        orbit1, orbit2 = one2two(orbit, mass1, mass2)
        waveform = h_22_strain_two_body(dt, orbit1, mass1, orbit2, mass2)
    else
        waveform = h_22_strain_one_body(dt, orbit)
    end
    return waveform
end

function interpolate_time_series(tsteps, tdata, fdata)

    @assert length(tdata) == length(fdata) "lengths of tdata and fdata must match"

    interp_fdata = zeros(length(tsteps))
    for j=1:length(tsteps)
        for i=1:length(tdata)-1
            if  tdata[i] <= tsteps[j] < tdata[i+1]
                weight = (tsteps[j] - tdata[i]) / (tdata[i+1] - tdata[i])
                interp_fdata[j] = (1-weight)*fdata[i] + weight*fdata[i+1]
                break
            end
        end
    end

    return interp_fdata
end

function file2waveform(tsteps, filename="waveform.txt")

    # read in file
    f = open(filename, "r")
    data = readdlm(f)
    tdata = data[:,1]
    wdata = data[:,2]

    # interpolate data to tsteps
    waveform = interpolate_time_series(tsteps, tdata, wdata)

    return waveform
end

function file2trajectory(tsteps, filename="trajectoryA.txt")

    # read in file
    f = open(filename, "r")
    data = readdlm(f)
    tdata = data[:,1]
    xdata = data[:,2]
    ydata = data[:,3]

    # interpolate data to tsteps
    x = interpolate_time_series(tsteps, tdata, xdata)
    y = interpolate_time_series(tsteps, tdata, ydata)

    return x, y
end
