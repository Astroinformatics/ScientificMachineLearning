#=
    ODE models for orbital mechanics
=#

function NewtonianOrbitModel(u, model_params, t)
    #=
        Defines system of odes which describes motion of
        point like particle with Newtonian physics, uses

        u[1] = χ
        u[2] = ϕ

        where, p, M, and e are constants
    =#
    χ, ϕ = u
    p, M, e  = model_params

    numer = (1+e*cos(χ))^2
    denom = M*(p^(3/2))

    χ̇ = numer / denom
    ϕ̇ = numer / denom

    return [χ̇, ϕ̇]

end

function RelativisticOrbitModel(u, model_params, t)
    #=
        Defines system of odes which describes motion of
        point like particle in schwarzschild background, uses

        u[1] = χ
        u[2] = ϕ

        where, p, M, and e are constants
    =#
    χ, ϕ = u
    p, M, e  = model_params

    numer = (p-2-2*e*cos(χ)) * (1+e*cos(χ))^2
    denom = sqrt( (p-2)^2-4*e^2 )

    χ̇ = numer * sqrt( p-6-2*e*cos(χ) )/( M*(p^2)*denom )
    ϕ̇ = numer / (M*(p^(3/2))*denom)

    return [χ̇, ϕ̇]

end

function AbstractNNOrbitModel(u, model_params, t; NN=nothing, NN_params=nothing)
    #=
        Defines system of odes which describes motion of
        point like particle with Newtonian physics, uses

        u[1] = χ
        u[2] = ϕ

        where, p, M, and e are constants
    =#
    χ, ϕ = u
    p, M, e  = model_params

    if isnothing(NN)
        nn = [1,1]
    else
        nn = 1 .+ NN(u, NN_params)
    end

    numer = (1 + e*cos(χ))^2
    denom = M*(p^(3/2))

    χ̇ = (numer / denom) * nn[1]
    ϕ̇ = (numer / denom) * nn[2]

    return [χ̇, ϕ̇]

end

function AbstractNROrbitModel(u, model_params, t;
                              NN_chiphi=nothing, NN_chiphi_params=nothing,
                              NN_pe=nothing, NN_pe_params=nothing)
    #=
        Defines system of odes which describes motion of
        point like particle with Newtonian physics, uses

        u[1] = χ
        u[2] = ϕ
        u[3] = p
        u[4] = e

        q is the mass ratio
    =#
    χ, ϕ, p, e = u
    q = model_params[1]
    M=1.0

    if p <= 0
        println("p = ", p)
    end

    if isnothing(NN_chiphi)
        nn_chiphi = [1,1]
    else
        nn_chiphi = 1 .+ NN_chiphi(u, NN_chiphi_params)
    end

    if isnothing(NN_pe)
        nn_pe = [0,0]
    else
        nn_pe = NN_pe(u, NN_pe_params)
    end

    numer = (1+e*cos(χ))^2
    denom = M*(abs(p)^(3/2))

    χ̇ = (numer / denom) * nn_chiphi[1]
    ϕ̇ = (numer / denom) * nn_chiphi[2]
    ṗ = nn_pe[1]
    ė = nn_pe[2]

    return [χ̇, ϕ̇, ṗ, ė]
end
