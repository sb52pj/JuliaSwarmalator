module rhs
using LinearAlgebra
include("./ViewFunction.jl")
using .ViewFunction

export rhs_unit_vector!
export rhs_unit_vector_noa0a1!
export rhs_unit_vector_AonPhase!
export rhs_unit_vector_recover_nophase!
export rhs_unit_vector_recover_noviewonphase!
export rhs_unit_vector_distance_α2_β3!
export rhs_unit_vector_VoA_a0isϕ!
export rhs_unit_vector_bool!
export rhs_unit_vector_recover!
export rhs_unit_vector_recover_3D!
export rhs_unit_vector_recover_lorenz!
export rhs_unit_vector_recover_Exp!
export rhs_unit_vector_distance_α0_β2!
export rhs_unit_vector_moveright!
export rhs_unit_vector_original!
export rhs_unit_vector_randinit!
export rhs_unit_vector_VoA!
export rhs_unit_vector_VoR!
export rhs_unit_vector_recover_reverse!

function rhs_unit_vector_distance_α2_β3!(dz, z, p, t)
    J, K, n, omega, a_0, a_1, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_xc_vavg
    VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)
    VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_1, t)

    α = 2
    β = 3

    x_rhs = - xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .^ α .- VoR .* inverse_dist .^ β)
    y_rhs = - yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .^ α .- VoR .* inverse_dist .^ β)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end

function rhs_unit_vector_distance_α0_β2!(dz, z, p, t)
    J, K, n, omega, a_0, a_1, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_xc_vavg
    VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)
    VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_1, t)

    α = 0
    β = 2

    x_rhs = - xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .^ α .- VoR .* inverse_dist .^ β)
    y_rhs = - yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .^ α .- VoR .* inverse_dist .^ β)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end



function rhs_unit_vector_VoA_a0isϕ!(dz, z, p, t)
    J, K, n, omega, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_xc_vavg_a0isϕ
    VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, theta, t)
    VoR = 1

    x_rhs = - xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs = - yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    # VoA here is not necessary
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end
function rhs_unit_vector_AonPhase!(dz, z, p, t)
    J, K, n, omega, a_0, a_1, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    #viewfun = I_xc_vavg
    viewfun = I_xc_vavg_a0isϕ_thetaa0
    #VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)
    #VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_1, t)
    Vo1 = 1
    VoA = I_t1k_theta(viewfun, x,y, xbuffer, ybuffer, theta, a_0, t)
    VoR = I_t1k_theta(viewfun, x,y, xbuffer, ybuffer, theta, a_1, t)

    x_rhs = - xd .* (Vo1 .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs = - yd .* (Vo1 .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2) 
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end
function rhs_unit_vector_noa0a1!(dz, z, p, t)
    #println(t)
    #@debug "Inside my_ode_function!"
    J, K, n, omega, T = p
    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]
    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    x_rhs = - xd .* ((1 .+ J .* cos.(theta_d)) .* inverse_dist .- inverse_dist_sq)
    y_rhs = - yd .* ((1 .+ J .* cos.(theta_d)) .* inverse_dist .- inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist 

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)

    return nothing
end



function rhs_unit_vector!(dz, z, p, t)
    J, K, n, omega, a_0, a_1, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x' .- x
    yd = y' .- y
    theta_d = theta' .- theta

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_xc_vavg
    #viewfun = I_xc_vavg_a0isϕ_thetaa0
    VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)
    VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_1, t)
    #VoA = I_t1k_theta(viewfun, x,y, xbuffer, ybuffer, theta, a_0, t)
    #VoR = I_t1k_theta(viewfun, x,y, xbuffer, ybuffer, theta, a_1, t)

    x_rhs =  xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs =  yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    theta_rhs = K * sin.(theta_d) .* inverse_dist .* VoA

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end
function rhs_unit_vector_recover_nophase!(dz, z, p, t)
    J, K, n, a_0, T = p

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    xbuffer = z[2n+1:3n]
    ybuffer = z[3n+1:4n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq = replace(inverse_dist_sq, NaN => 0.0)

    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_xc_vavg
    VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)

    x_rhs = - xd .* (1 .+  J .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs = - yd .* (1 .+  J .* inverse_dist .- VoR .* inverse_dist_sq)

    x_rhs = replace(x_rhs, NaN => 0.0)
    y_rhs = replace(y_rhs, NaN => 0.0)

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= x
    dz[3n+1:4n] .= y
    return nothing
end

function rhs_unit_vector_recover_noviewonphase!(dz, z, p, t)
    J, K, n, omega, a_0, a_1, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq = replace(inverse_dist_sq, NaN => 0.0)

    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_xc_vavg
    VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)
    VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_1, t)

    x_rhs = - xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs = - yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist 

    x_rhs = replace(x_rhs, NaN => 0.0)
    y_rhs = replace(y_rhs, NaN => 0.0)

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end
function rhs_unit_vector_recover_3D!(dz, zipvar, p, t)
  # undone
    J, K, n, omega, a_0, a_1, T = p

    # Instantiate -- set up
    x = zipvar[1:n]
    y = zipvar[n+1:2n]
    z = zipvar[2n+1:3n]
    theta = zipvar[3n+1:4n]

    xbuffer = zipvar[4n+1:5n]
    ybuffer = zipvar[5n+1:6n]
    zbuffer = zipvar[6n+1:7n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    zd = z .- z'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2 + (zd) .^ 2)
    inverse_dist_sq = replace(inverse_dist_sq, NaN => 0.0)

    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_3D
    VoA = I_3D_t1k(viewfun, x,y,z, xbuffer, ybuffer, zbuffer, a_0, t)
    VoR = I_3D_t1k(viewfun, x,y,z, xbuffer, ybuffer, zbuffer, a_1, t)

    x_rhs = - xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs = - yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    z_rhs = - zd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    x_rhs = replace(x_rhs, NaN => 0.0)
    y_rhs = replace(y_rhs, NaN => 0.0)
    z_rhs = replace(z_rhs, NaN => 0.0)

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= (1/n) * sum((1.0 .-  I(n)) .* z_rhs, dims=2)
    dz[3n+1:4n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[4n+1:5n] .= x
    dz[5n+1:6n] .= y
    dz[6n+1:7n] .= z

    return nothing
end

function  rhs_unit_vector_recover_Exp!(dz, z, p, t)
    J, K, n, omega, a_0, a_1, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq = replace(inverse_dist_sq, NaN => 0.0)

    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_Exp
    VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)
    VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_1, t)

    x_rhs = - xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs = - yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    x_rhs = replace(x_rhs, NaN => 0.0)
    y_rhs = replace(y_rhs, NaN => 0.0)

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end
function  rhs_unit_vector_recover_reverse!(dz, z, p, t)
    J, K, n, omega, a_0, a_1, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq = replace(inverse_dist_sq, NaN => 0.0)

    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_reverseVelocity
    VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)
    VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_1, t)

    x_rhs = - xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs = - yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    x_rhs = replace(x_rhs, NaN => 0.0)
    y_rhs = replace(y_rhs, NaN => 0.0)

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end



function  rhs_unit_vector_recover_lorenz!(dz, z, p, t)
    J, K, n, omega, a_0, a_1, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq = replace(inverse_dist_sq, NaN => 0.0)

    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_lorenz
    VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)
    VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_1, t)

    x_rhs = - xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs = - yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    x_rhs = replace(x_rhs, NaN => 0.0)
    y_rhs = replace(y_rhs, NaN => 0.0)

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end

function rhs_unit_vector_moveright!(dz, z, p, t)
    J, K, n, omega, a_0, a_1, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq = replace(inverse_dist_sq, NaN => 0.0)

    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_moveright
    VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)
    VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_1, t)

    x_rhs = - xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs = - yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    x_rhs = replace(x_rhs, NaN => 0.0)
    y_rhs = replace(y_rhs, NaN => 0.0)

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end


function rhs_unit_vector_original!(dz, z, p, t)
    J, K, n, omega, a_0, a_1, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq = replace(inverse_dist_sq, NaN => 0.0)

    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_xc_vavg
    VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)
    VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_1, t)

    x_rhs = - xd .* ((1 .+  J .* cos.(theta_d)) .* inverse_dist .-  inverse_dist_sq)
    y_rhs = - yd .* ((1 .+  J .* cos.(theta_d)) .* inverse_dist .-  inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    x_rhs = replace(x_rhs, NaN => 0.0)
    y_rhs = replace(y_rhs, NaN => 0.0)

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end
function rhs_unit_vector_randinit!(dz, z, p, t)
    J, K, n, omega, a_0, a_1, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq = replace(inverse_dist_sq, NaN => 0.0)

    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_xc_vavg
    VoA = viewfun( x,y, xbuffer, ybuffer, a_0)
    VoR = viewfun( x,y, xbuffer, ybuffer, a_1)

    x_rhs = - xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs = - yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    x_rhs = replace(x_rhs, NaN => 0.0)
    y_rhs = replace(y_rhs, NaN => 0.0)

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end

function rhs_unit_vector_VoA!(dz, z, p, t)
    J, K, n, omega, a_0,  T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq = replace(inverse_dist_sq, NaN => 0.0)

    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_xc_vavg
    VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)

    x_rhs = - xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .-  inverse_dist_sq)
    y_rhs = - yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .-  inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    x_rhs = replace(x_rhs, NaN => 0.0)
    y_rhs = replace(y_rhs, NaN => 0.0)

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end


function rhs_unit_vector_VoR!(dz, z, p, t)
    J, K, n, omega, a_0, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq = replace(inverse_dist_sq, NaN => 0.0)

    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    viewfun = I_xc_vavg
    #VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)
    VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)

    x_rhs = - xd .* ((1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs = - yd .* ((1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoR

    x_rhs = replace(x_rhs, NaN => 0.0)
    y_rhs = replace(y_rhs, NaN => 0.0)

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end





function rhs_unit_vector_bool!(dz, z, p, t)
    J, K, n, omega, a_0, a_1, T = p

    if (t + 1) % 2000 == 0
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
        prinln("%d has been down for %d %d %d, in %s\n", t, J, K, a_0, Dates.format(now(), "HH:MM:SS yy-mm-dd"))
    end

    # Instantiate -- set up
    x = z[1:n]
    y = z[n+1:2n]
    theta = z[2n+1:3n]

    xbuffer = z[3n+1:4n]
    ybuffer = z[4n+1:5n]

    # Set up as a matrix to make the computation faster
    xd = x .- x'
    yd = y .- y'
    theta_d = theta .- theta'

    inverse_dist_sq = 1.0 ./ ((xd).^2 .+ (yd).^2)
    inverse_dist_sq[1:n+1:end] .= 0.0  # correct 1 / d_ii = 1 / 0 error
    inverse_dist = sqrt.(inverse_dist_sq)

    #viewfun = I_xc_vavg
    #viewfun = I_xc_bool
    viewfun = I_xc_bool_simple
    #VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_0, t)
    #VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer, a_1, t)
    VoA = I_t1k(viewfun, x,y, xbuffer, ybuffer,  a_0, t)
    VoR = I_t1k(viewfun, x,y, xbuffer, ybuffer,  a_1, t)

    x_rhs = - xd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    y_rhs = - yd .* (VoA .*(1 .+  J .* cos.(theta_d)) .* inverse_dist .- VoR .* inverse_dist_sq)
    theta_rhs = -K * sin.(theta_d) .* inverse_dist .* VoA

    # The actual R.H.S.
    dz[1:n] .= (1/n) * sum((1.0 .- I(n)) .* x_rhs, dims=2)
    dz[n+1:2n] .= (1/n) * sum((1.0 .-  I(n)) .* y_rhs, dims=2)
    dz[2n+1:3n] .= omega .+ (1/n) * sum((1.0 .-  I(n)) .* theta_rhs, dims=2)
    dz[3n+1:4n] .= x
    dz[4n+1:5n] .= y

    return nothing
end

end
