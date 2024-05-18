module ViewFunction
using LinearAlgebra

export I_t1k
export I_t1k_theta
export I_xc_vavg
export I_moveright
export I_xc_vavg_a0isϕ
export I_xc_vavg_a0isϕ_moveright
export I_xc_vavg_a0isϕ_thetaa0
export I_xc_bool
export I_3D_t1k
export I_3D
export I_lorenz
export I_Exp
export I_xc_bool_simple
export I_reverseVelocity

function I_t1k(func, x_current, y_current, x_former, y_former, a_0, t)
    if t > 1000
        return func(x_current, y_current, x_former, y_former, a_0)
    else
        return ones(size(x_current))
    end
end
function I_t1k_theta(func, x_current, y_current, x_former, y_former, theta_current, a_0, t)
    if t > 1000
        return func(x_current, y_current, x_former, y_former, theta_current, a_0)
    else
        return ones(size(x_current))
    end
end
function I_reverseVelocity(x_current, y_current, x_former, y_former, a_0)
    # xd = x_current .- x_current'
    # yd = y_current .- y_current'
    xd = x_current' .- x_current
    yd = y_current' .- y_current
    vx = x_current .- x_former
    vy = y_current .- y_former
    vx = -vx
    vy = -vy
    n = size(x_current)[1]

    #println("xd vx vy:",xd,vx,vy)
    dist = sqrt.(xd.^2 .+ yd.^2)
    vector_len = sqrt.(vx.^2 .+ vy.^2)
    scala_product_array = xd .* vx' .+ yd .* vy'
    inverse_cos_normalize = 1.0 ./ (dist .* vector_len')
    inverse_cos_normalize = replace(inverse_cos_normalize, NaN => 0.0)
        
    inverse_cos_normalize .*= (1.0 .- I(n))     
    products = clamp.(scala_product_array .* inverse_cos_normalize, -1.0, 1.0) # omit computational error resulting in DomainError 
    phi = acos.(products)
    #println(phi)
    # different a
    Ifun = exp.(-phi .* phi / a_0)
    #Ifun = exp.(-(phi .- π) .^ 2 / a_0)
    Ifun .*= (1.0 .- I(n) ) 
    Ifun = replace(Ifun, NaN => 0.0)
    return Ifun
end


function I_moveright(x_current, y_current, x_former, y_former, a_0)
    # xd = x_current .- x_current'
    # yd = y_current .- y_current'
    xd = x_current' .- x_current
    yd = y_current' .- y_current
    vx = x_current .- x_former
    vy = y_current .- y_former
    n = size(x_current)[1]

    #println("xd vx vy:",xd,vx,vy)
    dist = sqrt.(xd.^2 .+ yd.^2)
    vector_len = sqrt.(vx.^2 .+ vy.^2)
    scala_product_array = xd .* vx' .+ yd .* vy'
    inverse_cos_normalize = 1.0 ./ (dist .* vector_len')
        
    inverse_cos_normalize .*= (1.0 .- I(n))     
    products = clamp.(scala_product_array .* inverse_cos_normalize, -1.0, 1.0) # omit computational error resulting in DomainError 
    phi = acos.(products)
    #println(phi)
    # different a
    #Ifun = 1.05 .- exp.(-phi .* phi / a_0)
    #Ifun = exp.(-(phi .- 2*π) .^ 2 / a_0)
    Ifun = exp.(-(phi .- 2*π) .^ 2 / a_0)
    Ifun .*= (1.0 .- I(n) ) 
    Ifun = replace(Ifun, NaN => 0.0)
    return Ifun
end
function I_xc_vavg(x_current, y_current, x_former, y_former, a_0)
    # xd = x_current .- x_current'
    # yd = y_current .- y_current'
    xd = x_current' .- x_current
    yd = y_current' .- y_current
    vx = x_current .- x_former
    vy = y_current .- y_former
    n = size(x_current)[1]

    #println("xd vx vy:",xd,vx,vy)
    dist = sqrt.(xd.^2 .+ yd.^2)
    vector_len = sqrt.(vx.^2 .+ vy.^2)
    scala_product_array = xd .* vx' .+ yd .* vy'
    inverse_cos_normalize = 1.0 ./ (dist .* vector_len')
    inverse_cos_normalize = replace(inverse_cos_normalize, NaN => 0.0)
        
    inverse_cos_normalize .*= (1.0 .- I(n))     
    products = clamp.(scala_product_array .* inverse_cos_normalize, -1.0, 1.0) # omit computational error resulting in DomainError 
    phi = acos.(products)
    #println(phi)
    # different a
    Ifun = exp.(-phi .* phi / a_0)
    #Ifun = exp.(-(phi .- π) .^ 2 / a_0)
    Ifun .*= (1.0 .- I(n) ) 
    Ifun = replace(Ifun, NaN => 0.0)
    return Ifun
end
function I_lorenz(x_current, y_current, x_former, y_former, a_0)
    # xd = x_current .- x_current'
    # yd = y_current .- y_current'
    xd = x_current' .- x_current
    yd = y_current' .- y_current
    vx = x_current .- x_former
    vy = y_current .- y_former
    n = size(x_current)[1]

    #println("xd vx vy:",xd,vx,vy)
    dist = sqrt.(xd.^2 .+ yd.^2)
    vector_len = sqrt.(vx.^2 .+ vy.^2)
    scala_product_array = xd .* vx' .+ yd .* vy'
    inverse_cos_normalize = 1.0 ./ (dist .* vector_len')
        
    inverse_cos_normalize .*= (1.0 .- I(n))     
    products = clamp.(scala_product_array .* inverse_cos_normalize, -1.0, 1.0) # omit computational error resulting in DomainError 
    phi = acos.(products)
    #println(phi)
    # different a
    Ifun = 1 ./ ( 1 .+ phi .* phi ./ a_0)
    #Ifun = exp.(-(phi .- π) .^ 2 / a_0)
    Ifun .*= (1.0 .- I(n) ) 
    Ifun = replace(Ifun, NaN => 0.0)
    return Ifun
end

function I_Exp(x_current, y_current, x_former, y_former, a_0)
    # xd = x_current .- x_current'
    # yd = y_current .- y_current'
    xd = x_current' .- x_current
    yd = y_current' .- y_current
    vx = x_current .- x_former
    vy = y_current .- y_former
    n = size(x_current)[1]

    #println("xd vx vy:",xd,vx,vy)
    dist = sqrt.(xd.^2 .+ yd.^2)
    vector_len = sqrt.(vx.^2 .+ vy.^2)
    scala_product_array = xd .* vx' .+ yd .* vy'
    inverse_cos_normalize = 1.0 ./ (dist .* vector_len')
        
    inverse_cos_normalize .*= (1.0 .- I(n))     
    products = clamp.(scala_product_array .* inverse_cos_normalize, -1.0, 1.0) # omit computational error resulting in DomainError 
    phi = acos.(products)
    #println(phi)
    # different a
    Ifun = exp.(-phi / a_0)
    Ifun .*= (1.0 .- I(n) ) 
    Ifun = replace(Ifun, NaN => 0.0)
    return Ifun
end



function I_xc_bool(x_current, y_current, x_former, y_former, a_0)
    # xd = x_current .- x_current'
    # yd = y_current .- y_current'
    xd = x_current' .- x_current
    yd = y_current' .- y_current
    vx = x_current .- x_former
    vy = y_current .- y_former
    n = size(x_current)[1]

    #println("xd vx vy:",xd,vx,vy)
    dist = sqrt.(xd.^2 .+ yd.^2)
    vector_len = sqrt.(vx.^2 .+ vy.^2)
    scala_product_array = xd .* vx' .+ yd .* vy'
    inverse_cos_normalize = 1.0 ./ (dist .* vector_len')
        
    inverse_cos_normalize .*= (1.0 .- I(n))     
    products = clamp.(scala_product_array .* inverse_cos_normalize, -1.0, 1.0) # omit computational error resulting in DomainError 
    phi = acos.(products)
    #println(phi)
    # different a
    phi[ phi .> sqrt(a_0)] .= 0
    phi[ phi .> sqrt(a_0)] .= 1
    Ifun = phi
    #Ifun = exp.(-(phi .- π) .^ 2 / a_0)
    Ifun .*= (1.0 .- I(n) ) 
    Ifun = replace(Ifun, NaN => 0.0)
    return Ifun
end
function I_xc_bool_simple(x_current, y_current, x_former, y_former, a_0)
    # xd = x_current .- x_current'
    # yd = y_current .- y_current'
    xd = x_current' .- x_current
    yd = y_current' .- y_current
    vx = x_current .- x_former
    vy = y_current .- y_former
    n = size(x_current)[1]

    #println("xd vx vy:",xd,vx,vy)
    dist = sqrt.(xd.^2 .+ yd.^2)
    vector_len = sqrt.(vx.^2 .+ vy.^2)
    scala_product_array = xd .* vx' .+ yd .* vy'
    inverse_cos_normalize = 1.0 ./ (dist .* vector_len')
        
    inverse_cos_normalize .*= (1.0 .- I(n))     
    products = clamp.(scala_product_array .* inverse_cos_normalize, -1.0, 1.0) # omit computational error resulting in DomainError 
    phi = acos.(products)
    #println(phi)
    # different a
    phi[ phi .> (a_0)] .= 0
    phi[ phi .> (a_0)] .= 1
    Ifun = phi
    #Ifun = exp.(-(phi .- π) .^ 2 / a_0)
    Ifun .*= (1.0 .- I(n) ) 
    Ifun = replace(Ifun, NaN => 0.0)
    return Ifun
end


function I_xc_vavg_a0isϕ(x_current, y_current, x_former, y_former, a_0)
    # xd = x_current .- x_current'
    # yd = y_current .- y_current'
    xd = x_current' .- x_current
    yd = y_current' .- y_current
    vx = x_current .- x_former
    vy = y_current .- y_former
    n = size(x_current)[1]

    #println("xd vx vy:",xd,vx,vy)
    dist = sqrt.(xd.^2 .+ yd.^2)
    vector_len = sqrt.(vx.^2 .+ vy.^2)
    scala_product_array = xd .* vx' .+ yd .* vy'
    inverse_cos_normalize = 1.0 ./ (dist .* vector_len')
        
    inverse_cos_normalize .*= (1.0 .- I(n))     
    products = clamp.(scala_product_array .* inverse_cos_normalize, -1.0, 1.0) # omit computational error resulting in DomainError 
    phi = acos.(products) 
    #println(phi)
    # different a
    #Ifun = exp.(-phi .* phi / a_0)
    
    a_0 = mod.(a_0, 2π)
    Ifun = exp.(-(phi) .^ 2 ./ a_0)
    Ifun .*= (1.0 .- I(n) ) 
    Ifun = replace(Ifun, NaN => 0.0)
    return Ifun
end

function I_xc_vavg_a0isϕ_moveright(x_current, y_current, x_former, y_former, a_0)
    # xd = x_current .- x_current'
    # yd = y_current .- y_current'
    xd = x_current' .- x_current
    yd = y_current' .- y_current
    vx = x_current .- x_former
    vy = y_current .- y_former
    n = size(x_current)[1]

    #println("xd vx vy:",xd,vx,vy)
    dist = sqrt.(xd.^2 .+ yd.^2)
    vector_len = sqrt.(vx.^2 .+ vy.^2)
    scala_product_array = xd .* vx' .+ yd .* vy'
    inverse_cos_normalize = 1.0 ./ (dist .* vector_len')
        
    inverse_cos_normalize .*= (1.0 .- I(n))     
    products = clamp.(scala_product_array .* inverse_cos_normalize, -1.0, 1.0) # omit computational error resulting in DomainError 
    phi = acos.(products) 
    #println(phi)
    # different a
    #Ifun = exp.(-phi .* phi / a_0)
    
    a_0 = mod.(a_0, 2π)
    Ifun = exp.(-(phi) .^ 2 ./ a_0)
    Ifun .*= (1.0 .- I(n) ) 
    Ifun = replace(Ifun, NaN => 0.0)
    return Ifun
end


function I_xc_vavg_a0isϕ_thetaa0(x_current, y_current, x_former, y_former, theta_current, a_0)
    # xd = x_current .- x_current'
    # yd = y_current .- y_current'
    xd = x_current' .- x_current
    yd = y_current' .- y_current
    vx = x_current .- x_former
    vy = y_current .- y_former
    n = size(x_current)[1]

    #println("xd vx vy:",xd,vx,vy)
    dist = sqrt.(xd.^2 .+ yd.^2)
    vector_len = sqrt.(vx.^2 .+ vy.^2)
    scala_product_array = xd .* vx' .+ yd .* vy'
    inverse_cos_normalize = 1.0 ./ (dist .* vector_len')
        
    inverse_cos_normalize .*= (1.0 .- I(n))     
    products = clamp.(scala_product_array .* inverse_cos_normalize, -1.0, 1.0) # omit computational error resulting in DomainError 
    phi = acos.(products) 
    #println(phi)
    # different a
    #Ifun = exp.(-phi .* phi / a_0)
    
    a_0 = mod.(theta_current, 2π) .+ a_0
    Ifun = exp.(-(phi) .^ 2 ./ a_0)
    Ifun .*= (1.0 .- I(n) ) 
    Ifun = replace(Ifun, NaN => 0.0)
    return Ifun
end


function I_3D_t1k(func, x_current,y_current,z_current, x_former,y_former,z_former, a_0, t)
    if t > 1000
        return func(x_current, y_current,z_current, x_former, y_former, z_former, a_0)
    else
        return ones(size(x_current))
    end
end

function I_3D(x_current, y_current,z_current, x_former, y_former, z_former, a_0)
    xd = x_current' .- x_current
    yd = y_current' .- y_current
    zd = z_current' .- z_current
    vx = x_current .- x_former
    vy = y_current .- y_former
    vz = z_current .- z_former
    n = size(x_current)[1]

    #println("xd vx vy:",xd,vx,vy)
    dist = sqrt.(xd.^2 .+ yd.^2 + zd.^2)
    vector_len = sqrt.(vx.^2 .+ vy.^2 + vz .^2)
    scala_product_array = xd .* vx' .+ yd .* vy' + zd .* vz'
    inverse_cos_normalize = 1.0 ./ (dist .* vector_len')

    inverse_cos_normalize .*= (1.0 .- I(n))     
    products = clamp.(scala_product_array .* inverse_cos_normalize, -1.0, 1.0) # omit computational error resulting in DomainError 
    phi = acos.(products)
    #println(phi)
    # different a
    Ifun = exp.(-phi .* phi / a_0)
    #Ifun = exp.(-(phi .- π) .^ 2 / a_0)
    Ifun .*= (1.0 .- I(n) ) 
    Ifun = replace(Ifun, NaN => 0.0)
    return Ifun

end

end

