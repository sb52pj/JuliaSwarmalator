using Distributed
const numofprocess = 10
#const numofprocess = 2
addprocs(numofprocess)
@everywhere using DifferentialEquations
@everywhere using Logging: global_logger
@everywhere using TerminalLoggers: TerminalLogger
@everywhere include("./rhs.jl")
@everywhere using .rhs
@everywhere using Random
@everywhere using JLD2
using ProgressMeter
@everywhere global_logger(TerminalLogger())
@everywhere using Dates
@everywhere using NPZ
# @everywhere global_logger(TerminalLoggers.TerminalLogger())


#@everywhere const a, dt, T, n, L = 1.0, 0.5, 21000, 100, 1.0
@everywhere const a, dt, T, n, L = 1.0, 0.5, 6000, 100, 1.0
@everywhere const tspan = (0, T)  # Replace t_start and t_end with actual time span
#@everywhere const scanstep = 10
#@everywhere const scanstep = 30
@everywhere const scanstep = 90
@everywhere sideLength = size(30:scanstep:900)[1]
#@everywhere const caption = "Julia_correctview_a0a1_32"
#@everywhere const caption = "Julia_correctview_a0a1_$sideLength"*"α2_β3"
#@everywhere const caption = "Julia_correctview_a0a1_$sideLength"*"moveright"
#@everywhere const caption = "Julia_correctview_a0a1_$sideLength"*"AonPhaseR"
#@everywhere const caption = "Julia_correctview_a0a1_$sideLength"*"originaleuler"
#@everywhere const caption = "Julia_correctview_a0a1_smallvalue"
@everywhere const caption = "Julia_correctview_a0a1_randominit"
@everywhere const savepath = "F:\\xj\\experimentresult\\$caption-n_$n-T_$T\\"
if !isdir(savepath)
    mkdir(savepath)
    mkdir(savepath *"data\\")
end

@everywhere completed_tasks::Int64 = 0

@everywhere  struct Timer
	tstart::Array{Float64,1}
	runtime::Float64
	Timer(t) = new([0.0],t)
    end
@everywhere  function (timer::Timer)(u,t,integrator) 
	if timer.tstart[1] == 0.0 
		    timer.tstart[1] = time()
	end
	time() - timer.tstart[1] >= timer.runtime && return true
    end



@everywhere function callback_func(sol, t, integrator)
    current_stepsize = integrator.h
    threshold = 0.001
    if current_stepsize < threshold
        println("Small timestep detected! Aborting...")
        println("Parameters: ", integrator.p)
        terminate!(integrator)
    end
end
@everywhere function swarm(z)
    # Generate random initial conditions
    x0 = rand(Float64, n) .* (2 * L) .- L
    y0 = rand(Float64, n) .* (2 * L) .- L
    theta0 = rand(Float64, n) .* (2 * π) .- π
    omega = zeros(Float64, n)
    scalefactor = 0.1
    #omega = rand(Float64, n) .* (2 * scalefactor * L) .- L * scalefactor
    
    J, K, a0, a1, loghandle = z
    println("calculating $J $K $a0 $a1 for $T")
    # Combine initial conditions into a single vector
    z0 = vcat(x0, y0, theta0,copy(x0), copy(y0))

    # Initial conditions and parameters
    p = (J, K, n, omega, a0, a1, T)  
    #prob = ODEProblem(rhs_unit_vector!, z0, tspan, p)
    #prob = ODEProblem(rhs_unit_vector_recover_3D!, z0, tspan, p)
    #prob = ODEProblem(rhs_unit_vector_recover_lorenz!, z0, tspan, p)
    #prob = ODEProblem(rhs_unit_vector_recover_Exp!, z0, tspan, p)
    #prob = ODEProblem(rhs_unit_vector_recover!, z0, tspan, p)
    prob = ODEProblem(rhs_unit_vector_randinit!, z0, tspan, p)
    #prob = ODEProblem(rhs_unit_vector_AonPhase!, z0, tspan, p)
    #prob = ODEProblem(rhs_unit_vector_distance_α2_β3!, z0, tspan, p)
    timetest = 180
    timelimit = 7200
    function callback_func!(integrator)
      if T / integrator.t * timetest > timelimit
	terminated = true
	terminate!(integrator)
      end
    end
    affect!(integrator) = error("tlns:$z ") #terminate!(integrator)
    # 先把20分钟能解决的扫了
    cb = DiscreteCallback(Timer(timetest),callback_func!)
    @time sols = solve(prob, Vern9(), alg_hints=[:stiff], progress = true, callback =cb)
    #@time sols = solve(prob, Vern9(), alg_hints=[:stiff], progress = true)
    #@time sols = solve(prob, Midpoint(), alg_hints=[:stiff], progress = true)
    #@time sols = solve(prob, AutoVern9(Rodas5()), alg_hints=[:stiff], progress = true, callback=cb)

    #save("$savepath\\solution$J-$K-$a0-$a1.jld2", "sol", sols)
    filename = "J_$J" * "__K_$K" * "__a0_$a0" * "__a1_$a1.npy"
    savedt = 1
    x = hcat([sol[1:n] for sol in sols.u[1:savedt:end]]...)
    y = hcat([sol[n+1:2n] for sol in sols.u[1:savedt:end]]...)
    theta = hcat([sol[2n+1:3n] for sol in sols.u[1:savedt:end]]...)
    npzwrite("$savepath\\data\\x_$filename", x)
    npzwrite("$savepath\\data\\y_$filename", y )
    npzwrite("$savepath\\data\\theta_$filename", theta)
    npzwrite("$savepath\\data\\t_$filename", sols.t) 
    global completed_tasks += 1
    current_time = now()
    @info("$current_time: saved $J $K $a0 $a1 for $T, $completed_tasks/$total_tasks")
    logfile = open("logfile.txt", "a")
    write(logfile, "$current_time: saved $J $K $a0 $a1 for $T, $completed_tasks/$total_tasks")
    close(logfile)

    x=y=theta=nothing
    # garbage collect
    sols = nothing
    prob = nothing    
end
@everywhere const JKs = [(0.1, 1.0), (1.0, 0.0), (1.0, -0.1), (1.0, -0.75), (0.1, -1.0)]
#@everywhere const JKs = [(1.0, -0.03),(1.0, -0.06),(1.0, -0.09)]
#@everywhere const JKs = [(0.1, 1.0), (1.0, 0.0), (1.0, -0.1), (1.0, -0.75)]
#@everywhere const JKs = [(0.1, 1.0), (1.0, 0.0)]

#@everywhere const a0s = [round(deg2rad(ϕ),digits=2) for ϕ in 30:scanstep:120]
#@everywhere const a1s = [round(deg2rad(ϕ),digits=2) for ϕ in 30:scanstep:120]
@everywhere const a0s = append!([round(deg2rad(ϕ),digits=2) for ϕ in 30:scanstep:900], [100, 1000, 5000] )
@everywhere const a1s = append!([round(deg2rad(ϕ),digits=2) for ϕ in 30:scanstep:900], [100, 1000, 5000] )
#@everywhere const a0s = [1000]
#@everywhere const a1s = [1000]
@everywhere const loghandle = ""
@everywhere const allparams = [(JK[1], JK[2], a0, a1, loghandle) for JK in JKs, a0 in a0s, a1 in a1s]
#@everywhere const iterparams= [(JK[1], JK[2], a0, a1, loghandle) for JK in JKs, a0 in a0s, a1 in a1s]

# filtering done tasks
@everywhere const donefile = filter(x->startswith(x, "x_"), readdir("$savepath\\data") )
#@everywhere const xfiles3d = [ "x_J_$J" * "__K_$K" * "__a0_$a0" * "__a1_$a1.npy" for (J,K,a0,a1,loghandle) in allparams]
#@everywhere const xfiles = reduce(vcat, xfiles3d)
#@everywhere const iterlist = setdiff(xfiles,donefile)

@everywhere const iterparams = [ (J,K,a0,a1,loghandle) for (J,K,a0,a1,loghandle) in allparams if "x_J_$J" * "__K_$K" * "__a0_$a0" * "__a1_$a1.npy" ∉ donefile]

#@everywhere const total_tasks = size(iterparams)[1] * size(iterparams)[2] * size(iterparams)[3]
@everywhere const remain_tasks = size(iterparams)[1]
@everywhere const total_tasks = size(allparams)[1] * size(allparams)[2] * size(allparams)[3]
println("computing $remain_tasks/$total_tasks tasks in $numofprocess processes")


@sync @showprogress @distributed for z in iterparams
    #if isfile("$savepath\\solution$J-$K-$a0-$a1.jld2")
    J, K, a0, a1, loghandle = z
    filename = "J_$J" * "__K_$K" * "__a0_$a0" * "__a1_$a1.npy"
    if isfile("$savepath\\data\\x_$filename")
        global completed_tasks += 1
        @info("skipped $J $K $a0 $a1 for $T")
    end
    #println("calculating $J $K $a0 $a1 for $T")
    swarm(z)
    GC.gc(true)
end
