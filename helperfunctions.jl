using LazySets, LinearAlgebra

using Random # For strategy 4

function rectangleFromHBox(res::AbstractVector{Shape}, corners)
    for i in 1:size(corners,1)
        tope = getindex(corners, i)
        res[i] = Shape(getindex.(tope, 1), getindex.(tope, 2))
    end
    return res
end

function rectangleFromHBox!(res::AbstractVector{Shape}, cornerss, timestepsize, dim, timescale::Float64)
    deltt = timestepsize
    currentTime = 0
    for i in 1:size(cornerss, 1)#ændrer til getindex -> max min 
        tope = getindex(cornerss, i)
        dimCoords = getindex.(tope, dim)
        maxcor = maximum(dimCoords)
        mincor = minimum(dimCoords)
        res[i] = Shape([currentTime, currentTime + deltt, currentTime + deltt, currentTime], [mincor, mincor, maxcor, maxcor])
        #res[i] = Shape([deltt*(i-1), deltt*i, deltt*i, deltt*(i-1)], [mincor, mincor, maxcor, maxcor])
        currentTime = currentTime + deltt
        deltt = deltt+timescale*timestepsize
    end
    return res
end

function rectangleFromHBox2Dims(res::AbstractVector{Shape}, corners, dim1, dim2)
    for i in 1:size(corners, 1)
        tope = getindex(corners, i)
        # Switch 3rd and 4th index
        tope[3], tope[4] = tope[4], tope[3]

        res[i] = Shape(getindex.(tope, dim1), getindex.(tope, dim2))
    end
end

function rectangleFromHBoxWithTimestepArray(res::AbstractVector{Shape}, cornerss, timesteps, startTime, dim)

    currentTime = startTime
    for i in 1:size(cornerss, 1)#
        timestep = timesteps[i]

        tope = getindex(cornerss, i)
        dimCoords = getindex.(tope, dim)
        maxcor = maximum(dimCoords)
        mincor = minimum(dimCoords)
        res[i] = Shape([currentTime, currentTime + timestep, currentTime + timestep, currentTime], [mincor, mincor, maxcor, maxcor])
        #res[i] = Shape([deltt*(i-1), deltt*i, deltt*i, deltt*(i-1)], [mincor, mincor, maxcor, maxcor])
        currentTime = currentTime + timestep
    end
    return res
end

function rectangleFromHBox!(res::AbstractVector{Shape}, cornerss, timestepsize, dim, timesteps)

    for i in 1:size(cornerss, 1)#ændrer til getindex -> max min 
        tope = getindex(cornerss, i)
        dimCoords = getindex.(tope, dim)
        maxcor = maximum(dimCoords)
        mincor = minimum(dimCoords)
        curtime = timesteps[i]
        res[i] = Shape([timestepsize*(curtime-1), timestepsize*curtime, timestepsize*curtime, timestepsize*(curtime-1)], [mincor, mincor, maxcor, maxcor])
    end
    return res
end

function reachsets(A, timestepsize, interval, X₀, μ, timescale)
    T = maximum(interval)
    N = floor(Int, T/timestepsize)
    
    ANorm = norm(A, Inf)
    α = (exp(ANorm*timestepsize)-1-timestepsize*ANorm)/norm(X₀, Inf)
    β = (exp(ANorm*timestepsize)-1)*μ/ANorm

    ϕ = exp(A*timestepsize)

    ϕp = (I+ϕ)/2
    ϕm = (I-ϕ)/2
    gens = hcat(ϕp*X₀.generators,ϕm*X₀.center, ϕm*X₀.generators)

    R₁ = minkowski_sum(Zonotope(ϕp*X₀.center, gens), Zonotope(zeros(2), (α+β)*I(2)))
    #res = [R₁]
    R = [R₁]

    ballβ = Zonotope(zeros(2), β*I(2))
    i = 2
    deltt=timestepsize
    time = timestepsize # We have already done one timestep
    #while i < N
    while time < T
        deltt = deltt+timescale*timestepsize
        #α = (exp(ANorm*timestepsize)-1-timestepsize*ANorm)/norm(X₀, Inf)
        β = (exp(ANorm*deltt)-1)*μ/ANorm

        ϕ = exp(A*deltt)

        ballβ = Zonotope(zeros(dim(X₀)), β*I(2))
        push!(R, minkowski_sum(linear_map(ϕ, R[i-1]), ballβ))
        #res[i] = R[i]



        #N = floor(Int, T/deltt) # Re-estimate amount of timesteps
        time = time + deltt
        i = i + 1
    end
    return R
end

function initialStep(A, ANorm, timestep, X₀, μ)
    α = (exp(ANorm*timestep)-1-timestep*ANorm)/norm(X₀, Inf)
    β = (exp(ANorm*timestep)-1)*μ/ANorm

    ϕ = exp(A*timestep)

    ϕp = (I+ϕ)/2
    ϕm = (I-ϕ)/2
    gens = hcat(ϕp*X₀.generators,ϕm*X₀.center, ϕm*X₀.generators)


    R = minkowski_sum(Zonotope(ϕp*X₀.center, gens), Zonotope(zeros(dim(X₀)), (α+β)*I(dim(X₀))))

    ballβ = Zonotope(zeros(dim(X₀)), β*I(dim(X₀)))

    return (R, ballβ, ϕ)
end

function initialStepNoInput(A, ANorm, timestep, X₀)
    α = (exp(ANorm*timestep)-1-timestep*ANorm)/norm(X₀, Inf)
    ϕ = exp(A*timestep)

    ϕp = (I+ϕ)/2
    ϕm = (I-ϕ)/2
    gens = hcat(ϕp*X₀.generators,ϕm*X₀.center, ϕm*X₀.generators)


    R = minkowski_sum(Zonotope(ϕp*X₀.center, gens), Zonotope(zeros(dim(X₀)), α*I(dim(X₀))))

    return (R, ϕ)
end


function reachSetsForTimesteps(A, timesteps, interval, X₀, μ)
    ANorm = norm(A, Inf)

    timestep = timesteps[1] # Use first element

    R₁, ballβ, ϕ = initialStep(A, ANorm, timestep, X₀, μ)

    
    R = [R₁]

    # We only further utilize these if we have an input
    Ω = []
    S = ballβ
    V = ballβ
    if (μ != 0)
        push!(Ω, minkowski_sum(R₁, ballβ))
        V = linear_map(ϕ, V)
    end

    totalLength = length(timesteps)
    currentTime = timestep

    for i in 2:totalLength
        timestep = timesteps[i]
        if timestep == timesteps[i-1] # Check if timestep is the same
            prevR = R[i-1] # Reuse previous result

            newR = linear_map(ϕ, prevR)
            push!(R, newR)

            # If we have an input
            if (μ != 0)
                S = minkowski_sum(S, V)
                V = linear_map(ϕ, V)

                push!(Ω, minkowski_sum(newR, S))
            end

        else # Recalculate values for new timestep
            newR, ballβ, ϕ = initialStep(A, ANorm, timestep, X₀, μ)
            
            # We move R to the current time
            ϕCurrentTime = exp(A*currentTime)
            newR = Zonotope(ϕCurrentTime * newR.center, ϕCurrentTime * newR.generators)

            push!(R, newR)

            if (μ != 0)
                # We have to retroactivily calculate the new S
                steps = floor(Int, currentTime / timestep)
                S = ballβ
                V = ballβ
                for j in 1:steps
                    S = minkowski_sum(S, V)
                    V = linear_map(ϕ, V)
                end
                
                push!(Ω, minkowski_sum(newR, S))
            end
        end
        currentTime = currentTime + timestep
    end
    # If no input just return R
    if (μ == 0)
        return R
    end
    # Else return Ω
    return Ω
end

# Move a zonotope forward in time, assumes no input
function forwardTimeNoInput(A, R, currentTime)
    ϕCurrentTime = exp(A*currentTime)
    return Zonotope(ϕCurrentTime * R.center, ϕCurrentTime * R.generators)
end

function forwardTime(A, R, currentTime, timestep, ballβ, ϕ)
    newR = forwardTimeNoInput(A, R, currentTime)

    # Note that steps should always work
    if currentTime % timestep != 0
        #throw(DomainError(currentTime % timestep, "Error in forwardTime, $currentTime % $timestep != 0"))
        println("Error in forwardTime, $currentTime % $timestep != 0")
    end

    steps = floor(Int, currentTime / timestep)
    S = ballβ
    V = ballβ
    for j in 1:steps
        S = minkowski_sum(S, V)
        V = linear_map(ϕ, V)
    end
    newΩ = minkowski_sum(newR, S)

    #println("Set S, V: $S, $V with steps $steps")
    return (newR, S, V, newΩ)
end

# CEGAR like approach
# Note that we assume no input

# We have different strategies for how to choose the timestep
# 1: We half the timestpe on fail, and reset each time

# 2: We half the timestep on fail. Upon consecutive fails we take 2^(x-1) - 1 steps before attempting reset, 
# where x is the number of consecutive fails. We reset to initial timestep

# 3: Never reset. If no errors in past x steps, double timestep

# 4: Randomly choose to retain or increase
# Adjust chance of increase based on previous success of increase
# Keep different profiles dependt on amount of recent fails
# Specifically 0-1, 2-3, 4-5, 6+


const DIGITS = 4

function reachSetsCegar(A, initialTimestep, interval, X0, constraint, μ, strategy = 1)
    ANorm = norm(A, Inf)
    R = []
    Ω = []
    changedTimeStep = true # Keepss track of whether timestep has changed

    # Dictionary to keep track of initial R and phi for each timestep
    previouslyCalculatedDict = Dict()

    time = minimum(interval)
    endtime = maximum(interval)

    currentTimeStep = initialTimestep
    timestepRecorder = []
    attemptsRecorder = []

    S = nothing # This will always be updated later
    V = nothing # This will always be updated later
    ϕ = nothing # This will always be updated later
    modelFails = false # Do this if we reach a constraint violation always
    i = 1 # Enumerator


    # Strategy unique variables, remove for more effecient code
    # Strategy 2
    consecutiveFails = 0 # Keeps track of consecutive fails at initialTimeStep
    skipResetCounter = 0 # Counter for skipping resets. 0 means no skip, >0 means skip
    skipCompute = false  # Whether to skip computation, useful for when we have just been skipping

    # Strategy 3
    requiredSuccessCountForAttempt = 4 # Number of required consecutive successes to try double timestep

    # Strategy 4
    currentAttempt = -1 # Keep track of current attempt, -1 means no attempt
    initialChance = 0.5 # Initial chance of increasing in each profile
    backTrackDebth = 10 # Amount of previous steps to look at when choosing profile
    profilesAmount = 5 # Amount of profiles
    profileInterval = 3 # Amount of values between each profile

    profiles = []
    for i in 1:profilesAmount
        push!(profiles, initialChance)
    end

    # Function to change timestep, also updates changedTimeStep
    function changeTimeStep(newTimestep)
        currentTimeStep = newTimestep
        changedTimeStep = true
    end


    # Main loop
    while time < endtime
        if i% 100 == 0
            println("Time at step $i: $time")
        end

        if modelFails
            break
        end

        attempts = 1 # Keep track of number of attempts
        approveFlag = false
        newR = nothing # This will always be updated later
        newS = nothing # This will always be updated later
        newV = nothing # This will always be updated later
        newΩ = nothing # This will always be updated later

        # If exceeds endtime, change timestep
        # if time + currentTimeStep > endtime
        #     currentTimeStep = endtime - time
        #     changedTimeStep = true
        # end

        while !approveFlag
            if currentTimeStep == 0
                # We have a model fail
                println("Error model fails at time $time, constraint is not satisfied")
                modelFails = true
                break
            elseif changedTimeStep
                # We make sure to round the timestep
                currentTimeStep = round(currentTimeStep, digits=DIGITS)

                # We check if the new timestep respects the interval
                if mod(time, currentTimeStep) != 0
                    println("Cannot use $currentTimeStep as timestep at time $time")
                    currentTimeStep = currentTimeStep/2
                    continue
                end
            end

            if changedTimeStep
                # Check if we have already calculated the initial step
                if haskey(previouslyCalculatedDict, currentTimeStep)
                    if μ == 0
                        newR, ϕ = previouslyCalculatedDict[currentTimeStep]
                    else
                        newR, ballβ, ϕ = previouslyCalculatedDict[currentTimeStep]
                    end
                else
                    # We calculate and store it
                    if μ == 0
                        newR, ϕ = initialStepNoInput(A, ANorm, currentTimeStep, X0)
                        previouslyCalculatedDict[currentTimeStep] = (newR, ϕ)
                    else 
                        newR, ballβ, ϕ = initialStep(A, ANorm, currentTimeStep, X0, μ)
                        previouslyCalculatedDict[currentTimeStep] = (newR, ballβ, ϕ)
                    end
                end
                # Forward in time
                if μ == 0
                    newR = forwardTimeNoInput(A, newR, time)
                else
                    newR, newS, newV, newΩ = forwardTime(A, newR, time, currentTimeStep, ballβ, ϕ)
                end
                changedTimeStep = false
            else
                # Do the next step
                newR = linear_map(ϕ, R[i-1])
                if (μ != 0)
                    newS = minkowski_sum(S, V)
                    newV = linear_map(ϕ, V)
                    newΩ = minkowski_sum(newR, newS)
                end
            end
            
            # Hoping this gets compiled effeciently as μ is a constant
            checkArea = nothing
            if μ == 0
                checkArea = newR
            else
                checkArea = newΩ
            end

            # Check if we intersect with constraint
            if !intersects(constraint, checkArea)
                approveFlag = true
            
            else # Reduce timestep
                changeTimeStep(currentTimeStep/2)
                attempts = attempts + 1
            end

        end
        # Found a valid timestep

        println("Sucess at time $time with timestep $currentTimeStep")

        push!(R, newR)
        if (μ != 0)
            S = newS
            V = newV
            push!(Ω, newΩ)
        end

        push!(timestepRecorder, currentTimeStep)
        push!(attemptsRecorder, attempts)
        i = i + 1
        time = round(time + currentTimeStep, digits=DIGITS)

        # Reset / apply strategy
    
        if strategy == 1
            # Check if should reset
            if currentTimeStep != initialTimestep
                changeTimeStep(initialTimestep)
            end
        elseif strategy == 2
            # Apply logic to consecutive bad resets
            if skipResetCounter == 0
                if skipCompute # We have just been skipping, and wanna do 1 computation with reset
                    skipCompute = false
                else
                    # # If initial timestep is good, set fails to 0
                    # if currentTimeStep == initialTimestep
                    #     consecutiveFails = 0

                    # If reset found a better timestep
                    prevTimestep = initialTimestep
                    if i > 2 # Attempt to get previous timestep
                        prevTimestep = timestepRecorder[i-2]
                    end

                    if prevTimestep < currentTimeStep
                        consecutiveFails = 0
                    else # If we encounter a fail
                        consecutiveFails = consecutiveFails + 1
                        if consecutiveFails > 1
                            skipResetCounter = 2^(consecutiveFails-1) - 1
                            #println("Adding $skipResetCounter to skip counter at time $time")
                        end
                    end
                end
            end


            # Check if should skip reset, we check twice as it could be changed above
            if 0 < skipResetCounter
                skipResetCounter = skipResetCounter - 1
                skipCompute = true # Make sure we do 1 reset computation, before we can add more skips
                continue
            end
            # Apply reset
            if currentTimeStep != initialTimestep
                changeTimeStep(initialTimestep)
            end

            
        elseif strategy == 3
            flag = true
            # We look for any errors in past x steps
            for j in 1:min(i-1, requiredSuccessCountForAttempt)
                if attemptsRecorder[i-j] != 1
                    flag = false
                    break
                end
            end
            if flag #If we have x success in a row, double timestep
                changeTimeStep(currentTimeStep*2)
            end
        
        elseif strategy == 4
            # Check if we need to penalize previous(!)
            if currentAttempt != -1

                attemptAmount = attemptsRecorder[i-1]
                if attemptAmount == 1 # It was a success
                    profiles[currentAttempt] = min(0.9, profiles[currentAttempt] + 0.2)
                else # It was a fail
                    profiles[currentAttempt] = max(0.1, profiles[currentAttempt] - 0.1)
                end
                currentAttempt = -1
            end

            # Figure out which strategy to useful
            backTrack = min(i-1, backTrackDebth)
            accumulatedAttempts = -backTrackDebth # We reduce by the atleast 1 attempt needed
            for j in 1:backTrack
                accumulatedAttempts = accumulatedAttempts + attemptsRecorder[(i-j)] # -1 due to i being 1 ahead
            end 
            profileToUse = 1
            while accumulatedAttempts > profileInterval && profileToUse < profilesAmount
                accumulatedAttempts = accumulatedAttempts - profileInterval
                profileToUse = profileToUse + 1
            end

            # Use randomness to determine to increase or decrease
            # According to https://discourse.julialang.org/t/random-number-in-0-1/16009/4
            randomValue = rand()
            if randomValue === 0.0
                randomValue = 1
            end

            # We try double timestep maybe
            if randomValue < profiles[profileToUse]
                changeTimeStep(currentTimeStep*2)
                # We note down which profile we used, and that we made an attempt
                currentAttempt = profileToUse
            end




        else 
            println("Unknown strategy, ending program")
            modelFails = true
        end
    end
    # End of main loop

    if strategy == 4
        println("These profiles were used: ", profiles)
    end

    # Input or no input
    if μ == 0
        return (R, timestepRecorder, attemptsRecorder)
    end
    return (Ω, timestepRecorder, attemptsRecorder)
end


function reachsetslog(A, timestepsize, interval, X₀, μ)
    T = maximum(interval)
    N = floor(Int, T/timestepsize)
    
    bits = ceil(Int, log2(N))
    println(N)
    println(bits)
    ϕs = Vector{AbstractMatrix}(undef, bits)

    ballβs = Vector{Zonotope}(undef, bits)



    ANorm = norm(A, Inf)
    α = (exp(ANorm*timestepsize)-1-timestepsize*ANorm)/norm(X₀, Inf)
    β = (exp(ANorm*timestepsize)-1)*μ/ANorm

    ϕ = exp(A*timestepsize)

    ϕp = (I+ϕ)/2
    ϕm = (I-ϕ)/2
    gens = hcat(ϕp*X₀.generators,ϕm*X₀.center, ϕm*X₀.generators)

    R₁ = minkowski_sum(Zonotope(ϕp*X₀.center, gens), Zonotope(zeros(2), (α+β)*I(2)))
    #res = [R₁]
    R = [R₁]
    ts = [1]

    ballβ = Zonotope(zeros(2), β*I(2))

    ϕs[1] = ϕ
    ballβs[1] = ballβ

    for i in 2:bits
        ϕs[i] = ϕs[i-1]*ϕs[i-1]
        ballβs[i] = minkowski_sum(ballβs[i-1], ballβs[i-1])
    end

    for i in 2:bits
        if ((N >> (i)) & 1) == 1
            push!(R, minkowski_sum(linear_map(ϕs[i], R[end]), ballβs[i]))
            push!(ts, ts[end]+(2^(i)))
        end
            #res[i] = R[i]
    end
    return R, ts
end

function approximateArea(zonotope::Zonotope)
    zonotope |> box_approximation |> area
end

function approximateArea(zonotope::AbstractVector{Zonotope{Float64, Vector{Float64}, Matrix{Float64}}})
    area(reduce(∪, map(box_approximation, zonotope)))
end

#intersects function is based upon "Set operations and order reductions for constrained zonotopes" - Vignesh Raghuraman, Justin P. Koeln
#Note that this only returns true if they intersect and does not depend on whether the intersection between them is empty.
function intersects(zonotope::Zonotope, halfspace::HalfSpace)
    h, f = tosimplehrep(halfspace)
    l = abs(only(f) - dot(h', zonotope.center))
    r = sum(abs.(zonotope.generators .* h'))
    return l <= r
end

function intersects(halfspace::HalfSpace, zonotope::Zonotope)
    intersects(zonotope, halfspace)
end

function intersects(z1::Zonotope, z2::Zonotope) #Should be possible to compute this from above paper by viewing the zonotopes as being constrained but with no constrains and whether their intersecting zonotope has empty 
    isdisjoint(z1, z2)
end