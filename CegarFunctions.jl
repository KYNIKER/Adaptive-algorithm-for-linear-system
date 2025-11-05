using LazySets, LinearAlgebra

using Random # For strategy 4
include("helperfunctions.jl")
"""
    reachSetsCegar(A, initialTimeStep, interval, X0, constraint, strategy = 1, DIGITS = 4)
# CEGAR like approach

Different strategies for how to choose the timestep
- 1: We half the timestpe on fail, and reset each time
- 2: We half the timestep on fail. Upon consecutive fails we take 2^(x-1) - 1 steps before attempting reset, where x is the number of consecutive fails. We reset to initial timestep
- 3: Never reset. If no errors in past x steps, double timestep
- 4: Randomly choose to retain or increase

Adjust chance of increase based on previous success of increase. 
Keep different profiles dependt on amount of recent fails. 
Specifically 0-1, 2-3, 4-5, 6+
"""
function reachSetsCegar(A, initialTimestep, interval, X0, constraint::HalfSpace, strategy::Int = 1, DIGITS = 4)
    if strategy == 0
        return strategy0(A, initialTimestep, interval, X0, constraint, DIGITS)
    elseif strategy == 1
        return strategy1(A, initialTimestep, interval, X0, constraint, DIGITS)
    elseif strategy == 2
        return strategy2(A, initialTimestep, interval, X0, constraint, DIGITS)
    elseif strategy == 3
        return strategy3(A, initialTimestep, interval, X0, constraint, DIGITS)
    elseif strategy == 4
        return strategy4(A, initialTimestep, interval, X0, constraint, DIGITS)
    else 
        throw(DomainError(strategy, "argument must be 0, 1, 2, 3 or 4"))
    end
end

# Naive solution. Half timestep on fail and reset completely
function strategy0(A, initialTimestep, interval, X0, constraint, DIGITS = 4)
    startTime = minimum(interval)
    endTime = maximum(interval)
    finishFlag = false
    timestep = initialTimestep
    ANorm = norm(A, Inf)
    i = 1
    # Run till finished, or fail
    R = []
    while !finishFlag
        intersectFlag = false
        newR, ϕ = initialStepNoInput(A, ANorm, timestep, X0)
        # Update R to fit startTime
        push!(R, forwardTimeNoInput(A, newR, startTime))

        currentTime = startTime + timestep
        while currentTime < endTime && intersectFlag == false
            i = i + 1
            newR = linear_map(ϕ, R[i-1])

            if intersects(constraint, newR)
                intersectFlag = true
            else
                push!(R, newR)
                currentTime = currentTime + timestep
            end
        end
        if intersectFlag
            timestep = timestep / 2
            R = [] # Reset timer
            i = 1
        else
            finishFlag = true
        end
    end

    # Set timestep recorder
    amountOfPoints = length(R)

    timestepRecorder = collect(Iterators.repeated(timestep, amountOfPoints))
    attemptsRecorder = collect(Iterators.repeated(1, amountOfPoints))

    # End of main loop
    return (R, timestepRecorder, attemptsRecorder)
end

function strategy1(A, initialTimestep, interval, X0, constraint, DIGITS = 4)
    ANorm = norm(A, Inf)
    R = []
    changedTimeStep = true # Keeps track of whether timestep has changed

    # Dictionary to keep track of initial R and phi for each timestep
    previouslyCalculatedDict = Dict()

    time = minimum(interval)
    endtime = maximum(interval)

    currentTimeStep = initialTimestep
    timestepRecorder = []
    attemptsRecorder = []

    ϕ = nothing # This will always be updated later
    i = 1 # Enumerator

    # Main loop
    while time < endtime
        if i% 100 == 0
            println("Time at step $i: $time")
        end

        attempts = 1 # Keep track of number of attempts
        approveFlag = false
        newR = nothing # This will always be updated later

        currentTimeStep, newR, ϕ, changedTimeStep, attempts = fitTimeStep(currentTimeStep, changedTimeStep, previouslyCalculatedDict, A, ANorm, X0, R, ϕ, attempts, time, i, constraint, DIGITS)

        # Found a valid timestep
        push!(R, newR)

        push!(timestepRecorder, currentTimeStep)
        push!(attemptsRecorder, attempts)
        i = i + 1
        time = round(time + currentTimeStep, digits=DIGITS)

        # Reset / apply strategy
        # Check if should reset
        if currentTimeStep != initialTimestep
            currentTimeStep = initialTimestep
            changedTimeStep = true
        end
    end
    # End of main loop
    return (R, timestepRecorder, attemptsRecorder)
end

function strategy2(A, initialTimestep, interval, X0, constraint, DIGITS = 4)
    ANorm = norm(A, Inf)
    R = []
    changedTimeStep = true # Keeps track of whether timestep has changed

    # Dictionary to keep track of initial R and phi for each timestep
    previouslyCalculatedDict = Dict()

    time = minimum(interval)
    endtime = maximum(interval)

    currentTimeStep = initialTimestep
    timestepRecorder = []
    attemptsRecorder = []

    ϕ = nothing # This will always be updated later
    i = 1 # Enumerator

    consecutiveFails = 0 # Keeps track of consecutive fails at initialTimeStep
    skipResetCounter = 0 # Counter for skipping resets. 0 means no skip, >0 means skip
    skipCompute = false  # Whether to skip computation, useful for when we have just been skipping

    # Main loop
    while time < endtime
        if i% 100 == 0
            println("Time at step $i: $time")
        end

        attempts = 1 # Keep track of number of attempts
        newR = nothing # This will always be updated later

        currentTimeStep, newR, ϕ, changedTimeStep, attempts = fitTimeStep(currentTimeStep, changedTimeStep, previouslyCalculatedDict, A, ANorm, X0, R, ϕ, attempts, time, i, constraint, DIGITS)

        # Found a valid timestep
        push!(R, newR)

        push!(timestepRecorder, currentTimeStep)
        push!(attemptsRecorder, attempts)
        i = i + 1
        time = round(time + currentTimeStep, digits=DIGITS)

        # Reset / apply strategy
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
            currentTimeStep = initialTimestep
            changedTimeStep = true
        end
    end
    # End of main loop
    return (R, timestepRecorder, attemptsRecorder)
end

function strategy3(A, initialTimestep, interval, X0, constraint, DIGITS = 4)
    ANorm = norm(A, Inf)
    R = []
    changedTimeStep = true # Keeps track of whether timestep has changed

    # Dictionary to keep track of initial R and phi for each timestep
    previouslyCalculatedDict = Dict()

    time = minimum(interval)
    endtime = maximum(interval)

    currentTimeStep = initialTimestep
    timestepRecorder = []
    attemptsRecorder = []

    ϕ = nothing # This will always be updated later
    i = 1 # Enumerator

    requiredSuccessCountForAttempt = 4 # Number of required consecutive successes to try double timestep

    # Main loop
    while time < endtime
        if i% 100 == 0
            println("Time at step $i: $time")
        end

        attempts = 1 # Keep track of number of attempts
        newR = nothing # This will always be updated later

        currentTimeStep, newR, ϕ, changedTimeStep, attempts = fitTimeStep(currentTimeStep, changedTimeStep, previouslyCalculatedDict, A, ANorm, X0, R, ϕ, attempts, time, i, constraint, DIGITS)

        # Found a valid timestep
        push!(R, newR)

        push!(timestepRecorder, currentTimeStep)
        push!(attemptsRecorder, attempts)
        i = i + 1
        time = round(time + currentTimeStep, digits=DIGITS)

        # Reset / apply strategy
        flag = true
        # We look for any errors in past x steps
        for j in 1:min(i-1, requiredSuccessCountForAttempt)
            if attemptsRecorder[i-j] != 1
                flag = false
                break
            end
        end
        if flag #If we have x success in a row, double timestep
            currentTimeStep = currentTimeStep*2
            changedTimeStep = true
        end
    end
    # End of main loop
    return (R, timestepRecorder, attemptsRecorder)
end

function strategy4(A, initialTimestep, interval, X0, constraint, DIGITS = 4)
    ANorm = norm(A, Inf)
    R = []
    changedTimeStep = true # Keeps track of whether timestep has changed

    # Dictionary to keep track of initial R and phi for each timestep
    previouslyCalculatedDict = Dict()

    time = minimum(interval)
    endtime = maximum(interval)

    currentTimeStep = initialTimestep
    timestepRecorder = []
    attemptsRecorder = []

    ϕ = nothing # This will always be updated later
    i = 1 # Enumerator

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
    
    # Main loop
    while time < endtime
        if i% 100 == 0
            println("Time at step $i: $time")
        end

        attempts = 1 # Keep track of number of attempts
        newR = nothing # This will always be updated later

        currentTimeStep, newR, ϕ, changedTimeStep, attempts = fitTimeStep(currentTimeStep, changedTimeStep, previouslyCalculatedDict, A, ANorm, X0, R, ϕ, attempts, time, i, constraint, DIGITS)

        # Found a valid timestep
        push!(R, newR)

        push!(timestepRecorder, currentTimeStep)
        push!(attemptsRecorder, attempts)
        i = i + 1
        time = round(time + currentTimeStep, digits=DIGITS)

        # Reset / apply strategy
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
            currentTimeStep = currentTimeStep * 2
            changedTimeStep = true
            # We note down which profile we used, and that we made an attempt
            currentAttempt = profileToUse
        end
    end
    # End of main loop
    return (R, timestepRecorder, attemptsRecorder)
end


function reachSetsCegarInput(A, initialTimestep, interval, X0, constraint::HalfSpace, μ::Float64, strategy::Int = 1, DIGITS = 4, REUSE = false)
    if strategy == 1
        return strategy1Input(A, initialTimestep, interval, X0, constraint, μ, DIGITS, REUSE)
    elseif strategy == 2
        return strategy2Input(A, initialTimestep, interval, X0, constraint, μ, DIGITS, REUSE)
    elseif strategy == 3
        return strategy3Input(A, initialTimestep, interval, X0, constraint, μ, DIGITS, REUSE)
    elseif strategy == 4
        return strategy4Input(A, initialTimestep, interval, X0, constraint, μ, DIGITS, REUSE)
    else 
        throw(DomainError(strategy, "argument must be 1, 2, 3 or 4"))
    end
end

function strategy1Input(A, initialTimestep, interval, X0, constraint, μ, DIGITS::Int, REUSE::Bool)
    ANorm = norm(A, Inf)
    R = []
    Ω = []
    changedTimeStep = true # Keeps track of whether timestep has changed

    # Dictionary to keep track of initial R and phi for each timestep
    previouslyCalculatedDict = Dict()
    
    previousInputValuesDict = Dict()

    time = minimum(interval)
    endtime = maximum(interval)

    currentTimeStep = initialTimestep
    timestepRecorder = []
    attemptsRecorder = []

    S = REUSE ? Zonotope(zeros(dim(X0)), [zeros(dim(X0))]) : nothing # This will always be updated later
    V = nothing # This will always be updated later
    ϕ = nothing # This will always be updated later
    modelFails = false # Do this if we reach a constraint violation always
    i = 1 # Enumerator

    # Main loop
    while time < endtime
        if i% 100 == 0
            println("Time at step $i: $time")
        end

        if modelFails
            break
        end

        attempts = 1 # Keep track of number of attempts
        newR = nothing # This will always be updated later
        newΩ = nothing # This will always be updated later
        prevTime = 0

        currentTimeStep,  newR, newΩ, ϕ, changedTimeStep, attempts, newS, newV = REUSE ?
            fitTimeStepInputReUse(currentTimeStep, changedTimeStep, previouslyCalculatedDict, previousInputValuesDict, A, ANorm, X0, R, ϕ, attempts, time, i,μ, S, V, constraint, DIGITS) :
            fitTimeStepInput(currentTimeStep, changedTimeStep, previouslyCalculatedDict, previousInputValuesDict, A, ANorm, X0, R, ϕ, attempts, time, i,μ, S, V, constraint, DIGITS)

        # Found a valid timestep
        push!(R, newR)
    
        S = newS
        V = newV
        push!(Ω, newΩ)


        push!(timestepRecorder, currentTimeStep)
        push!(attemptsRecorder, attempts)
        i = i + 1
        time = round(time + currentTimeStep, digits=DIGITS)

        # Reset / apply strategy
        if currentTimeStep != initialTimestep
            currentTimeStep = initialTimeStep
            changedTimeStep = true
        end
        
    end
    # End of main loop
    return (Ω, timestepRecorder, attemptsRecorder)
end

function strategy2Input(A, initialTimestep, interval, X0, constraint, μ, DIGITS::Int, REUSE::Bool)
    ANorm = norm(A, Inf)
    R = []
    Ω = []
    changedTimeStep = true # Keeps track of whether timestep has changed

    # Dictionary to keep track of initial R and phi for each timestep
    previouslyCalculatedDict = Dict()
    
    previousInputValuesDict = Dict()

    time = minimum(interval)
    endtime = maximum(interval)

    currentTimeStep = initialTimestep
    timestepRecorder = []
    attemptsRecorder = []

    S = REUSE ? Zonotope(zeros(dim(X0)), [zeros(dim(X0))]) : nothing # This will always be updated later
    V = nothing # This will always be updated later
    ϕ = nothing # This will always be updated later
    modelFails = false # Do this if we reach a constraint violation always
    i = 1 # Enumerator

    consecutiveFails = 0 # Keeps track of consecutive fails at initialTimeStep
    skipResetCounter = 0 # Counter for skipping resets. 0 means no skip, >0 means skip
    skipCompute = false  # Whether to skip computation, useful for when we have just been skipping


    # Main loop
    while time < endtime
        if i% 100 == 0
            println("Time at step $i: $time")
        end

        if modelFails
            break
        end

        attempts = 1 # Keep track of number of attempts
        newR = nothing # This will always be updated later
        newΩ = nothing # This will always be updated later
        prevTime = 0

        currentTimeStep,  newR, newΩ, ϕ, changedTimeStep, attempts, newS, newV = REUSE ?
            fitTimeStepInputReUse(currentTimeStep, changedTimeStep, previouslyCalculatedDict, previousInputValuesDict, A, ANorm, X0, R, ϕ, attempts, time, i,μ, S, V, constraint, DIGITS) : 
            fitTimeStepInput(currentTimeStep, changedTimeStep, previouslyCalculatedDict, previousInputValuesDict, A, ANorm, X0, R, ϕ, attempts, time, i,μ, S, V, constraint, DIGITS)

        # Found a valid timestep
        push!(R, newR)
    
        S = newS
        V = newV
        push!(Ω, newΩ)


        push!(timestepRecorder, currentTimeStep)
        push!(attemptsRecorder, attempts)
        i = i + 1
        time = round(time + currentTimeStep, digits=DIGITS)

        # Reset / apply strategy
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
            currentTimeStep = initialTimeStep
            changedTimeStep = true
        end
        
    end
    # End of main loop
    return (Ω, timestepRecorder, attemptsRecorder)
end

function strategy3Input(A, initialTimestep, interval, X0, constraint, μ, DIGITS::Int, REUSE::Bool)
    ANorm = norm(A, Inf)
    R = []
    Ω = []
    changedTimeStep = true # Keeps track of whether timestep has changed

    # Dictionary to keep track of initial R and phi for each timestep
    previouslyCalculatedDict = Dict()
    
    previousInputValuesDict = Dict()

    time = minimum(interval)
    endtime = maximum(interval)

    currentTimeStep = initialTimestep
    timestepRecorder = []
    attemptsRecorder = []

    S = REUSE ? Zonotope(zeros(dim(X0)), [zeros(dim(X0))]) : nothing # This will always be updated later
    V = nothing # This will always be updated later
    ϕ = nothing # This will always be updated later
    modelFails = false # Do this if we reach a constraint violation always
    i = 1 # Enumerator

    requiredSuccessCountForAttempt = 4 # Number of required consecutive successes to try double timestep

    # Main loop
    while time < endtime
        if i% 100 == 0
            println("Time at step $i: $time")
        end

        if modelFails
            break
        end

        attempts = 1 # Keep track of number of attempts
        newR = nothing # This will always be updated later
        newΩ = nothing # This will always be updated later

        currentTimeStep, newR, newΩ, ϕ, changedTimeStep, attempts, newS, newV = REUSE ? 
            fitTimeStepInputReUse(currentTimeStep, changedTimeStep, previouslyCalculatedDict, previousInputValuesDict, A, ANorm, X0, R, ϕ, attempts, time, i,μ, S, V, constraint, DIGITS) : 
            fitTimeStepInput(currentTimeStep, changedTimeStep, previouslyCalculatedDict, previousInputValuesDict, A, ANorm, X0, R, ϕ, attempts, time, i,μ, S, V, constraint, DIGITS)

        # Found a valid timestep
        push!(R, newR)
    
        S = newS
        V = newV
        push!(Ω, newΩ)


        push!(timestepRecorder, currentTimeStep)
        push!(attemptsRecorder, attempts)
        i = i + 1
        time = round(time + currentTimeStep, digits=DIGITS)

        # Reset / apply strategy
        flag = true
        # We look for any errors in past x steps
        for j in 1:min(i-1, requiredSuccessCountForAttempt)
            if attemptsRecorder[i-j] != 1
                flag = false
                break
            end
        end
        if flag #If we have x success in a row, double timestep
            currentTimeStep = currentTimeStep * 2
            changedTimeStep = true
        end
        
    end
    # End of main loop
    return (Ω, timestepRecorder, attemptsRecorder)
end

function strategy4Input(A, initialTimestep, interval, X0, constraint, μ, DIGITS::Int, REUSE::Bool)
    ANorm = norm(A, Inf)
    R = []
    Ω = []
    changedTimeStep = true # Keeps track of whether timestep has changed

    # Dictionary to keep track of initial R and phi for each timestep
    previouslyCalculatedDict = Dict()
    
    previousInputValuesDict = Dict()

    time = minimum(interval)
    endtime = maximum(interval)

    currentTimeStep = initialTimestep
    timestepRecorder = []
    attemptsRecorder = []

    S = REUSE ? Zonotope(zeros(dim(X0)), [zeros(dim(X0))]) : nothing # This will always be updated later
    V = nothing # This will always be updated later
    ϕ = nothing # This will always be updated later
    modelFails = false # Do this if we reach a constraint violation always
    i = 1 # Enumerator

    currentAttempt = -1 # Keep track of current attempt, -1 means no attempt
    initialChance = 0.5 # Initial chance of increasing in each profile
    backTrackDebth = 10 # Amount of previous steps to look at when choosing profile
    profilesAmount = 5 # Amount of profiles
    profileInterval = 3 # Amount of values between each profile

    profiles = []
    for i in 1:profilesAmount
        push!(profiles, initialChance)
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
        newR = nothing # This will always be updated later
        newΩ = nothing # This will always be updated later
        prevTime = 0

        currentTimeStep,  newR, newΩ, ϕ, changedTimeStep, attempts, newS, newV = REUSE ? 
            fitTimeStepInputReUse(currentTimeStep, changedTimeStep, previouslyCalculatedDict, previousInputValuesDict, A, ANorm, X0, R, ϕ, attempts, time, i,μ, S, V, constraint, DIGITS) :
            fitTimeStepInput(currentTimeStep, changedTimeStep, previouslyCalculatedDict, previousInputValuesDict, A, ANorm, X0, R, ϕ, attempts, time, i,μ, S, V, constraint, DIGITS)

        # Found a valid timestep
        push!(R, newR)
    
        S = newS
        V = newV
        push!(Ω, newΩ)


        push!(timestepRecorder, currentTimeStep)
        push!(attemptsRecorder, attempts)
        i = i + 1
        time = round(time + currentTimeStep, digits=DIGITS)

        # Reset / apply strategy
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
            currentTimeStep = currentTimeStep * 2
            changedTimeStep = true
            # We note down which profile we used, and that we made an attempt
            currentAttempt = profileToUse
        end
        
    end
    # End of main loop
    return (Ω, timestepRecorder, attemptsRecorder)
end

function fitTimeStep(currentTimeStep, changedTimeStep, previouslyCalculatedDict, A, ANorm, X0, R, ϕ, attempts, time, i, constraint, DIGITS::Int)
    approveFlag = false
    newR = nothing
    while !approveFlag
        if changedTimeStep
            # We make sure to round the timestep
            currentTimeStep = round(currentTimeStep, digits=DIGITS)

            if currentTimeStep == 0
                # We have a model fail
                throw(AssertionError("Error model fails at time $time, constraint is not satisfied"))
            end

            # Check if we have already calculated the initial step, else calculate it
            if !haskey(previouslyCalculatedDict, currentTimeStep)
                previouslyCalculatedDict[currentTimeStep] = initialStepNoInput(A, ANorm, currentTimeStep, X0)
            end
            newR, ϕ = previouslyCalculatedDict[currentTimeStep]

            # Forward in time
            newR = forwardTimeNoInput(A, newR, time)
            changedTimeStep = false
        else
            # Do the next step
            newR = linear_map(ϕ, R[i-1])
        end

        # Check if we intersect with constraint
        if !intersects(constraint, newR)
            approveFlag = true
        else # Reduce timestep
            currentTimeStep = currentTimeStep / 2
            changedTimeStep = true
            attempts = attempts + 1
        end
    end
    return currentTimeStep, newR, ϕ, changedTimeStep, attempts
end

function fitTimeStepInput(currentTimeStep, changedTimeStep, previouslyCalculatedDict,previousInputValuesDict, A, ANorm, X0, R, ϕ, attempts, time, i,μ, S, V, constraint, DIGITS::Int)
    approveFlag = false
    newR = nothing # This will always be updated later
    newS = nothing # This will always be updated later
    newV = nothing # This will always be updated later
    newΩ = nothing # This will always be updated later
    prevTime = 0
    while !approveFlag
        if changedTimeStep
            # We make sure to round the timestep
            currentTimeStep = round(currentTimeStep, digits=DIGITS)

            if currentTimeStep == 0
                # We have a model fail
                throw(AssertionError("Error model fails at time $time, constraint is not satisfied"))
            end

            # We check if the new timestep respects the interval. Only do this if we have an input
            if mod(rationalize(time), rationalize(currentTimeStep)) != 0
                #println("Cannot use $currentTimeStep as timestep at time $time")
                currentTimeStep = currentTimeStep/2
                continue
            end
        end

        if changedTimeStep
            # Check if we have already calculated the initial step
            if haskey(previouslyCalculatedDict, currentTimeStep)
                newR, ϕ = previouslyCalculatedDict[currentTimeStep]
                S, V, prevTime = previousInputValuesDict[currentTimeStep]
            else
                # We calculate and store it
                newR, ballβ, ϕ = initialStep(A, ANorm, currentTimeStep, X0, μ)
                previouslyCalculatedDict[currentTimeStep] = (newR, ϕ)
                previousInputValuesDict[currentTimeStep] = (S, V, time)
                # These should initially be ballβ
                S = ballβ
                V = ballβ
            end
            # Forward in time
            newR, newS, newV, newΩ = forwardTime(A, newR, time, currentTimeStep, ϕ, prevTime, S, V)
            previousInputValuesDict[currentTimeStep] = (newS, newV, time)
        
            changedTimeStep = false
        else
            # Do the next step
            newR = linear_map(ϕ, R[i-1])
            newS = minkowski_sum(S, V)
            newV = linear_map(ϕ, V)
            newΩ = minkowski_sum(newR, newS)
            previousInputValuesDict[currentTimeStep] = (newS, newV, time)

        end

        # Check if we intersect with constraint
        if !intersects(constraint, newΩ)
            approveFlag = true
        else # Reduce timestep
            currentTimeStep = (currentTimeStep/2)
            changedTimeStep = true
            attempts = attempts + 1
        end
    end
    return currentTimeStep, newR, newΩ, ϕ, changedTimeStep, attempts, newS, newV
end

function fitTimeStepInputReUse(currentTimeStep, changedTimeStep, previouslyCalculatedDict,previousInputValuesDict, A, ANorm, X0, R, ϕ, attempts, time, i,μ, S, V, constraint, DIGITS::Int)
    approveFlag = false
    newR = nothing # This will always be updated later
    newS = nothing # This will always be updated later
    newV = nothing # This will always be updated later
    newΩ = nothing # This will always be updated later
    ballβ = nothing
    while !approveFlag
        if changedTimeStep
            # We make sure to round the timestep
            currentTimeStep = round(currentTimeStep, digits=DIGITS)

            if currentTimeStep == 0
                # We have a model fail
                throw(AssertionError("Error model fails at time $time, constraint is not satisfied"))
            end
        end

        if changedTimeStep
            # Check if we have already calculated the initial step
            if haskey(previouslyCalculatedDict, currentTimeStep)
                newR, ϕ = previouslyCalculatedDict[currentTimeStep]
                ballβ = previousInputValuesDict[currentTimeStep]
            else
                # We calculate and store it
                newR, ballβ, ϕ = initialStep(A, ANorm, currentTimeStep, X0, μ)
                previouslyCalculatedDict[currentTimeStep] = (newR, ϕ)
                previousInputValuesDict[currentTimeStep] = ballβ
            end
            # Forward in time
            newR, newS, newV, newΩ = forwardTimeReUse(A, newR, time, ϕ, ballβ, S)
            changedTimeStep = false
        else
            # Do the next step
            newR = linear_map(ϕ, R[i-1])
            newS = minkowski_sum(S, V)
            newV = linear_map(ϕ, V)
            newΩ = minkowski_sum(newR, newS)
        end

        # Check if we intersect with constraint
        if !intersects(constraint, newΩ)
            approveFlag = true
        else # Reduce timestep
            currentTimeStep = (currentTimeStep/2)
            changedTimeStep = true
            attempts = attempts + 1
        end
    end
    return currentTimeStep, newR, newΩ, ϕ, changedTimeStep, attempts, newS, newV
end