
levelStartFile = "levelStart.state"

populationSize = 200
maxTimeout = 500
maxNeurons = 10000


--new species constructor
function newSpecies()
    local species = {}
    species.genomes = {}
    species.topFitness = 0
    species.averageFitness = 0
    species.staleness = 0
    return species
end

--new genome constructor
function newGenome()
	local genome = {}
	genome.genes = {}
	genome.fitness = 0
	genome.adjustedFitness = 0
	genome.network = {}
	genome.maxneuron = 0
	genome.globalRank = 0
	genome.mutationRates = {}
	genome.mutationRates["connections"] = MutateConnectionsChance
	genome.mutationRates["link"] = LinkMutationChance
	genome.mutationRates["bias"] = BiasMutationChance
	genome.mutationRates["node"] = NodeMutationChance
	genome.mutationRates["enable"] = EnableMutationChance
	genome.mutationRates["disable"] = DisableMutationChance
	genome.mutationRates["step"] = StepSize
	
	return genome
end

--create a new genome with slight mutation
function newSimpleGenome()
    local genome = newGenome()
    genome.maxneuron = numInputs
    mutate(genome)
    return genome
end

--new gene constructor
function newGene()
    local gene = {}
    gene.into = 0
    gene.out = 0
    gene.w = 0.0
    gene.enabled = true
    gene.innovation = 0
    return gene
end

--new neuron constructor
function newNeuron()
   local neuron = {}
   neuron.genesInto = {}
   neuron.value = 0.0
   return neuron
end

--new session constructor
function newSession()
    local session = {}
    --current species, etc storage
    session.curGeneration = 0
    session.species = {}
    session.curSpeciesNum = 1
    session.curGenomeNum = 1
    session.innovationNum = numOutputs

    --fitness evals
    session.curFrame = 0
    session.highestPoint = 0
    session.maxFitness = 0

    return session
end

--checks if a genome belongs in the same species as another genome
function fitsSpecies(newGenome, speciesGenome)
    local a = DeltaDisjoint*disjoint(newGenome.genes, speciesGenome.genes)
    local b = DeltaWeights*weights(newGenome.genes, speciesGenome.genes) 
    return a + b < DeltaThreshold
end

--group the genome to an existing or new species
function groupGenome(genomeToGroup)

    --try to find a species the genome belongs in
    for _, species in pairs(session.species) do
        if fitsSpecies(genomeToGroup, species.genomes[1]) then
            table.insert(species.genomes, genomeToGroup) 
            return nil
        end
    end

    --create a new species for the genome
    local species = newSpecies()
    table.insert(species.genomes, genomeToGroup)
    table.insert(session.species, species)

end

--create the network for the specified genome
function buildGenomeNetwork(genome)
    local network = {}
    network.neurons = {}

    --create input neurons
    for i=1,numInputs do
        network.neurons[i] = newNeuron()
    end
    --create output neurons
    for i=1,numOutputs do
        network.neurons[i + maxNeurons] = newNeuron()
    end
    --sort genes in the genome
    table.sort(genome.genes, function(a,b) return a.out < b.out end)

    --create the hidden layer neurons and connections
    for _, gene in pairs(genome.genes) do
        if gene.enabled then
            --if the neurons that the gene connects to are not yet created, create them
            if network.neurons[gene.into] == nil then
                network.neurons[gene.into] = newNeuron()
            end
            if network.neurons[gene.out] == nil then
                network.neurons[gene.out] = newNeuron()
            end
            --insert the gene as an incoming edge to the out neuron of that gene
            table.insert(network.neurons[gene.out].genesInto, gene)
        end
    end
    genome.network = network
end

--start a new learning session
function startSession()

    session = newSession()

    --create a the new population and group them
    for p=1,populationSize do
        genome = newSimpleGenome()
        groupGenome(genome)  
    end

    --start running the first genome
    runCurrentGenome()

end

--set controller to press no buttons
function resetController()
    controller = {}
    for _, button in pairs(controllerButtons) do
        controller["P1"..button] = false
    end
    joypad.set(controller)
end

--get the inputs and apply them to the controller with a given genome, return outputs
function evaluateGenome(genome)
    inputs = getInputs()
    local controller = evaluateNetwork(genome.network, inputs)
    joypad.set(controller)
    return controller
end

--evaluates a genomes network with the game inputs to determine buttons to press
function evaluateNetwork(network, inputs)
    --insert the bias neuron
    table.insert(inputs, 1)

    --set input neurons to inputs
    for i=1,numInputs do
        netowrk.neurons[i].value = inputs[i]
    end

    --evaluate hidden layer neurons
    for _, neuron in pairs(network.neurons) do
        local sum = 0.0
        for _,geneIn in pairs(neuron.genesInto) do
            --find the value of the gene coming into the neuron
            local geneValue = geneIn.weight * network.neurons[geneIn.into]
            sum = sum + geneValue
        end
        if #neuron.genesInto > 0 then
            neuron.value = sigmoid(sum)
        end
    end

    --evaluate outputs
    local outputs = {}

    --eval jump seperately
    local jumpButton = "P1 " .. controllerButtons[1]
    if network.neurons[maxNeurons + 1].value > 0 then
            outputs[jumpButton] = true
    else
            outputs[button] = false
    end
    
    
    --set the maximum direction output to true since no direction button can be pressed at the same time
    --on a real controller
    local maxIndex = 2
    local maxValue = -100000
    for i=2,numOutputs do
            local button = "P1 " .. controllerButtons[i]
            outputs[button] = false
            if network.neurons[maxNeurons+i].value >= maxValue then
                    maxValue = network.neurons[maxNeurons+i].value
                    maxIndex = i
            end
    end
    local button = "P1 " .. controllerButtons[maxIndex]
    outputs[button] = true
    
    return outputs
end

--get the next genome that has not been evaluated, if none exist, evolve a new set of species
function findNextGenome()
    local nextSpecies = 1
    local nextGenome = 1

    local nextGenome = session.species[nextSpecies].genomes[nextGenome]
    while nextGenome.fitness ~= 0 do
        nextGenome = nextGenome + 1
        if nextGenome > #session.species[nextSpecies].genomes then
            nextGenome = 1
            nextSpecies = nextSpecies + 1
            if nextSpecies > #session.species then
                evolveNewGeneration()
                nextSpecies = 1
            end
        end
    end

    return {["s"]=nextSpecies, ["g"]=nextGenome}

end

--reload the save state point and reinitialize the network being used
function runCurrentGenome()
    savestate.load(levelStartFile) 
    session.highestPoint = 0
    session.currentFrame = 0
    session.timeout = maxTimeout
    resetController()

    --get the next genome
    local nextGenome = session.species[session.curSpeciesNum].genomes[session.curGenomeNum]
    buildGenomeNetwork(nextGenome)
    evaluateGenome(nextGenome)

end

function getFitness(marioLoc)
    local isJumping = memory.readbyte(0x96) == 0x04
    local isOnBadLadder = marioX == 84 or marioX == 100 or marioX == 108 or marioX == 188
    --valid fitness increase situation
    if not isJumping and not isDead() and not isOnBadLadder then
        --new highest point, reset timout
        local normalizedHeight = 255-marioLoc['y']-51
        if normalizedHeight > session.highestPoint then
            session.highestPoint = normalizedHeight
            session.timeout = maxTimeout
        end
    end
    --punish for going to higher points that are not on the ideal path
    local punishment = 0
    if marioLoc['x'] <= 35 or marioLoc['x'] >= 215 then
        punishment = 20
    end
    local fitness = session.highestPoint*10 - session.currentFrame / 150.0 - punishment 
    if fitness < 0 then
        fitness = -1
    end
    return fitness
end

function getPercentEvaluated()
    local total = 0.0
    local evaluated = 0.0
    for _,species in pairs(session.species) do
        for _,genome in pairs(species.genomes) do
            total = total + 1
            if genome.fitness ~= 0 then
                evaluated = evaluated + 1
            end
        end
    end
    return math.floor(evaluated * 100.0 / total)
end


--initalize a new session
startSession()
controller = {}

--main game loop
while true do
    local curGenome = session.species[session.curSpeciesNum].genomes[session.curGenomeNum]
    
    --set the controller according to the current genome 
    if session.currentFrame%5 == 0 then
        controller = evaluateGenome(curGenome)
    end
    joypad.set(controller)

    local fitness = getFitness()
    local percDone = getPercentEvaluated()

    --check for game over or timeout
    session.timeout = session.timeout - 1
    local livingTimeBonus = session.currentFrame / 4
    if session.timeout + livingTimeBonus <= 0 or isDead() then

        --set final fitness value
        curGenome.fitness = fitness
        if fitness > session.maxFitness then
            session.maxFitness = fitness
        end

        --get the next genome to evaluate
        local nextGenome = findNextGenome()
        session.curGenomeNum = nextGenome["g"]
        session.curSpeciesNum = nextGenome["s"]
        --run the next genome
        runCurrentGenome()

    end

    --step up the session frame count and advance the emulator
    session.currentFrame = session.currentFrame + 1
    emu.frameadvance()

end
