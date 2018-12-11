
levelStartFile = "levelStart.state"

populationSize = 200
maxTimeout = 500
maxNeurons = 10000

--new species constructor
function newSpecies()
    local species = {}
    species.genomes = {}
    species.topFitness = 0
    species.averageRank = 0
    species.numGensStuck = 0
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
    genome.rank = 0
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

function duplicateGenome(dup)
    local genome = newGenome()
    for _,gene in pairs(dup.genes) do
            table.insert(genome.genes, copyGene(gene))
    end
    genome.maxneuron = dup.maxneuron
    genome.mutationRates["connections"] = dup.mutationRates["connections"]
    genome.mutationRates["link"] = dup.mutationRates["link"]
    genome.mutationRates["bias"] = dup.mutationRates["bias"]
    genome.mutationRates["node"] = dup.mutationRates["node"]
    genome.mutationRates["enable"] = dup.mutationRates["enable"]
    genome.mutationRates["disable"] = dup.mutationRates["disable"]
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


function duplicateGene(dup)
	local gene = newGene()
	gene.into = dup.into
	gene.out = dup.out
	gene.weight = dup.weight
	gene.enabled = dup.enabled
	gene.innovation = dup.innovation
	return gene
end

--new neuron constructor
function newNeuron()
   local neuron = {}
   neuron.genesInto = {}
   neuron.value = 0.0
   return neuron
end

function getRandomNeuron(genes, includeInputs)
    local choices = {}
    --inputs
    if includeInputs then 
        for i=1,numInputs do
            choices[i] = true
        end
    end
    --outputs
    for i=1,numOutputs do
        choices[maxNeurons+i] = true
    end
    --hidden
    for i,gene in pairs(genes) do
        if includeInputs or gene.into > numInputs then
            neurons[gene.into] = true
        end
        if includeInputs or gene.out > numInputs then
            neurons[gene.out] = true
        end
    end

    --choose a random number neuron from the options and return it
    local choice = math.random(1,#choices)
    for neuron,_ in pairs(choices) do
        choice = choice - 1
        if choice == 0 then
            return neuron
        end
    end
end

--new session constructor
function newSession()
    local session = {}
    --current species, etc storage
    session.curGeneration = 0
    session.species = {}
    session.curSpeciesNum = 1
    session.curGenomeNum = 1
    session.innovation = numOutputs
    --fitness evals
    session.curFrame = 0
    session.highestPoint = 0
    session.maxFitness = 0
    return session
end

function newInnovation()
   session.innovation = session.innovation + 1
   return session.innovation
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
    local found = false
    for _, species in pairs(session.species) do
        if fitsSpecies(genomeToGroup, species.genomes[1]) then
            table.insert(species.genomes, genomeToGroup) 
        end
    end

    --create a new species for the genome
    if not found then
        local species = newSpecies()
        table.insert(species.genomes, genomeToGroup)
        table.insert(session.species, species)
    end

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
            --if ns that the gene connects to are not yet created, create them
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
    
    
    --set the maximum direction output to true
    --since no direction button can be pressed at the same time
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

--get the next genome that has not been evaluated, 
--if none exist, evolve a new set of species
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

--kill off lower genomes in a generation
function removeLowPerformingGenomes(onlyKeepTop)
    for _,species in pairs(session.species) do
        --sort genomes in species by fitness
        table.sort(species.genomes, function(a,b) return a.fitness > b.fitness end )
        local numToKeep = math.ceil(#species.genomes/2)
        if onlyKeepTop then
            numToKeep = 1
        end
        while #species.genomes > remaining do
            table.remove(species.genomes)
        end
    end

end

--rank the genomes in terms of their fitness values
function rankGenomes()
    local sorted = {}
    for _,species in pairs(session.species) do
        for _,genome in pairs(species.genomes) do
            table.insert(sorted, genome)
        end
    end
    table.sort(sorted, function(a,b) return (a.fitness < b.fitness) end)
    for rank,genome in pairs(sorted) do
        genome.rank = rank
    end
end

--remove species that have been stuck for too long
function removeStuckSpecies()
    local notStuck = {}
    for _,species in pairs(session.species) do
        table.sort(species.genomes, function(a,b) return (a.fitness > b.fitness) end)
        species.numGensStuck = numGensStuck + 1
        if species.genome[1] > species.topFitness then
            species.topFitness = species.genomes[1].fitness
            species.numGensStuck = 0
        end
        if species.numGensStuck < maxGensStuck then
            table.insert(notStuck, species)
        end
    end
    session.species = notStuck
end

--remove species that have less than their share of the total fitness
function removeWeakSpecies()
    local strong = {}

    --get total rank sum
    local totalAverageRank = 0
    for _,species in pairs(session.species) do
        total = total + species.averageRank
    end

    --if a species contributes at least an even share of the expected fitness, it lives
    for _,species in pairs(session.species) do
        local rankShare = math.floor(species.averageRank / totalAverageRank * populationSize)
        if rankShare >= 1 then
            table.insert(strong, species)
        end
    end
    session.species = strong
end

function crossover(genome1, genome2)
    --make genome 1 the highest
    if genome2.fitness > genome1.fitness then
       temp = genome1
       genome1 = genome2
       genome2 = temp
    end

    --take paramters from stronger genome
    local crossoverChild = newGenome()
    crossoverChild.maxneuron = math.max(genome1.maxneuron, genome2.maxneuron)
    for mutationType, rate in pairs(genome1.mutationRates) do
        crossoverChild.mutationRates[mutationType] = rate
    end

    --map innovations to genes for lookup
    local weakInnovations = {}
    for _,gene in genome2.genes do
        weakInnovations[gene.innovation] = gene
    end

    --randomly choose matching genomes, or take stronger if not matching
    for _,gene in genome1.genes do
        local matchingGene = weakInnovations[gene.innovation]
        --choose the stronger genome
        if matchingGene == nil or not matchingGene.enabled or math.random(2) == 1 then
            table.insert(crossoverChild.genes, duplicateGene(gene))

        --choose the weaker genome
        else
            table.insert(crossoverChild.genes, duplicateGene(matchingGene))
        end
    end
    return crossoverChild
end

--mutate gene weights
function connectionMutate(genome)
    local maxStep = genome.mutationRates['step']
    for _,gene in pairs(genome.genes) do
        if math.random() < perturbChance then
            gene.weight = gene.weight + math.random()*maxStep*2 - step
        else
            gene.weight = math.random()*4-2
        end
    end
end


--randomly add a gene between two neurons
function addGene(genome, fromBias)
    --get two random neurons
    local randomIn = randomNeuron(genome.genes, true)
    local randomOut = randomNeuron(genome.genes, false)

    --connect the neurons with a gene, conditionally from bias
    local newGene = newGene()
    newGene.into = randomIn
    newGene.out = randomOut
    if fromBias then
        newGene.into = numInputs
    end

    --gene already exists check
    for _,gene in pairs(genome.genes) do
        if gene.into == newGene.into and gene.out == newGene.out then
            return
        end
    end

    newGene.innovation = nextInnovation()
    newGene.weight = math.random()*4-2
    
    table.insert(genome.genes, newGene)
end

--mutate a neuron
function neuronMutate(genome)
    --no genes to mutate from or disabled
    if #genome.genes == 0 then
        return
    end
    local randomGene = genome.genes[math.random(1, #genome.genes)]
    if not randomGene.enabled then
        return
    end
    randomGene.enabled = false

    local g1 = duplicateGene(randomGene)
    g1.out = genome.maxneuron
    g1.weight = 1.0
    g1.innovation = newInnovation()
    g1.enabled = true
    table.insert(genome.genes, g1)

    local g2 = duplicateGene(randomGene)
    g1.into = genome.maxneuron
    g1.innovation = newInnovation()
    g1.enabled = true
    table.insert(genome.genes, g1)

end

--randomly toggle a enabled
function enableMutate(genome, doEnable)
    local toggleChoices = {}
    for _,gene in pairs(genome.genes) do
        if not gene.enabled == doEnable then
            table.insert(toggleChoices, gene)
        end
    end

    if #toggleChoices == 0 then
        return
    end

    local randomGene = toggleChoices[math.random(1,#toggleChoices)]
    randomGene.enabled = not randomGene.enabled

end

--handles all the possible mutations and whether or not they should happen
function mutate(genome)
    for mutationType, rate in pairs(genome.mutationRates) do
        if math.random(2) == 1 then
            genome.mutationRates[mutationType] = rate*0.95
        else
            genome.mutationRates[mutationType] = rate/0.95
        end
    end

    --mutate connections
    if math.random() < genome.mutationRates["connections"] then
        connectionMutate(genome)
    end

    --mutate to add gene
    local mutateChance = genome.mutationRates["link"]
    while mutateChance > 0 do
        if math.random() < mutateChance then
            addGene(genome, false)
        end
        mutateChance = mutateChance - 1
    end

    --mutate to add gene from bias
    mutateChance = genome.mutationRates["bias"]
    while mutateChance > 0 do
        if math.random() < mutateChance then
            addGene(genome, true)
        end
        mutateChance = mutateChance - 1
    end

    --mutate neurons 
    mutateChance = genome.mutationRates["node"]
    while mutateChance > 0 do
        if math.random() < mutateChance then
            neuronMutate(genome)
        end
        mutateChance = mutateChance - 1
    end

    --mutate enables
    mutateChance = genome.mutationRates["enable"]
    while mutateChance > 0 do
        if math.random() < mutateChance then
            enableMutate(genome, true)
        end
        mutateChance = mutateChance - 1
    end

    --mutate disables
    mutateChance = genome.mutationRates["disable"]
    while mutateChance > 0 do
        if math.random() < mutateChance then
            enableMutate(genome, false)
        end
        mutateChance = mutateChance - 1
    end

end

--creates a child with roots to the species arg
function createChild(parentSpecies)
    local child = {}
    if math.random() < crossoverChance then
        genome1 = species.genomes[math.random(1,#species.genomes)]
        genome2 = species.genomes[math.random(1,#species.genomes)]
        child = crossover(genome1, genome2)
    else
        genome = species.genomes[math.random(1,#species.genomes)]
        child = duplicateGenome(genome)
    end

    mutate(child)
    return child
end

--evolve the next genration of genomes
function evolveNewGeneration()
    --cull bottom half of each species
    removeLowPerformingGenomes(false)
    --remove stuck species
    removeStuckSpecies()
    --rank the species
    rankGenomes() 
    --calculate the average fitness of each species
    for _,species in session.species do
        local speciesTotal = 0
        for _,genome in species.genomes do
            speciesTotal = speciesTotal + genome.rank
        end
        species.averageRank = speciesTotal / #species.genomes
    end
    --remove weaker species
    removeWeakSpecies()
    
    --get total fitness value average on remaining strong species
    local totalAverageRank = 0
    for _,species in session.species do
        total = total + species.averageRank
    end

    --breed new children
    local children = {}
    for _,species in session.species do
        numToBreed = math.floor(species.averageRank * populationSize / totalAverageRank) - 1
        for _=1,numToBreed do
            table.insert(children, createChild(species))
        end
    end
    --cull all but top of each species
    removeLowPerformingGenomes(true)
    while #children+#session.species < populationSize do
        local randomSpecies = session.species[math.random(1,#session.species)]
        table.insert(children, createChild(randomSpecies))
    end
    for _,child in children do
        groupGenome(child)
    end
    session.generationNum = session.generationNum + 1

end

--get the fitness value given mario's location
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

--get the total percentage of genomes that have already been evaluated
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

