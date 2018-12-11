--------------------------------------
--drawing constants ------------------
--------------------------------------
viewDist = 10
scale = 8
inputDisplayOffset = 10

--------------------------------------
--reward values ----------------------
--------------------------------------
framesElapsed = 0 
--------------------------------------
function round(num)
    lower = math.floor(num)
    upper = math.floor(num) + 1
    lowerV = -(lower - num)
    upperV = upper - num
    if (upperV > lowerV) then
        return lower
    else
        return upper
    end
end


function outputsToController(outputs)
	controller = {}
	for button, action in pairs(outputs) do
		controller["P1 " .. button] = action
	end
	joypad.set(controller)
end

function getMarioLocation()
	local marioX = memory.readbyte(0x0203)+6
	local marioY = memory.readbyte(0x0204)-3
	gui.drawPixel(marioX, marioY, 0xFF00FF00)
	return {['x']=marioX, ['y']=marioY}
end

function getSpriteLocations()
	local sprites = {}
	--barrels
	for i=0x0233, 0x02C3, 0x10 do
		x = memory.readbyte(i)
		y = memory.readbyte(i+1)
		sprites[#sprites+1] = {['x']=x+2, ['y']=y-4}
	end
	--fire dudes
	for i=0x0213, 0x0232, 0x10 do
		x = memory.readbyte(i)
		y = memory.readbyte(i+1)
		sprites[#sprites+1] = {['x']=x+2, ['y']=y-4}
	end
	--hammers
	for i=0x02CB, 0x02DB, 0x08 do
		x = memory.readbyte(i)
		y = memory.readbyte(i+1)
		--sprites[#sprites+1] = {['x']=x+2, ['y']=y-4}
	end
	--barrel fire location
	sprites[#sprites+1] = {['x']=32, ['y']=192}

	return sprites
end


function hasWon()
	local marioZone = memory.readbyte(0x0400)
	return marioZone ~= 1
end


function isDead()
	local marioStatus = memory.readbyte(0X0096)
	return marioStatus == 0xFF
end

function getPlatformLocations()
	local platforms = {}
	for i=0x808C, 0x80BB, 0x06 do
		x = memory.readbyte(i)
		y = memory.readbyte(i+1)
		offset = memory.readbyte(i+2)
		platSize = memory.readbyte(0x80C3 + 2 + offset)
		hsd = memory.readbyte(i+3)
		numSections = memory.readbyte(i+4) - 1
		if i==0x808C+6 then
			numSections = numSections - 1
		end
		for j=0,numSections do
			if hsd == 0x18 then
				xp = x + platSize*j
			else
				xp = x - platSize*j
			end
			yp = y - j
			for k=0,(platSize/0x8)-1 do
				platforms[#platforms+1] = {['x']=xp+k*8,['y']=yp}
			end
		end
	end
	
	return platforms
end
	
function getInputs()
	local inputGrid = {}
	for i=0,viewDist*2 do
		inputGrid[i] = {}
		for j=0, viewDist*2 do
			inputGrid[i][j] = 0
		end
	end
	
	mario = getMarioLocation()
	platforms = getPlatformLocations()
	for _, plat in pairs(platforms) do
		x = plat['x']
		y = plat['y']
		x_diff = x-mario['x']
		y_diff = y-mario['y']
		if math.abs(x_diff) <= viewDist*0x8 and
		   math.abs(y_diff) <= viewDist*0x8 then
		   ix = round(x_diff/0x8) + viewDist
		   iy = round(y_diff/0x8) + viewDist
		   inputGrid[ix][iy] = 1
		end
	end
	ladders = getLadderLocations()
	for _, ladder in pairs(ladders) do
		x = ladder['x']
		for y=ladder['y'], ladder['y']+ladder['sizeY']+1, scale do
			x_diff = x-mario['x']
			y_diff = y-mario['y']
			if math.abs(x_diff) <= viewDist*0x8 and
			   math.abs(y_diff) <= viewDist*0x8 then
			   ix = round(x_diff/0x8) + viewDist
			   iy = round(y_diff/0x8) + viewDist
			   inputGrid[ix][iy] = 1
			end
		end
	end

	sprites = getSpriteLocations()
	for _, sprite in pairs(sprites) do
		x = sprite['x']
		y = sprite['y']
		if y < 240 then --check if the sprite is on screen (enabled)
			x_diff = x-mario['x']
			y_diff = y-mario['y']
			if math.abs(x_diff) <= viewDist*0x8 and
			   math.abs(y_diff) <= viewDist*0x8 then
			   ix = round(x_diff/0x8) + viewDist
			   iy = round(y_diff/0x8) + viewDist
			   inputGrid[ix][iy] = -1
			end	
		end
	end 
	inputs = {}
	for i=0,viewDist*2 do
		for j=0, viewDist*2 do
			inputs[#inputs+1] = inputGrid[j][i]
		end
	end
	
	return inputs
end

function drawInputScreen(inputs)
	gui.drawBox(inputDisplayOffset-1, inputDisplayOffset-1,
	           (viewDist*0x8*2)/scale+1+inputDisplayOffset, 
			   (viewDist*0x8*2)/scale+1+inputDisplayOffset,
			   0xFFFF0000, 0x80000000)
	gui.drawPixel(inputDisplayOffset + viewDist, 
				  inputDisplayOffset + viewDist,
				  0xFFFF0000)
	for x=0,viewDist*2 do
		for y=0, viewDist*2 do
			value = inputs[x][y]
			if value == 1 then
				gui.drawPixel(x+inputDisplayOffset,
							  y+inputDisplayOffset,
							  0xFFFFFFFF)
			elseif value == -1 then
				gui.drawPixel(x+inputDisplayOffset,
							  y+inputDisplayOffset,
							  0xFF0000FF)
			end
		end
	end
end

function getLadderLocations()
	local ladders = {}
	local cnt = 0
	for i=0x80E3, 0x810A, 0x03 do
		if cnt ~= 5 and cnt ~= 8 and cnt ~= 2 and cnt ~= 3 then
			x = memory.readbyte(i)
			if x < 240 then
				y = memory.readbyte(i+1)
				--gui.drawText(x, y, #ladders, 0xFFFFFFFF, 9)
				offset = memory.readbyte(i+2)
				ladderSizeX = memory.readbyte(0x810B + 2 + offset)
				ladderSizeY = memory.readbyte(0x810B + 3 + offset)
				ladders[#ladders + 1] = {['x']=x,['y']=y, ['sizeX']=ladderSizeX, ['sizeY']=ladderSizeY}
			end
		end
		cnt = cnt + 1
	end
	return ladders
end

function drawLadders(ladders)
	for _,l in pairs(ladders) do
		gui.drawBox(l['x'], l['y'] ,l['x']+l['sizeX'], l['y']+l['sizeY'], 0xFF00FF00, 0xFF00FF00)
	end
end

--------------------------------------
--neural network constants
--------------------------------------
saveStateFile = "DP1.state"

numPixels = (viewDist*2+1)*(viewDist*2+1)

numInputs = numPixels+1
controllerButtons = {"A", "Left", "Right", "Up", "Down"}
numOutputs = #controllerButtons

Population = 300
DeltaDisjoint = 2.0
DeltaWeights = 0.4
DeltaThreshold = 1.0

StaleSpecies = 15

MutateConnectionsChance = 0.25
PerturbChance = 0.90
CrossoverChance = 0.75
LinkMutationChance = 2.0
NodeMutationChance = 0.50
BiasMutationChance = 0.40
StepSize = 0.1
DisableMutationChance = 0.4
EnableMutationChance = 0.2

TimeoutConstant = 300

MaxNodes = 1000000

function sigmoid(x)
	return 2/(1+math.exp(-4.9*x))-1
end

function newInnovation()
	pool.innovation = pool.innovation + 1
	return pool.innovation
end

function newPool()
	local pool = {}
	pool.species = {}
	pool.generation = 0
	pool.innovation = numOutputs
	pool.currentSpecies = 1
	pool.currentGenome = 1
	pool.currentFrame = 0
	pool.maxFitness = 0
	
	return pool
end

function newSpecies()
	local species = {}
	species.topFitness = 0
	species.staleness = 0
	species.genomes = {}
	species.averageFitness = 0
	
	return species
end

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

function copyGenome(genome)
	local genome2 = newGenome()
	for g=1,#genome.genes do
		table.insert(genome2.genes, copyGene(genome.genes[g]))
	end
	genome2.maxneuron = genome.maxneuron
	genome2.mutationRates["connections"] = genome.mutationRates["connections"]
	genome2.mutationRates["link"] = genome.mutationRates["link"]
	genome2.mutationRates["bias"] = genome.mutationRates["bias"]
	genome2.mutationRates["node"] = genome.mutationRates["node"]
	genome2.mutationRates["enable"] = genome.mutationRates["enable"]
	genome2.mutationRates["disable"] = genome.mutationRates["disable"]
	
	return genome2
end

function basicGenome()
	local genome = newGenome()
	local innovation = 1

	genome.maxneuron = numInputs
	mutate(genome)
	
	return genome
end

function newGene()
	local gene = {}
	gene.into = 0
	gene.out = 0
	gene.weight = 0.0
	gene.enabled = true
	gene.innovation = 0
	
	return gene
end

function copyGene(gene)
	local gene2 = newGene()
	gene2.into = gene.into
	gene2.out = gene.out
	gene2.weight = gene.weight
	gene2.enabled = gene.enabled
	gene2.innovation = gene.innovation
	
	return gene2
end

function newNeuron()
	local neuron = {}
	neuron.incoming = {}
	neuron.value = 0.0
	
	return neuron
end

function generateNetwork(genome)
	local network = {}
	network.neurons = {}
	
	for i=1,numInputs do
		network.neurons[i] = newNeuron()
	end
	
	for o=1,numOutputs do
		network.neurons[MaxNodes+o] = newNeuron()
	end
	
	table.sort(genome.genes, function (a,b)
		return (a.out < b.out)
	end)
	for i=1,#genome.genes do
		local gene = genome.genes[i]
		if gene.enabled then
			if network.neurons[gene.out] == nil then
				network.neurons[gene.out] = newNeuron()
			end
			local neuron = network.neurons[gene.out]
			table.insert(neuron.incoming, gene)
			if network.neurons[gene.into] == nil then
				network.neurons[gene.into] = newNeuron()
			end
		end
	end
	
	genome.network = network
end

function evaluateNetwork(network, inputs)
	table.insert(inputs, 1)
	if #inputs ~= numInputs then
		console.writeline("Incorrect number of neural network inputs.")
		return {}
	end
	
	for i=1,numInputs do
		network.neurons[i].value = inputs[i]
	end
	
	for _,neuron in pairs(network.neurons) do
		local sum = 0
		for j = 1,#neuron.incoming do
			local incoming = neuron.incoming[j]
			local other = network.neurons[incoming.into]
			sum = sum + incoming.weight * other.value
		end
		
		if #neuron.incoming > 0 then
			neuron.value = sigmoid(sum)
		end
	end
	
	local outputs = {}
	for o=1,numOutputs do
		local button = "P1 " .. controllerButtons[o]
		if network.neurons[MaxNodes+o].value > 0 then
			outputs[button] = true
		else
			outputs[button] = false
		end
	end
	
	return outputs
end

function crossover(g1, g2)
	-- Make sure g1 is the higher fitness genome
	if g2.fitness > g1.fitness then
		tempg = g1
		g1 = g2
		g2 = tempg
	end

	local child = newGenome()
	
	local innovations2 = {}
	for i=1,#g2.genes do
		local gene = g2.genes[i]
		innovations2[gene.innovation] = gene
	end
	
	for i=1,#g1.genes do
		local gene1 = g1.genes[i]
		local gene2 = innovations2[gene1.innovation]
		if gene2 ~= nil and math.random(2) == 1 and gene2.enabled then
			table.insert(child.genes, copyGene(gene2))
		else
			table.insert(child.genes, copyGene(gene1))
		end
	end
	
	child.maxneuron = math.max(g1.maxneuron,g2.maxneuron)
	
	for mutation,rate in pairs(g1.mutationRates) do
		child.mutationRates[mutation] = rate
	end
	
	return child
end

function randomNeuron(genes, nonInput)
	local neurons = {}
	if not nonInput then
		for i=1,numInputs do
			neurons[i] = true
		end
	end
	for o=1,numOutputs do
		neurons[MaxNodes+o] = true
	end
	for i=1,#genes do
		if (not nonInput) or genes[i].into > numInputs then
			neurons[genes[i].into] = true
		end
		if (not nonInput) or genes[i].out > numInputs then
			neurons[genes[i].out] = true
		end
	end

	local count = 0
	for _,_ in pairs(neurons) do
		count = count + 1
	end
	local n = math.random(1, count)
	
	for k,v in pairs(neurons) do
		n = n-1
		if n == 0 then
			return k
		end
	end
	
	return 0
end

function containsLink(genes, link)
	for i=1,#genes do
		local gene = genes[i]
		if gene.into == link.into and gene.out == link.out then
			return true
		end
	end
end

function pointMutate(genome)
	local step = genome.mutationRates["step"]
	
	for i=1,#genome.genes do
		local gene = genome.genes[i]
		if math.random() < PerturbChance then
			gene.weight = gene.weight + math.random() * step*2 - step
		else
			gene.weight = math.random()*4-2
		end
	end
end

function linkMutate(genome, forceBias)
	local neuron1 = randomNeuron(genome.genes, false)
	local neuron2 = randomNeuron(genome.genes, true)
	 
	local newLink = newGene()
	if neuron1 <= numInputs and neuron2 <= numInputs then
		--Both input nodes
		return
	end
	if neuron2 <= numInputs then
		-- Swap output and input
		local temp = neuron1
		neuron1 = neuron2
		neuron2 = temp
	end

	newLink.into = neuron1
	newLink.out = neuron2
	if forceBias then
		newLink.into = numInputs
	end
	
	if containsLink(genome.genes, newLink) then
		return
	end
	newLink.innovation = newInnovation()
	newLink.weight = math.random()*4-2
	
	table.insert(genome.genes, newLink)
end

function nodeMutate(genome)
	if #genome.genes == 0 then
		return
	end

	genome.maxneuron = genome.maxneuron + 1

	local gene = genome.genes[math.random(1,#genome.genes)]
	if not gene.enabled then
		return
	end
	gene.enabled = false
	
	local gene1 = copyGene(gene)
	gene1.out = genome.maxneuron
	gene1.weight = 1.0
	gene1.innovation = newInnovation()
	gene1.enabled = true
	table.insert(genome.genes, gene1)
	
	local gene2 = copyGene(gene)
	gene2.into = genome.maxneuron
	gene2.innovation = newInnovation()
	gene2.enabled = true
	table.insert(genome.genes, gene2)
end

function enableDisableMutate(genome, enable)
	local candidates = {}
	for _,gene in pairs(genome.genes) do
		if gene.enabled == not enable then
			table.insert(candidates, gene)
		end
	end
	
	if #candidates == 0 then
		return
	end
	
	local gene = candidates[math.random(1,#candidates)]
	gene.enabled = not gene.enabled
end

function mutate(genome)
	for mutation,rate in pairs(genome.mutationRates) do
		if math.random(1,2) == 1 then
			genome.mutationRates[mutation] = 0.95*rate
		else
			genome.mutationRates[mutation] = 1.05263*rate
		end
	end

	if math.random() < genome.mutationRates["connections"] then
		pointMutate(genome)
	end
	
	local p = genome.mutationRates["link"]
	while p > 0 do
		if math.random() < p then
			linkMutate(genome, false)
		end
		p = p - 1
	end

	p = genome.mutationRates["bias"]
	while p > 0 do
		if math.random() < p then
			linkMutate(genome, true)
		end
		p = p - 1
	end
	
	p = genome.mutationRates["node"]
	while p > 0 do
		if math.random() < p then
			nodeMutate(genome)
		end
		p = p - 1
	end
	
	p = genome.mutationRates["enable"]
	while p > 0 do
		if math.random() < p then
			enableDisableMutate(genome, true)
		end
		p = p - 1
	end

	p = genome.mutationRates["disable"]
	while p > 0 do
		if math.random() < p then
			enableDisableMutate(genome, false)
		end
		p = p - 1
	end
end

function disjoint(genes1, genes2)
	local i1 = {}
	for i = 1,#genes1 do
		local gene = genes1[i]
		i1[gene.innovation] = true
	end

	local i2 = {}
	for i = 1,#genes2 do
		local gene = genes2[i]
		i2[gene.innovation] = true
	end
	
	local disjointGenes = 0
	for i = 1,#genes1 do
		local gene = genes1[i]
		if not i2[gene.innovation] then
			disjointGenes = disjointGenes+1
		end
	end
	
	for i = 1,#genes2 do
		local gene = genes2[i]
		if not i1[gene.innovation] then
			disjointGenes = disjointGenes+1
		end
	end
	
	local n = math.max(#genes1, #genes2)
	
	return disjointGenes / n
end

function weights(genes1, genes2)
	local i2 = {}
	for i = 1,#genes2 do
		local gene = genes2[i]
		i2[gene.innovation] = gene
	end

	local sum = 0
	local coincident = 0
	for i = 1,#genes1 do
		local gene = genes1[i]
		if i2[gene.innovation] ~= nil then
			local gene2 = i2[gene.innovation]
			sum = sum + math.abs(gene.weight - gene2.weight)
			coincident = coincident + 1
		end
	end
	
	return sum / coincident
end

function sameSpecies(genome1, genome2)
	local dd = DeltaDisjoint*disjoint(genome1.genes, genome2.genes)
	local dw = DeltaWeights*weights(genome1.genes, genome2.genes) 
	return dd + dw < DeltaThreshold
end

function rankGlobally()
	local global = {}
	for s = 1,#pool.species do
		local species = pool.species[s]
		for g = 1,#species.genomes do
			table.insert(global, species.genomes[g])
		end
	end
	table.sort(global, function (a,b)
		return (a.fitness < b.fitness)
	end)
	
	for g=1,#global do
		global[g].globalRank = g
	end
end

function calculateAverageFitness(species)
	local total = 0
	
	for g=1,#species.genomes do
		local genome = species.genomes[g]
		total = total + genome.globalRank
	end
	
	species.averageFitness = total / #species.genomes
end

function totalAverageFitness()
	local total = 0
	for s = 1,#pool.species do
		local species = pool.species[s]
		total = total + species.averageFitness
	end

	return total
end

function cullSpecies(cutToOne)
	for s = 1,#pool.species do
		local species = pool.species[s]
		
		table.sort(species.genomes, function (a,b)
			return (a.fitness > b.fitness)
		end)
		
		local remaining = math.ceil(#species.genomes/2)
		if cutToOne then
			remaining = 1
		end
		while #species.genomes > remaining do
			table.remove(species.genomes)
		end
	end
end

function breedChild(species)
	local child = {}
	if math.random() < CrossoverChance then
		g1 = species.genomes[math.random(1, #species.genomes)]
		g2 = species.genomes[math.random(1, #species.genomes)]
		child = crossover(g1, g2)
	else
		g = species.genomes[math.random(1, #species.genomes)]
		child = copyGenome(g)
	end
	
	mutate(child)
	
	return child
end

function removeStaleSpecies()
	local survived = {}

	for s = 1,#pool.species do
		local species = pool.species[s]
		
		table.sort(species.genomes, function (a,b)
			return (a.fitness > b.fitness)
		end)
		
		if species.genomes[1].fitness > species.topFitness then
			species.topFitness = species.genomes[1].fitness
			species.staleness = 0
		else
			species.staleness = species.staleness + 1
		end
		if species.staleness < StaleSpecies or species.topFitness >= pool.maxFitness then
			table.insert(survived, species)
		end
	end

	pool.species = survived
end

function removeWeakSpecies()
	local survived = {}

	local sum = totalAverageFitness()
	for s = 1,#pool.species do
		local species = pool.species[s]
		breed = math.floor(species.averageFitness / sum * Population)
		if breed >= 1 then
			table.insert(survived, species)
		end
	end

	pool.species = survived
end


function addToSpecies(child)
	local foundSpecies = false
	for s=1,#pool.species do
		local species = pool.species[s]
		if not foundSpecies and sameSpecies(child, species.genomes[1]) then
			table.insert(species.genomes, child)
			foundSpecies = true
		end
	end
	
	if not foundSpecies then
		local childSpecies = newSpecies()
		table.insert(childSpecies.genomes, child)
		table.insert(pool.species, childSpecies)
	end
end

function newGeneration()
	cullSpecies(false) -- Cull the bottom half of each species
	rankGlobally()
	removeStaleSpecies()
	rankGlobally()
	for s = 1,#pool.species do
		local species = pool.species[s]
		calculateAverageFitness(species)
	end
	removeWeakSpecies()
	local sum = totalAverageFitness()
	local children = {}
	for s = 1,#pool.species do
		local species = pool.species[s]
		breed = math.floor(species.averageFitness / sum * Population) - 1
		for i=1,breed do
			table.insert(children, breedChild(species))
		end
	end
	cullSpecies(true) -- Cull all but the top member of each species
	while #children + #pool.species < Population do
		local species = pool.species[math.random(1, #pool.species)]
		table.insert(children, breedChild(species))
	end
	for c=1,#children do
		local child = children[c]
		addToSpecies(child)
	end
	
	pool.generation = pool.generation + 1
	
	--writeFile("backup." .. pool.generation .. "." .. forms.gettext(saveLoadFile))
end
	
function initializePool()
	pool = newPool()

	for i=1,Population do
		basic = basicGenome()
		addToSpecies(basic)
	end

	initializeRun()
end

function clearJoypad()
	controller = {}
	for b = 1,#controllerButtons do
		controller["P1 " .. controllerButtons[b]] = false
	end
	joypad.set(controller)
end

function initializeRun()
	savestate.load(saveStateFile);
	uppermost = 0
	pool.currentFrame = 0
	timeout = TimeoutConstant
	clearJoypad()
	
	local species = pool.species[pool.currentSpecies]
	local genome = species.genomes[pool.currentGenome]
	generateNetwork(genome)
	evaluateCurrent()
	print("I did it mom!")
end

function evaluateCurrent()
	local species = pool.species[pool.currentSpecies]
	local genome = species.genomes[pool.currentGenome]

	inputs = getInputs()
	controller = evaluateNetwork(genome.network, inputs)
	
	if controller["P1 Left"] and controller["P1 Right"] then
		controller["P1 Left"] = false
		controller["P1 Right"] = false
	end
	if controller["P1 Up"] and controller["P1 Down"] then
		controller["P1 Up"] = false
		controller["P1 Down"] = false
	end

	joypad.set(controller)
end

if pool == nil then
	initializePool()
end


function nextGenome()
	pool.currentGenome = pool.currentGenome + 1
	if pool.currentGenome > #pool.species[pool.currentSpecies].genomes then
		pool.currentGenome = 1
		pool.currentSpecies = pool.currentSpecies+1
		if pool.currentSpecies > #pool.species then
			newGeneration()
			pool.currentSpecies = 1
		end
	end
end

function fitnessAlreadyMeasured()
	local species = pool.species[pool.currentSpecies]
	local genome = species.genomes[pool.currentGenome]
	
	return genome.fitness ~= 0
end

function displayGenome(genome)
	local network = genome.network
	local cells = {}
	local i = 1
	local cell = {}
	for dy=-viewDist,viewDist do
		for dx=-viewDist,viewDist do
			cell = {}
			cell.x = 50+5*dx
			cell.y = 70+5*dy
			cell.value = network.neurons[i].value
			cells[i] = cell
			i = i + 1
		end
	end
	local biasCell = {}
	biasCell.x = 80
	biasCell.y = 110
	biasCell.value = network.neurons[numInputs].value
	cells[numInputs] = biasCell
	
	for o = 1,numOutputs do
		cell = {}
		cell.x = 220
		cell.y = 30 + 8 * o
		cell.value = network.neurons[MaxNodes + o].value
		cells[MaxNodes+o] = cell
		local color
		if cell.value > 0 then
			color = 0xFF0000FF
		else
			color = 0xFF000000
		end
		gui.drawText(223, 24+8*o, controllerButtons[o], color, 9)
	end
	
	for n,neuron in pairs(network.neurons) do
		cell = {}
		if n > numInputs and n <= MaxNodes then
			cell.x = 140
			cell.y = 40
			cell.value = neuron.value
			cells[n] = cell
		end
	end
	
	for n=1,4 do
		for _,gene in pairs(genome.genes) do
			if gene.enabled then
				local c1 = cells[gene.into]
				local c2 = cells[gene.out]
				if gene.into > numInputs and gene.into <= MaxNodes then
					c1.x = 0.75*c1.x + 0.25*c2.x
					if c1.x >= c2.x then
						c1.x = c1.x - 40
					end
					if c1.x < 90 then
						c1.x = 90
					end
					
					if c1.x > 220 then
						c1.x = 220
					end
					c1.y = 0.75*c1.y + 0.25*c2.y
					
				end
				if gene.out > numInputs and gene.out <= MaxNodes then
					c2.x = 0.25*c1.x + 0.75*c2.x
					if c1.x >= c2.x then
						c2.x = c2.x + 40
					end
					if c2.x < 90 then
						c2.x = 90
					end
					if c2.x > 220 then
						c2.x = 220
					end
					c2.y = 0.25*c1.y + 0.75*c2.y
				end
			end
		end
	end
	
	gui.drawBox(50-viewDist*5-3,70-viewDist*5-3,50+viewDist*5+2,70+viewDist*5+2,0xFF000000, 0x80808080)
	for n,cell in pairs(cells) do
		if n > numInputs or cell.value ~= 0 then
			local color = math.floor((cell.value+1)/2*256)
			if color > 255 then color = 255 end
			if color < 0 then color = 0 end
			local opacity = 0xFF000000
			if cell.value == 0 then
				opacity = 0x50000000
			end
			color = opacity + color*0x10000 + color*0x100 + color
			gui.drawBox(cell.x-2,cell.y-2,cell.x+2,cell.y+2,opacity,color)
		end
	end
	for _,gene in pairs(genome.genes) do
		if gene.enabled then
			local c1 = cells[gene.into]
			local c2 = cells[gene.out]
			local opacity = 0xA0000000
			if c1.value == 0 then
				opacity = 0x20000000
			end
			
			local color = 0x80-math.floor(math.abs(sigmoid(gene.weight))*0x80)
			if gene.weight > 0 then 
				color = opacity + 0x8000 + 0x10000*color
			else
				color = opacity + 0x800000 + 0x100*color
			end
			gui.drawLine(c1.x+1, c1.y, c2.x-3, c2.y, color)
		end
	end
	
	gui.drawBox(49,71,51,78,0x00000000,0x80FF0000)
	
	if forms.ischecked(showMutationRates) then
		local pos = 100
		for mutation,rate in pairs(genome.mutationRates) do
			gui.drawText(100, pos, mutation .. ": " .. rate, 0xFF000000, 10)
			pos = pos + 8
		end
	end
end

--************************************
--**---main state execution --------**
--************************************


form = forms.newform(200, 260, "Fitness")
maxFitnessLabel = forms.label(form, "Max Fitness: " .. math.floor(pool.maxFitness), 5, 8)
showNetwork = forms.checkbox(form, "Show Map", 5, 30)
showMutationRates = forms.checkbox(form, "Show M-Rates", 5, 52)
restartButton = forms.button(form, "Restart", initializePool, 5, 77)
--saveButton = forms.button(form, "Save", savePool, 5, 102)
--loadButton = forms.button(form, "Load", loadPool, 80, 102)
--saveLoadFile = forms.textbox(form, Filename .. ".pool", 170, 25, nil, 5, 148)
--saveLoadLabel = forms.label(form, "Save/Load:", 5, 129)
--playTopButton = forms.button(form, "Play Top", playTop, 5, 170)
hideBanner = forms.checkbox(form, "Hide Banner", 5, 190)

--outputs = {["A"]=false, 
--		   ["Left"]=false,
--		   ["Right"]=false,
--		   ["Up"]=false,
--		   ["Down"]=false}

function onExit()
	forms.destroy(form)
end
event.onexit(onExit)

while true do
	local backgroundColor = 0xD0FFFFFF
	if not forms.ischecked(hideBanner) then
		gui.drawBox(0, 10, 300, 36, backgroundColor, backgroundColor)
	end

	local species = pool.species[pool.currentSpecies]
	local genome = species.genomes[pool.currentGenome]
	
	if forms.ischecked(showNetwork) then
		displayGenome(genome)
	end
	
	if pool.currentFrame%5 == 0 then
		evaluateCurrent()
	end

	joypad.set(controller)

	marioLoc = getMarioLocation()
	isJumping = memory.readbyte(0x96) == 0x04 or memory.readbyte(0x96) == 0xFF
	if (255-marioLoc['y']-51)*10 > uppermost and not isJumping then
		uppermost = (255-marioLoc['y']-51)*10
		timeout = TimeoutConstant
	end
	
	timeout = timeout - 1
	
	local timeoutBonus = pool.currentFrame / 4
	if timeout + timeoutBonus <= 0 or isDead() then
		local fitness = uppermost - pool.currentFrame / 150.0
		if gameinfo.getromname() == "Super Mario World (USA)" and uppermost > 4816 then
			fitness = fitness + 1000
		end
		if fitness < 0 then
			fitness = -1
		end
		genome.fitness = fitness
		
		if fitness > pool.maxFitness then
			pool.maxFitness = fitness
			forms.settext(maxFitnessLabel, "Max Fitness: " .. math.floor(pool.maxFitness))
			--writeFile("backup." .. pool.generation .. "." .. forms.gettext(saveLoadFile))
		end
		
		console.writeline("Gen " .. pool.generation .. " species " .. pool.currentSpecies .. " genome " .. pool.currentGenome .. " fitness: " .. fitness)
		pool.currentSpecies = 1
		pool.currentGenome = 1
		while fitnessAlreadyMeasured() do
			nextGenome()
		end
		initializeRun()
	end

	local measured = 0
	local total = 0
	for _,species in pairs(pool.species) do
		for _,genome in pairs(species.genomes) do
			total = total + 1
			if genome.fitness ~= 0 then
				measured = measured + 1
			end
		end
	end
	if not forms.ischecked(hideBanner) then
		gui.drawText(0, 10, "Gen " .. pool.generation .. " species " .. pool.currentSpecies .. " genome " .. pool.currentGenome .. " (" .. math.floor(measured/total*100) .. "%)", 0xFF000000, 11)
		gui.drawText(0, 22, "Fitness: " .. math.floor(uppermost - (pool.currentFrame) / 150.0 - (timeout + timeoutBonus)*0/45), 0xFF000000, 11)
		gui.drawText(100, 22, "Max Fitness: " .. math.floor(pool.maxFitness), 0xFF000000, 11)
	end
		
	pool.currentFrame = pool.currentFrame + 1

	emu.frameadvance();
end





-----GAME LOOP --------
while true do
	
	inputs = getInputs()
	drawInputScreen(inputs)
	
	if forms.ischecked(DEBUG) then
		gui.drawText(170, 5, "Frame: " .. framesElapsed, 0xFFFFFFFF, 4)
		for _, s in pairs(sprites) do
			gui.drawText(s['x']+2, s['y'] - 4, "X",0xFF0000FF, 4)
		end
	end
	
	--outputsToController(outputs)
	
	if hasWon() or isDead() then
		savestate.load(saveStateFile);
		framesElapsed = 0
	end
	
	nextFrame()
end
