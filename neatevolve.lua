-- MarI/O by SethBling
-- Feel free to use this code, but please do not redistribute it.
-- Intended for use with the BizHawk emulator and Super Mario World or Super Mario Bros. ROM.
-- For SMW, make sure you have a save state named "DP1.state" at the beginning of a level,
-- and put a copy in both the Lua folder and the root directory of BizHawk.
--
-- Enhanced MarI/O by Kevin Sun and Lior Vansteenskie
-- The power up block of memory is 2 bits, which means it can hold 4 states.
-- The power up data is located at 0x0019. The following are the values
-- corresponding to the powerup that mario gets when these values are
-- set.
PowerUpType_Normal = 0
PowerUpType_Mushroom = 1
PowerUpType_Cape = 2
PowerUpType_Fire = 3

PowerUpType_Curr = PowerUpType_Normal;

function SetTestname(testname)
    Testname = testname
end

function TogglePowerUp(PowerUpType)
	PowerUpType_Curr = PowerUpType;
end

if gameinfo.getromname() == "Super Mario World (USA)" then
	Filename = "DP1.state"
	ButtonNames = {
		"A",
		"B",
		"X",
		"Y",
		"Up",
		"Down",
		"Left",
		"Right",
	}
elseif gameinfo.getromname() == "Super Mario Bros." then
	Filename = "SMB1-1.state"
	ButtonNames = {
		"A",
		"B",
		"Up",
		"Down",
		"Left",
		"Right",
	}
end

--initialize later. boolean array of all the sprites. True if the sprite is
--an enemy. Goes from 0x00 to 0xFF (0-255)
Sprites = {}
--List of Sprites which we do not consider harmful.
--make sure this list is sorted before initialization.
GoodSprites = {0x0E, 0x21, 0x2C, 0x2D, 0x2F, 0x35, 0x3E, 0x41, 0x42, 0x43, 0x45,
                0x47, 0x48, 0x49, 0x4A, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58,
                0x59, 0x5A, 0x5B, 0x5C, 0x5D, 0x5E, 0x5F, 0x60, 0x61, 0x62, 0x63, 
                0x64, 0x6A, 0x6B, 0x6C, 0x6D, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79,
                0x7B, 0x7C, 0x7D, 0x80, 0x81, 0x83, 0x84, 0x87, 0x8A, 0x8B, 0xA3,
                0xBA, 0xC0, 0xC1, 0xC4, 0xC7, 0xC8, 0xC9, 0xDA, 0xDB, 0xDC, 0xDD,
                0xDF, 0xE0}

--Turns graphics on or off
function ToggleGraphics(GraphicsOn)
    snes.setlayer_bg_1(GraphicsOn)
    snes.setlayer_bg_2(GraphicsOn)
    snes.setlayer_bg_3(GraphicsOn)
    snes.setlayer_bg_4(GraphicsOn)
    snes.setlayer_obj_1(GraphicsOn)
    snes.setlayer_obj_2(GraphicsOn)
    snes.setlayer_obj_3(GraphicsOn)
    snes.setlayer_obj_4(GraphicsOn)
end

BoxRadius = 6
InputSize = (BoxRadius*2+1)*(BoxRadius*2+1)

Inputs = InputSize+1
Outputs = #ButtonNames * 2

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

TimeoutConstant = 20

MaxNodes = 1000000

--sets values for where mario is in world space
--and relative to the screen
function getPositions()
	if gameinfo.getromname() == "Super Mario World (USA)" then
		marioX = memory.read_s16_le(0x94)
		marioY = memory.read_s16_le(0x96)
		
		local layer1x = memory.read_s16_le(0x1A);
		local layer1y = memory.read_s16_le(0x1C);
		
		screenX = marioX-layer1x
		screenY = marioY-layer1y
	elseif gameinfo.getromname() == "Super Mario Bros." then
		marioX = memory.readbyte(0x6D) * 0x100 + memory.readbyte(0x86)
		marioY = memory.readbyte(0x03B8)+16
	
		screenX = memory.readbyte(0x03AD)
		screenY = memory.readbyte(0x03B8)
	end
end
--Gets whatever is in that tile. Returns 1, -1, or 0
function getTile(dx, dy)
	if gameinfo.getromname() == "Super Mario World (USA)" then
		x = math.floor((marioX+dx+8)/16)
		y = math.floor((marioY+dy)/16)
		return memory.readbyte(0x1C800 + math.floor(x/0x10)*0x1B0 + y*0x10 + x%0x10)
	elseif gameinfo.getromname() == "Super Mario Bros." then
		local x = marioX + dx + 8
		local y = marioY + dy - 16
		local page = math.floor(x/256)%2

		local subx = math.floor((x%256)/16)
		local suby = math.floor((y - 32)/16)
		local addr = 0x500 + page*13*16+suby*16+subx
		
		if suby >= 13 or suby < 0 then
			return 0
		end
		
		if memory.readbyte(addr) ~= 0 then
			return 1
		else
			return 0
		end
	end
end
--Gets the 12 sprites that the SNES can hold, or 4 if we're on NES.
function getSprites()
	if gameinfo.getromname() == "Super Mario World (USA)" then
		local sprites = {}
		for slot=0,11 do
			local status = memory.readbyte(0x14C8+slot)
			if status ~= 0 then
				spritex = memory.readbyte(0xE4+slot) + memory.readbyte(0x14E0+slot)*256
				spritey = memory.readbyte(0xD8+slot) + memory.readbyte(0x14D4+slot)*256
				sprites[#sprites+1] = {["x"]=spritex, ["y"]=spritey, ["good"]=(Sprites[memory.readbyte(0x009e + slot) - 1])}
			end
		end
		return sprites
	elseif gameinfo.getromname() == "Super Mario Bros." then
		local sprites = {}
		for slot=0,4 do
			local enemy = memory.readbyte(0xF+slot)
			if enemy ~= 0 then
				local ex = memory.readbyte(0x6E + slot)*0x100 + memory.readbyte(0x87+slot)
				local ey = memory.readbyte(0xCF + slot)+24
				sprites[#sprites+1] = {["x"]=ex,["y"]=ey}
			end
		end
		
		return sprites
	end
end
--most likely used for different versions of the game, maybe had
--more memory. Don't worry about this for now.
function getExtendedSprites()
	if gameinfo.getromname() == "Super Mario World (USA)" then
		local extended = {}
		for slot=0,11 do
			local number = memory.readbyte(0x170B+slot)
			if number ~= 0 then
				spritex = memory.readbyte(0x171F+slot) + memory.readbyte(0x1733+slot)*256
				spritey = memory.readbyte(0x1715+slot) + memory.readbyte(0x1729+slot)*256
				extended[#extended+1] = {["x"]=spritex, ["y"]=spritey, ["good"] = Sprites[memory.readbyte(0x009e + slot) - 1]}
			end
		end		
		
		return extended
	elseif gameinfo.getromname() == "Super Mario Bros." then
		return {}
	end
end

--Grabs input data for the neural network.
--Currently the inputs are simply to check each tile in Super Mario World (16x16)
--and see if there is nothing (0), a normal tile (1), or a sprite (-1).
function getInputs()
	getPositions()
	
	sprites = getSprites()
	extended = getExtendedSprites()
	
	local inputs = {}
	
	for dy=-BoxRadius*16,BoxRadius*16,16 do
		for dx=-BoxRadius*16,BoxRadius*16,16 do
			inputs[#inputs+1] = 0
			
			tile = getTile(dx, dy)
			if tile == 1 and marioY+dy < 0x1B0 then
				inputs[#inputs] = 1
			end
			
			for i = 1,#sprites do
				distx = math.abs(sprites[i]["x"] - (marioX+dx))
				disty = math.abs(sprites[i]["y"] - (marioY+dy))
                if distx <= 8 and disty <= 8 then
                    if sprites[i]["good"] == 1 then
                        inputs[#inputs] = 2
                    else
                        inputs[#inputs] = -1
                    end
				end
			end

			for i = 1,#extended do
				distx = math.abs(extended[i]["x"] - (marioX+dx))
				disty = math.abs(extended[i]["y"] - (marioY+dy))
                if distx < 8 and disty < 8 then
                    if extended[i]["good"] == 1 then
                        inputs[#inputs] = 2
                    else
                        inputs[#inputs] = -1
                    end
				end
			end
		end
	end

	return inputs
end
--normal sigmoid implementation
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
	pool.innovation = Outputs
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

	genome.maxneuron = Inputs
	mutate(genome)
	
	return genome
end
--Generates empty gene
function newGene()
	local gene = {}
	gene.into = 0
	gene.out = 0
	gene.weight = 0.0
	gene.enabled = true
	gene.innovation = 0
	
	return gene
end
--Copies gene into another one, gene2, which has the same exact properties.
function copyGene(gene)
	local gene2 = newGene()
	gene2.into = gene.into
	gene2.out = gene.out
	gene2.weight = gene.weight
	gene2.enabled = gene.enabled
	gene2.innovation = gene.innovation
	
	return gene2
end
--Creates empty neuron.
function newNeuron()
	local neuron = {}
	neuron.incoming = {}
	neuron.value = 0.0
	
	return neuron
end

function generateNetwork(genome)
	local network = {}
	network.neurons = {}
	
	for i=1,Inputs do
		network.neurons[i] = newNeuron()
	end
	
	for o=1,Outputs do
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
--Returns the outputs resulting from the inputs.
--Output is formatted to get passed into SetJoystick.
--May not have to worry about that detail too much though.
function evaluateNetwork(network, inputs)
	table.insert(inputs, 1)
    --check length
	if #inputs ~= Inputs then
		console.writeline("Incorrect number of neural network inputs.")
		return {}
	end
    --set values of network
	for i=1,Inputs do
		network.neurons[i].value = inputs[i]
	end
    --Passes the values of the neurons through then network and multiplies
    --each node by the weights of it's incoming neurons.
	for _,neuron in pairs(network.neurons) do
		local sum = 0
		for j = 1,#neuron.incoming do
			local incoming = neuron.incoming[j]
			local other = network.neurons[incoming.into]
			sum = sum + incoming.weight * other.value
		end
        --sigmoiiiiiid
		if #neuron.incoming > 0 then
			neuron.value = sigmoid(sum)
		end
	end
    --stores the result of all the output nodes.
	local outputs = {}
    local powerup = memory.readbyte(0x0019)
    local offset = 0
    if powerup > 0 then
        offset = Outputs / 2
    end
	for o=1,#ButtonNames do
		local button = "P1 " .. ButtonNames[o]
		if network.neurons[MaxNodes + o + offset].value > 0 then
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
    --Add input neurons if nonInput is true
	if not nonInput then
		for i=1,Inputs do
			neurons[i] = true
		end
	end
    --Add output neurons
	for o=1,Outputs do
		neurons[MaxNodes+o] = true
	end
    --For each gene, add into nodes and out nodes.
	for i=1,#genes do
		if (not nonInput) or genes[i].into > Inputs then
			neurons[genes[i].into] = true
		end
		if (not nonInput) or genes[i].out > Inputs then
			neurons[genes[i].out] = true
		end
	end
    --r u fkn srs
    --counts number of neurons added lmao
	local count = 0
	for _,_ in pairs(neurons) do
		count = count + 1
	end
    --random index to choose from
	local n = math.random(1, count)
    --is this really be the best way jesus	
	for k,v in pairs(neurons) do
		n = n-1
		if n == 0 then
			return k
		end
	end
    --return 0 if neuron not found
	return 0
end
--Checks if link exists in genes.
function containsLink(genes, link)
	for i=1,#genes do
		local gene = genes[i]
		if gene.into == link.into and gene.out == link.out then
			return true
		end
	end
end
--mutate connection weights for genes. Why does it do it this way? idk.
function pointMutate(genome)
    --I don't even
	local step = genome.mutationRates["step"]
    --some random ass function lmao
	for i=1,#genome.genes do
		local gene = genome.genes[i]
		if math.random() < PerturbChance then
			gene.weight = gene.weight + math.random() * step*2 - step
		else
			gene.weight = math.random()*4-2
		end
	end
end
--
function linkMutate(genome, forceBias)
    --Gets 2 random neurons, one input and one nonInput.
	local neuron1 = randomNeuron(genome.genes, false)
	local neuron2 = randomNeuron(genome.genes, true)
	--Get a new gene
	local newLink = newGene()
    --Return if both nodes are input nodes.
	if neuron1 <= Inputs and neuron2 <= Inputs then
		return
	end
    --uhhhhh....lmao
    --Make sure neuron 1 is an input neuron if necessary
	if neuron2 <= Inputs then
		local temp = neuron1
		neuron1 = neuron2
		neuron2 = temp
	end
    --set into and out to links
	newLink.into = neuron1
	newLink.out = neuron2
    --idk mane
	if forceBias then
		newLink.into = Inputs
	end
    --return if this link already exists
	if containsLink(genome.genes, newLink) then
		return
	end
    --add another output
	newLink.innovation = newInnovation()
    --give this new link a random weight [-2, 2)
	newLink.weight = math.random()*4-2
    --finally, add this new link into the gene pool
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
--random enable or disable mutation for a gene
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
--Mutate various features of the genome once again based on random probability.
function mutate(genome)
    --For each mutatable attribute, ether raise the rate by 5% or lower it by
    --5%.
	for mutation,rate in pairs(genome.mutationRates) do
		if math.random(1,2) == 1 then
			genome.mutationRates[mutation] = 0.95*rate
		else
			genome.mutationRates[mutation] = 1.05263*rate
		end
	end
    --Mutate weights of genes
	if math.random() < genome.mutationRates["connections"] then
		pointMutate(genome)
	end
    --Add another link/node
	local p = genome.mutationRates["link"]
	while p > 0 do
		if math.random() < p then
			linkMutate(genome, false)
		end
		p = p - 1
	end
    --Add another link 
	p = genome.mutationRates["bias"]
	while p > 0 do
		if math.random() < p then
			linkMutate(genome, true)
		end
		p = p - 1
	end
    --Add another link with forced bias
	p = genome.mutationRates["node"]
	while p > 0 do
		if math.random() < p then
			nodeMutate(genome)
		end
		p = p - 1
	end
    --randomly enable mutatations for a genome
	p = genome.mutationRates["enable"]
	while p > 0 do
		if math.random() < p then
			enableDisableMutate(genome, true)
		end
		p = p - 1
	end
    --randomly enable mutatations for a genome
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

--Takes all genomes in the species in the current pool and
--ranks them. Stores the rankings in global.
--Does this even do anything????
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

--Selection process.
--Grabs the highest fitness genomes. If cutToOne is true,
--then it selects the highest. If not, then it selects
--the top half.
function cullSpecies(cutToOne)
	for s = 1,#pool.species do
        --sort a species genomes by fitness
		local species = pool.species[s]
		table.sort(species.genomes, function (a,b)
			return (a.fitness > b.fitness)
		end)
        --Leave half the species remaining.
		local remaining = math.ceil(#species.genomes/2)
        --if cutToOne is true, we only want one remaining.
		if cutToOne then
			remaining = 1
		end
		while #species.genomes > remaining do
			table.remove(species.genomes)
		end
	end
end
--randomly crossover child genomes with two parents. Else just use one.
--Also mutate the child.
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

--learn what staleness is??????
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

--Removes species below average fitness.
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
--Searches the species in the gene pool until we find one that matches the
--child and then adds the child to that species. If the species
--does not exist then create a new species.
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
--Begins selection process for what species to remove and creation of new
--genomes for the next generation.
function newGeneration()
    --gets top half of species
	cullSpecies(false)
	rankGlobally()
    --literally wat
	removeStaleSpecies()
	rankGlobally()
	for s = 1,#pool.species do
		local species = pool.species[s]
		calculateAverageFitness(species)
	end
	removeWeakSpecies()
	local sum = totalAverageFitness()
	local children = {}
    --for each species, bread some number of children
	for s = 1,#pool.species do
		local species = pool.species[s]
		breed = math.floor(species.averageFitness / sum * Population) - 1
		for i=1,breed do
			table.insert(children, breedChild(species))
		end
	end
    --cuts all but the top member of each species
	cullSpecies(true)
    --keep making children until we match the population
	while #children + #pool.species < Population do
		local species = pool.species[math.random(1, #pool.species)]
		table.insert(children, breedChild(species))
	end
    --for each child, add to species
	for c=1,#children do
		local child = children[c]
		addToSpecies(child)
	end
    --increment generation
	pool.generation = pool.generation + 1
    --save data
    if Testname == nil then
        writeFile("backup." .. pool.generation .. "." .. forms.gettext(saveLoadFile))
    else
        writeFile(Testname)
        writeData(Testname)
    end
end
--creates a new pool for testing and creates a population of the same species
function initializePool()
	pool = newPool()

	for i=1,Population do
		basic = basicGenome()
		addToSpecies(basic)
	end

	initializeRun()
end
--sets all inputs to false
function clearJoypad()
	controller = {}
	for b = 1,#ButtonNames do
		controller["P1 " .. ButtonNames[b]] = false
	end
	joypad.set(controller)
end
--loads the save state, clears inputs, and grabs a genome to start testing on
function initializeRun()
	savestate.load(Filename);
	rightmost = 0
	pool.currentFrame = 0
	timeout = TimeoutConstant
	clearJoypad()
	
	local species = pool.species[pool.currentSpecies]
	local genome = species.genomes[pool.currentGenome]
	generateNetwork(genome)
	evaluateCurrent()
end
--Is called every 12 frames. Gets the output of the network.
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
--r u fkn srs
if pool == nil then
	initializePool()
end
--Sets the current genome to the next one in the generation. If there are none,
--then it moves on to the next generation.
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
--If fitness > 0, we have already done the calculations for it.
function fitnessAlreadyMeasured()
	local species = pool.species[pool.currentSpecies]
	local genome = species.genomes[pool.currentGenome]

	return genome.fitness ~= 0
end
--gui stuff
function displayGenome(genome)
    --get relevant data
	local network = genome.network
	local cells = {}
	local i = 1
	local cell = {}
    --get value of the cell as well as the position, then add the cell to the
    --list of cells
	for dy=-BoxRadius,BoxRadius do
		for dx=-BoxRadius,BoxRadius do
			cell = {}
			cell.x = 50+5*dx
			cell.y = 70+5*dy
			cell.value = network.neurons[i].value
			cells[i] = cell
			i = i + 1
		end
	end
    --gets bias cell
	local biasCell = {}
	biasCell.x = 80
	biasCell.y = 110
	biasCell.value = network.neurons[Inputs].value
	cells[Inputs] = biasCell
    --For each output node, create another cell with the cell position and
    --the value. Set the color to either blue for output is activated or
    --black for output is not.
	for o = 1,Outputs do
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
        local button = ButtonNames[o]
        if o > #ButtonNames then
            button = ButtonNames[o - #ButtonNames]
        end
		gui.drawText(223, 24+8*o, button, color, 9)
	end
    --For add each neuron to the thing
	for n,neuron in pairs(network.neurons) do
		cell = {}
		if n > Inputs and n <= MaxNodes then
			cell.x = 140
			cell.y = 40
			cell.value = neuron.value
			cells[n] = cell
		end
	end
    --set up stuff for line drawing
	for n=1,4 do
		for _,gene in pairs(genome.genes) do
			if gene.enabled then
				local c1 = cells[gene.into]
				local c2 = cells[gene.out]
				if gene.into > Inputs and gene.into <= MaxNodes then
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
				if gene.out > Inputs and gene.out <= MaxNodes then
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
    --draw map box
	gui.drawBox(50-BoxRadius*5-3,70-BoxRadius*5-3,50+BoxRadius*5+2,70+BoxRadius*5+2,0xFF000000, 0x80808080)
    --draw input and output cells
	for n,cell in pairs(cells) do
		if n > Inputs or cell.value ~= 0 then
            --initialize as white and opaque
            local opacity = 0xFF000000
            local color = 0xFFFFFFFF
            --red if mushroom. not possible with current network implementation.
            --black if sprite
            --white otherwise
            if cell.value < 0 then
                color = 0xFF000000
            end
			gui.drawBox(cell.x-2,cell.y-2,cell.x+2,cell.y+2,opacity,color)
		end
	end
    --draw lines now
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

function writeFile(filename)
    local file = io.open(filename, "w")
	file:write(pool.generation .. "\n")
	file:write(pool.maxFitness .. "\n")
	file:write(#pool.species .. "\n")
    for n,species in pairs(pool.species) do
        file:write(species.topFitness .. "\n")
        file:write(species.staleness .. "\n")
        file:write(#species.genomes .. "\n")
        for m,genome in pairs(species.genomes) do
            file:write(genome.fitness .. "\n")
            file:write(genome.maxneuron .. "\n")
            for mutation,rate in pairs(genome.mutationRates) do
                file:write(mutation .. "\n")
                file:write(rate .. "\n")
            end
            file:write("done\n")
            
            file:write(#genome.genes .. "\n")
            for l,gene in pairs(genome.genes) do
                file:write(gene.into .. " ")
                file:write(gene.out .. " ")
                file:write(gene.weight .. " ")
                file:write(gene.innovation .. " ")
                if(gene.enabled) then
                    file:write("1\n")
                else
                    file:write("0\n")
                end
            end
        end
    end
    file:close()
end

function writeData(filename)
    local file = io.open(filename .. "." .. pool.generation .. ".data", "w")
    file:write(pool.generation .. "\n")
    file:write(pool.maxFitness .. "\n")
    local totalfitness = 0
    local numgenes = 0
    for n,species in pairs(pool.species) do
        for m,genome in pairs(species.genomes) do
            totalfitness = totalfitness + genome.fitness
            numgenes = numgenes + 1
        end
    end
    local avgfitness = totalfitness / numgenes
    file:write(avgfitness .. "\n")
    file:close()
end

function savePool()
	local filename = forms.gettext(saveLoadFile)
    Testname = filename
	writeFile(filename)
end

function loadFile(filename)
        local file = io.open(filename, "r")
    Testname = filename
	pool = newPool()
	pool.generation = file:read("*number")
	pool.maxFitness = file:read("*number")
	forms.settext(maxFitnessLabel, "Max Fitness: " .. math.floor(pool.maxFitness))
        local numSpecies = file:read("*number")
        for s=1,numSpecies do
		local species = newSpecies()
		table.insert(pool.species, species)
		species.topFitness = file:read("*number")
		species.staleness = file:read("*number")
		local numGenomes = file:read("*number")
		for g=1,numGenomes do
			local genome = newGenome()
			table.insert(species.genomes, genome)
			genome.fitness = file:read("*number")
			genome.maxneuron = file:read("*number")
			local line = file:read("*line")
			while line ~= "done" do
				genome.mutationRates[line] = file:read("*number")
				line = file:read("*line")
			end
			local numGenes = file:read("*number")
			for n=1,numGenes do
				local gene = newGene()
				table.insert(genome.genes, gene)
				local enabled
				gene.into, gene.out, gene.weight, gene.innovation, enabled = file:read("*number", "*number", "*number", "*number", "*number")
				if enabled == 0 then
					gene.enabled = false
				else
					gene.enabled = true
				end
				
			end
		end
	end
        file:close()
	
	while fitnessAlreadyMeasured() do
		nextGenome()
	end
	initializeRun()
	pool.currentFrame = pool.currentFrame + 1
end
 
function loadPool()
	local filename = forms.gettext(saveLoadFile)
	loadFile(filename)
end

function playTop()
	local maxfitness = 0
	local maxs, maxg
	for s,species in pairs(pool.species) do
		for g,genome in pairs(species.genomes) do
			if genome.fitness > maxfitness then
				maxfitness = genome.fitness
				maxs = s
				maxg = g
			end
		end
	end
	
	pool.currentSpecies = maxs
	pool.currentGenome = maxg
	pool.maxFitness = maxfitness
	forms.settext(maxFitnessLabel, "Max Fitness: " .. math.floor(pool.maxFitness))
	initializeRun()
	pool.currentFrame = pool.currentFrame + 1
	return
end

function onExit()
	forms.destroy(form)
end

function InitSpriteList()
    local k = 1
    for i=1, 256 do
        local isGood = (k <= #GoodSprites) and (GoodSprites[k] == i - 1)
        if isGood then
            k = k + 1
            Sprites[#Sprites + 1] = 1
        else
            Sprites[#Sprites + 1] = 0
        end
    end
end

writeFile("temp.pool")

event.onexit(onExit)

form = forms.newform(200, 260, "Fitness")
maxFitnessLabel = forms.label(form, "Max Fitness: " .. math.floor(pool.maxFitness), 5, 8)
showNetwork = forms.checkbox(form, "Show Map", 5, 30)
showMutationRates = forms.checkbox(form, "Show M-Rates", 5, 52)
restartButton = forms.button(form, "Restart", initializePool, 5, 77)
saveButton = forms.button(form, "Save", savePool, 5, 102)
loadButton = forms.button(form, "Load", loadPool, 80, 102)
saveLoadFile = forms.textbox(form, Filename .. ".pool", 170, 25, nil, 5, 148)
saveLoadLabel = forms.label(form, "Save/Load:", 5, 129)
playTopButton = forms.button(form, "Play Top", playTop, 5, 170)
hideBanner = forms.checkbox(form, "Hide Banner", 5, 190)

InitSpriteList()

while true do
	local backgroundColor = 0xD0FFFFFF
	if not forms.ischecked(hideBanner) then
		gui.drawBox(0, 0, 300, 26, backgroundColor, backgroundColor)
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

	getPositions()
	if marioX > rightmost then
		rightmost = marioX
		timeout = TimeoutConstant
	end
	--memory.writebyte(0x0019, PowerUpType_Curr)

	timeout = timeout - 1

	local timeoutBonus = pool.currentFrame / 4
	if timeout + timeoutBonus <= 0 then
		local fitness = rightmost
		if gameinfo.getromname() == "Super Mario World (USA)" and rightmost > 4816 then
			fitness = fitness + 1000
		end
		if gameinfo.getromname() == "Super Mario Bros." and rightmost > 3186 then
			fitness = fitness + 1000
		end
		if fitness == 0 then
			fitness = -1
		end
		genome.fitness = fitness
		
		if fitness > pool.maxFitness then
			pool.maxFitness = fitness
			forms.settext(maxFitnessLabel, "Max Fitness: " .. math.floor(pool.maxFitness))
			writeFile("backup." .. pool.generation .. "." .. forms.gettext(saveLoadFile))
		end
        if Testname ~= nil then
            writeData(Testname)
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
		gui.drawText(0, 0, "Gen " .. pool.generation .. " species " .. pool.currentSpecies .. " genome " .. pool.currentGenome .. " (" .. math.floor(measured/total*100) .. "%)", 0xFF000000, 11)
		gui.drawText(0, 12, "Fitness: " .. math.floor(rightmost - (pool.currentFrame) / 2 - (timeout + timeoutBonus)*2/3), 0xFF000000, 11)
		gui.drawText(100, 12, "Max Fitness: " .. math.floor(pool.maxFitness), 0xFF000000, 11)
	end
		
	pool.currentFrame = pool.currentFrame + 1

	emu.frameadvance();
end
