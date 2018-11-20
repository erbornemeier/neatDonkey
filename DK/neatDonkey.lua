
--memory.writebyte(0x043E, 2) sweet jumps
--memory.writebyte(0x0055,3) --change lives

--------------------------------------
--drawing constants ------------------
--------------------------------------
viewDist = 8
scale = 8
inputDisplayOffset = 10
--------------------------------------

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
--------------------------------------
--emulator controls and access--------
--------------------------------------
function nextFrame()
	framesElapsed = framesElapsed + 1
	emu.frameadvance()
end

function getTimeElapsed()
	return framesElapsed / 60.0
end

function outputsToController(outputs)
	controller = {}
	for button, action in pairs(outputs) do
		controller["P1 " .. button] = action
	end
	joypad.set(controller)
end

function getMarioLocation()
	local marioX = memory.readbyte(0x0046)+2
	local marioY = memory.readbyte(0x0047)
	--gui.drawPixel(marioX, marioY, 0xFF00FF00)
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
		sprites[#sprites+1] = {['x']=x+2, ['y']=y-4}
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
	inputs = {}
	for i=0,viewDist*2 do
		inputs[i] = {}
		for j=0, viewDist*2 do
			inputs[i][j] = 0
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
		   inputs[ix][iy] = 1
		end
	end
	sprites = getSpriteLocations()
	for _, sprite in pairs(sprites) do
		x = sprite['x']
		y = sprite['y']
		x_diff = x-mario['x']
		y_diff = y-mario['y']
		if math.abs(x_diff) <= viewDist*0x8 and
		   math.abs(y_diff) <= viewDist*0x8 then
		   ix = round(x_diff/0x8) + viewDist
		   iy = round(y_diff/0x8) + viewDist
		   inputs[ix][iy] = -1
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
	for i=0x80E3, 0x810A, 0x03 do
		x = memory.readbyte(i)
		y = memory.readbyte(i+1)
		offset = memory.readbyte(i+2)
		ladderSizeX = memory.readbyte(0x810B + 2 + offset)
		ladderSizeY = memory.readbyte(0x810B + 3 + offset)
		--ladder 8 is funky for some reason
		if #ladders == 8 then
			y = y + 0xA
			ladderSizeY = ladderSizeY - 0xB
		end
		ladders[#ladders + 1] = {['x']=x,['y']=y,['sizeX']=ladderSizeX, ['sizeY']=ladderSizeY}
		
	end
	return ladders
end

function drawLadders(ladders)
	for _,l in pairs(ladders) do
		gui.drawBox(l['x'], l['y'] ,l['x']+l['sizeX'], l['y']+l['sizeY'], 0xFF00FF00, 0xFF00FF00)
	end
end

--************************************
--**---main state execution --------**
--************************************
saveStateFile = "DP1.state"
savestate.load(saveStateFile);
controllerButtons = {"A", "LEFT", "RIGHT", "UP", "DOWN"}

form = forms.newform(300,400, "My Form")
DEBUG = forms.checkbox(form, "Debug", 5, 5)
outputs = {["A"]=false, 
		   ["Left"]=false,
		   ["Right"]=true,
		   ["Up"]=false,
		   ["Down"]=false}
		   
-----GAME LOOP --------
while true do

	--mario = getMarioLocation()
	--sprites = getSpriteLocations()
	--platforms = getPlatformLocations()
	--drawPlatforms(platforms)
	--ladders = getLadderLocations()
	
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