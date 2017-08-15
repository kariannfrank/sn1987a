pro makegrid

common grid_com, xpos, ypos

xpos = (float(indgen(200)) - 100.0) - 0.5
xpos = rebin(xpos,200,200,/sample)
ypos = transpose(xpos)

end
