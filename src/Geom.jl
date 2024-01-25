function abs_angle_cartesian(relD, batiD)
    ###########################################################################    
    # Absolute angle in cartesian notation, angle between [180,-180], 
    # 0 is in EAST & positive counterclockwise.
    # From a relative angle from wave & bathymetry.
    # The same as rel_angle_cartesian[relD,-1*batiD]
    # INPUT:
    # relD:     relative wave angle between wave and bathymetry; 0 is the bathymetry & positive counterclockwise.
    # batiD:    bathymetry angle (normal to the shoreline) in Cartesian notation.
    #
    # OUTPUT:
    # waveD:    wave angle in Cartesian notation.
    ###########################################################################    
    
    waveD = relD + batiD
    replace!(waveD,)
    waveD[waveD.>180]=waveD[waveD.>180].-360
    waveD[waveD.<-180]=waveD[waveD.<-180].+360
        
    return waveD
end

function abs_pos(X0,Y0,phi,dn)
    #####################    
    # INPUT:
    #
    # X0 : x coordinate; origin of the transect
    # Y0 : y coordinate; origin of the transect
    # phi : transect orientation in radians
    # dn : position on the transect
    #
    # OUTPUT:
    # XN : x coordinate
    # YN : y coordinate
    #####################

    XN = X0 .+ dn.*cos.(phi)
    YN = Y0 .+ dn.*sin.(phi)
       
    return XN, YN
end

function cartesianDir2nauticalDir(cDir)
    ###########################################################################    
    # Cartesian convention with 0 in East & positive counterclockwise TO
    # Nautical convention with 0 in North & positive clockwise. 
    ###########################################################################    
    
    nDir = 90.0 .- cDir
    nDir[nDir.<0] = 360.0 .+ nDir[nDir.<0]
    
    return nDir
end

function interp_lon(x,lon,xq,varargin)

    # INTERP_LON interpolates a set of longitude angles [in deg]
    #
    # Usage: out = interp_lon(x,lon,xq)
    #
    # x & lon are vectors of length N.  function evalutes longitude 
    # (in deg -180..180) at points xq using unwrap & interp1()
    #
    # to specify interpolation method used in interp1; use
    # out = interp_lon(x,lon,xq,METHOD)
    #
    # Written by D.G. Long; 27 Nov 2017 
    
    ulon=unwrap(lon*pi/180)*180/pi
    if nargin>3
      out=interp1(x,ulon,xq,varargin[1])
    else()
      out=interp1(x,ulon,xq)
    end
    out=out%360
    out[out.>180]=out[out.>180]-360
    return out
end

function nauticalDir2cartesianDir(nDir)

    ###########################################################################    
    # Nautical convention with 0 in North & positive clockwise TO 
    # Cartesian convention with 0 in East & positive counterclockwise.
    ###########################################################################    
   
    cDir = 90.0 .- nDir
    for i in eachindex(cDir)
        if cDir[i] < -180.0
            cDir[i] = 360.0 .+ cDir[i]
        end
    end
   
    
    return cDir
end

function nauticalDir2cartesianDirP(nDir)

    ###########################################################################    
    # Nautical convention with 0 in North & positive clockwise TO 
    # Cartesian convention with 0 in East & positive counterclockwise.
    ###########################################################################    

    cDir = 90.0 - nDir
    
    if cDir < -180.0
        cDir = cDir + 360.0
    end
 
    
    return cDir
end

function pol2cart(rho,phi)

    x = rho .* cos(phi)
    y = rho .* sin(phi)
    return x,y
end

function rel_angle_cartesian(waveD, batiD)

    ###########################################################################    
    # Relative angle (in degrees) between wave direction & bathymetry with 
    # angles in cartesian coordinates, angle between [180,-180], 
    # 0 is in EAST & positive counterclockwise.
    #
    # INPUT:
    # waveD:    wave angle in Cartesian notation.
    # batiD:    bathymetry angle (normal to the shoreline) in Cartesian notation.
    #
    # OUTPUT:
    # relD:     relative wave angle between wave and bathymetry; 0 is the bathymetry & positive counterclockwise.
    ###########################################################################    
    
    
    relD = waveD .- batiD
    relD[relD .> 180.] .= relD[relD .> 180.] .- 360.0
    relD[relD .< -180.] .= relD[relD .< -180.] .+ 360.0
    
    return relD
end

function rel_angle_cartesianP(waveD, batiD)

    ###########################################################################    
    # Relative angle (in degrees) between wave direction & bathymetry with 
    # angles in cartesian coordinates, angle between [180,-180], 
    # 0 is in EAST & positive counterclockwise.
    #
    # INPUT:
    # waveD:    wave angle in Cartesian notation.
    # batiD:    bathymetry angle (normal to the shoreline) in Cartesian notation.
    #
    # OUTPUT:
    # relD:     relative wave angle between wave and bathymetry; 0 is the bathymetry & positive counterclockwise.
    ###########################################################################    

    
    relD = waveD - batiD
    if relD > 180
        relD = relD - 360.0
    elseif relD < -180.0
        relD = relD + 360.0
    end
 
    
    return relD
end

function shore_angle(XN,YN,wave_angle)
    #####################    
    # INPUT:
    #
    # XN : x coordinate
    # YN : y coordinate
    # wave_angle : wave angle()
    #
    # OUTPUT:
    # shoreAng : angle of the shoreline to compute sediment transport
    #####################
    
    shoreAng = zeros(1,length(XN)-1)
    # method = zeros(length(XN)-1)
    
    for i in eachindex(shoreAng)
        alfa = atan(YN[i+1]-YN[i],XN[i+1]-XN[i]).*180.0./pi
        relang = rel_angle_cartesianP(nauticalDir2cartesianDirP(wave_angle[i+1]),alfa)
    #         # check the relative angle between the waves & the shoreline orientation [not the normal to the shoreline]
        if abs(relang).>=45 && abs(relang).<=180-45
            shoreAng[i] = alfa
    #         method[i] = "CenDiff"
        elseif abs(relang).<45  && i.<length(shoreAng)-1
            try
                shoreAng[i] = atan(YN[i+2]-YN[i+1],XN[i+2]-XN[i+1]).*180.0./pi
    #             method[i] = "UpWind"
            catch
                shoreAng[i] = alfa
    #             method[i] = "CenDiff"
            end
                    
        elseif abs(relang).>180 - 45 && i.>1
            try
                shoreAng[i] = atan(YN[i]-YN[i-1],XN[i]-XN[i-1]).*180.0./pi
    #             method[i] = "UpWind"
            catch
                shoreAng[i] = alfa
    #             method[i] = "CenDiff"
            end
        else
                shoreAng[i] = alfa
    #             method[i] = "CenDiff"
        end
    end
                
    return shoreAng
    
end
