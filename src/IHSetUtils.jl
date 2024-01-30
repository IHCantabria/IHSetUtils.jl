module IHSetUtils

using Polynomials
# using Optim
using SciPy
# using NLsolve
export abs_angle_cartesian, abs_pos, cartesianDir2nauticalDir, interp_lon, nauticalDir2cartesianDir, nauticalDir2cartesianDirP, pol2cart, rel_angle_cartesian, rel_angle_cartesianP, shore_angle
export ADEAN, ALST, Bruun_rule, depthOfClosure, Hs12Calc, wast, wMOORE, deanSlope
export BreakingPropagation, GroupCelerity, hunt, LinearShoal, LinearShoalBreak_Residual, LinearShoalBreak_ResidualVOL, RelDisp, RU2_Stockdon2006, Snell_Law

include("Geom.jl")
include("Morfo.jl")
include("Waves.jl")

end
