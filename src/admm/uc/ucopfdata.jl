using Random

mutable struct UCBus
  bus_i::Int
  bustype::Int
  # Pd::Vector{Float64}
  # Qd::Vector{Float64}
  Pd::Float64
  Qd::Float64
  Gs::Float64
  Bs::Float64
  area::Int
  Vm::Float64
  Va::Float64
  baseKV::Float64
  zone::Int
  Vmax::Float64
  Vmin::Float64
end

mutable struct UCGener
  # .gen fields
  bus::Int
  Pg::Float64
  Qg::Float64
  Qmax::Float64
  Qmin::Float64
  Vg::Float64
  mBase::Float64
  status::Int
  Pmax::Float64
  Pmin::Float64
  Pc1::Float64
  Pc2::Float64
  Qc1min::Float64
  Qc1max::Float64
  Qc2min::Float64
  Qc2max::Float64
  ramp_agc::Float64
  ru::Float64
  rd::Float64
  minUp::Int
  minDown::Int
  # .gencost fields
  gentype::Int
  startup::Float64
  shutdown::Float64
  n::Int
  coeff::Array
end

mutable struct UCGenerVec{VI, VD}
  # .gen fields
  bus::VI
  Pg::VD
  Qg::VD
  Qmax::VD
  Qmin::VD
  Vg::VD
  mBase::VD
  status::VI
  Pmax::VD
  Pmin::VD
  Pc1::VD
  Pc2::VD
  Qc1min::VD
  Qc1max::VD
  Qc2min::VD
  Qc2max::VD
  ramp_agc::VD
  ru::VD
  rd::VD
  minUp::VI
  minDown::VI
  # .gencost fields
  gentype::VI
  startup::VD
  shutdown::VD
  n::VI
#  coeff::Array
  coeff2::VD
  coeff1::VD
  coeff0::VD

  function UCGenerVec{VI, VD}(num_on) where {VI, VD}
    genvec = new{VI, VD}()
    genvec.bus = VI(undef, num_on)
    genvec.Pg = VD(undef, num_on)
    genvec.Qg = VD(undef, num_on)
    genvec.Qmax = VD(undef, num_on)
    genvec.Qmin = VD(undef, num_on)
    genvec.Vg = VD(undef, num_on)
    genvec.mBase = VD(undef, num_on)
    genvec.status = VI(undef, num_on)
    genvec.Pmax = VD(undef, num_on)
    genvec.Pmin = VD(undef, num_on)
    genvec.Pc1 = VD(undef, num_on)
    genvec.Pc2 = VD(undef, num_on)
    genvec.Qc1min = VD(undef, num_on)
    genvec.Qc1max = VD(undef, num_on)
    genvec.Qc2min = VD(undef, num_on)
    genvec.Qc2max = VD(undef, num_on)
    genvec.ramp_agc = VD(undef, num_on)
    genvec.ru = VD(undef, num_on)
    genvec.rd = VD(undef, num_on)
    genvec.minUp = VI(undef, num_on)
    genvec.minDown = VI(undef, num_on)
    genvec.gentype = VI(undef, num_on)
    genvec.startup = VD(undef, num_on)
    genvec.shutdown = VD(undef, num_on)
    genvec.n = VI(undef, num_on)
    genvec.coeff0 = VD(undef, num_on)
    genvec.coeff1 = VD(undef, num_on)
    genvec.coeff2 = VD(undef, num_on)
    return genvec
  end
end

struct UCOPFData
  buses::Array{UCBus}
  lines::Array{Line}
  generators::Array{UCGener}
  bus_ref::Int
  baseMVA::Float64
  BusIdx::Dict{Int,Int}    #map from bus ID to bus index
  FromLines::Array         #From lines for each bus (Array of Array)
  ToLines::Array           #To lines for each bus (Array of Array)
  BusGenerators::Array     #list of generators for each bus (Array of Array)
end

function uc_initialize(gen::UCGener)::Nothing
  typenum = rand(1:3)
  if typenum == 1
    gen.ru        = max(gen.Pmin, gen.Pmax / 2)
    gen.rd        = max(gen.Pmin, gen.Pmax / 2)
    gen.minUp     = 2
    gen.minDown   = 2
  elseif typenum == 2
    gen.ru        = max(gen.Pmin, gen.Pmax / 3)
    gen.rd        = max(gen.Pmin, gen.Pmax / 3)
    gen.minUp     = 3
    gen.minDown   = 3
  else
    gen.ru        = max(gen.Pmin, gen.Pmax / 5)
    gen.rd        = max(gen.Pmin, gen.Pmax / 5)
    gen.minUp     = 4
    gen.minDown   = 4
  end
  return
end

function ucopf_loaddata(case_name, lineOff=Line())
  Random.seed!(35)
  #
  # load buses
  #
  # bus_arr = readdlm("data/" * case_name * ".bus")
  bus_arr = readdlm(case_name * ".bus")
  num_buses = size(bus_arr,1)
  buses = Array{UCBus}(undef, num_buses)
  bus_ref=-1
  for i in 1:num_buses
    @assert bus_arr[i,1]>0  #don't support nonpositive bus ids
    buses[i] = UCBus(bus_arr[i,1:13]...) # TODO: CHANGE THIS SO THAT PD AND QD VECTORS ARE FED INTO INITIALIZATION
    buses[i].Va *= pi/180
    if buses[i].bustype==3
      if bus_ref>0
        error("More than one reference bus present in the data")
      else
         bus_ref=i
      end
    end
    #println("bus ", i, " ", buses[i].Vmin, "      ", buses[i].Vmax)
  end

  #
  # load branches/lines
  #
  # branch_arr = readdlm("data/" * case_name * ".branch")
  branch_arr = readdlm(case_name * ".branch")
  num_lines = size(branch_arr,1)
  lines_on = findall((branch_arr[:,11].>0) .& ((branch_arr[:,1].!=lineOff.from) .| (branch_arr[:,2].!=lineOff.to)) )
  num_on   = length(lines_on)

  if lineOff.from>0 && lineOff.to>0
    println("opf_loaddata: was asked to remove line from,to=", lineOff.from, ",", lineOff.to)
    #println(lines_on, branch_arr[:,1].!=lineOff.from, branch_arr[:,2].!=lineOff.to)
  end
  if length(findall(branch_arr[:,11].==0))>0
    println("opf_loaddata: ", num_lines-length(findall(branch_arr[:,11].>0)), " lines are off and will be discarded (out of ", num_lines, ")")
  end



  lines = Array{Line}(undef, num_on)

  lit=0
  for i in lines_on
    @assert branch_arr[i,11] == 1  #should be on since we discarded all other
    lit += 1
    lines[lit] = Line(branch_arr[i, 1:13]...)
    #=
    if (lines[lit].angmin != 0 || lines[lit].angmax != 0) && (lines[lit].angmin>-360 || lines[lit].angmax<360)
      println("Voltage bounds on line ", i, " with angmin ", lines[lit].angmin, " and angmax ", lines[lit].angmax)
      error("Bounds of voltage angles are still to be implemented.")
    end
    =#

  end
  @assert lit == num_on

  #
  # load generators
  #
  # gen_arr = readdlm("data/" * case_name * ".gen")
  gen_arr = readdlm(case_name * ".gen")
  # costgen_arr = readdlm("data/" * case_name * ".gencost")
  costgen_arr = readdlm(case_name * ".gencost")
  num_gens = size(gen_arr,1)

  baseMVA=100

  @assert num_gens == size(costgen_arr,1)

  gens_on=findall(x->x!=0, gen_arr[:,8]); num_on=length(gens_on)
  if num_gens-num_on>0
    println("loaddata: ", num_gens-num_on, " generators are off and will be discarded (out of ", num_gens, ")")
  end

  generators = Array{UCGener}(undef, num_on)
  i=0
  for git in gens_on
    i += 1
    generators[i] = UCGener(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, Array{Int}(undef, 0)) #gen_arr[i,1:end]...)

    generators[i].bus      = gen_arr[git,1]
    generators[i].Pg       = gen_arr[git,2] / baseMVA
    generators[i].Qg       = gen_arr[git,3] / baseMVA
    generators[i].Qmax     = gen_arr[git,4] / baseMVA
    generators[i].Qmin     = gen_arr[git,5] / baseMVA
    generators[i].Vg       = gen_arr[git,6]
    generators[i].mBase    = gen_arr[git,7]
    generators[i].status   = gen_arr[git,8]
    @assert generators[i].status==1
    generators[i].Pmax     = gen_arr[git,9]  / baseMVA
    generators[i].Pmin     = gen_arr[git,10] / baseMVA
    generators[i].Pc1      = gen_arr[git,11]
    generators[i].Pc2      = gen_arr[git,12]
    generators[i].Qc1min   = gen_arr[git,13]
    generators[i].Qc1max   = gen_arr[git,14]
    generators[i].Qc2min   = gen_arr[git,15]
    generators[i].Qc2max   = gen_arr[git,16]
    uc_initialize(generators[i])
    generators[i].gentype  = costgen_arr[git,1]
    generators[i].startup  = costgen_arr[git,2]
    generators[i].shutdown = costgen_arr[git,3]
    generators[i].n        = costgen_arr[git,4]
    @assert(generators[i].n <= 3 && generators[i].n >= 2)
    if generators[i].gentype == 1
      generators[i].coeff = costgen_arr[git,5:end]
      error("Piecewise linear costs remains to be implemented.")
    else
      if generators[i].gentype == 2
        generators[i].coeff = costgen_arr[git,5:end]
        #println(generators[i].coeff, " ", length(generators[i].coeff), " ", generators[i].coeff[2])
      else
        error("Invalid generator cost model in the data.")
      end
    end
  end

  # build a dictionary between buses ids and their indexes
  busIdx = mapBusIdToIdx(buses)

  # set up the FromLines and ToLines for each bus
  FromLines,ToLines = mapLinesToBuses(buses, lines, busIdx)

  # generators at each bus
  BusGeners = mapGenersToBuses(buses, generators, busIdx)

  #println(generators)
  #println(bus_ref)
  return UCOPFData(buses, lines, generators, bus_ref, baseMVA, busIdx, FromLines, ToLines, BusGeners)
end