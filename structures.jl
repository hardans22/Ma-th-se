
mutable struct NodeSets 
    dummy_nodes::Vector{String}
    fl_nodes::Vector{String}                   #Flight nodes
    ac_nodes::Vector{String}                   #Aircraft nodes
    mtn_nodes::Vector{String}                  #Maintenance nodes   
    mtn_nodes_stn::Dict{String,Vector{String}}              #Maintenance nodes per station
    H_aircraft::Vector{String}                 #High time aircraft nodes
    NH_aircraft::Vector{String}                #Non High Time aircraft nodes
end

mutable struct ArcSets 
    arcs_S::Vector{Tuple{String,String}}        #s-to-aircraft arcs
    arcs_K::Vector{Tuple{String,String}}        #aircraft-to-flight arcs
    arcs_F::Vector{Tuple{String,String}}        #flight-to-flight arcs
    arcs_T::Vector{Tuple{String,String}}        #flight-to-t arcs
    arcs_M::Vector{Tuple{String,String}}        #maintenance arcs
    arcs_M_bar::Vector{Tuple{String,String}}    #non maintenance arcs 
end

mutable struct ArcInfo
    arc::Tuple{String,String}
    type::String
end

mutable struct Graph
    nodes ::Vector{String}                      #All graph nodes
    arcs::Vector{Tuple{String,String}}          #All graph arcs
    node_sets ::NodeSets                        #Subsets of nodes
    arc_sets ::ArcSets                          #Subsets of arcs   
    arc_info:: Dict{Tuple{String,String},String}                         #Information on each arc 
end

mutable struct FlightData
    fl_legs::Dict{String,Vector{String}}    #Flight legs (origin, destination)
    d::Dict{String,Int64}                       #Flight durations   
    b::Dict{Tuple{String,String},Int64}
    tk::Dict{String,Int64}                      #Flight takeoff
    DT::Dict{String,Int64}                      #Departure times
    AT::Dict{String,Int64}                      #Arrival times
    init_airport::Dict{String,String}           #Initial airport of each aircraft
    init_flying_time::Dict{String,Int64}     #Initial flying time of each aircraft
    init_takeoff::Dict{String,Int64}         #Initial takeoff time of each aircraft
    init_flying_day::Dict{String,Int64}     #Initial flying days of each aircraft
end



mutable struct InstanceData
    graph::Graph
    fl_data::FlightData
    other_data::Dict{String,Any}
    TRT::Int64                                  #Turnaround times
    mtn_time::Int64                             #Maintenance time
    max_flying_time::Int64                      #Maximum flying time
    max_takeoff::Int64                          #Maximum takeoff time
    max_flying_day::Int64                      #Maximum flying days
    nbr_ac::Int64                               #Number of aircrafts
end


mutable struct Solution
    obj::Float64
    x::Dict{Tuple{String,String,String}, Float64}
    y::Dict{Tuple{String,String},Float64}
    u::Dict{String,Float64}
    v::Dict{String,Float64}
    w::Dict{String,Float64}
    rho::Dict{String,Float64}
    lambda::Dict{String,Float64}
    phi::Dict{String,Float64}
    other_info::Dict{String,Any}   #Any other information to store
    
end

mutable struct BendersSolution
    obj::Float64
    time::Float64
    mp_time::Float64
    sp_time::Float64
    prep1_time::Float64
    prep2_time::Float64
    nbr_feasibility_cuts::Int64
    nbr_optimality_cuts::Int64
    nbr_cover_cuts::Int64
    nbr_cuts::Int64
    nbr_iterations::Int64
    x_vals::Dict{Tuple{String,String,String},Float64}
    y_vals::Dict{Tuple{String,String},Float64}
    
    BendersSolution() = new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 
                           Dict{Tuple{String,String}, Float64}(), 
                           Dict{Tuple{String,String}, Float64}())
end