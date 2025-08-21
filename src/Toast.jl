module Toast

using Hexagons, Delica, Tessen, Unitful, Statistics, Meshes
using LinearAlgebra, Ahmes
import Base:convert
import Tessen:HatchLine, pointalong, intersections

export createconfig, scaffold, repjob, psweep, arrangescaffolds

"""
```julia
createconfig([filename])
```
Write an example config file to `filename`. If `filename` is omitted, the file
will be written to `"config.jl"`.

# Configuration Parameters
- dfield: calibrated FOV for the objective being used
- hbottom: distance from the bottom of the post to the bottom of the beams
- htop: distance from the bottom of the beams to the top of the bumpers
- chamfertop: angle at which the tops of the posts and beams should be chamfered
- cutangle: angles at which blocks should be cut to avoid shadowing
- overlap: amount that neighboring blocks should overlap to ensure they are connected
- nhammock: number of hammock slices to write
- dhammockslice: distance between hammock slices
- dslice: slicing distance
- wbumper: width of the bumpers
- fillet: fillet radius on xz and yz crossections of the posts and bumpers
- wpost: post width
- hbeam: beam height
- lbeammax: maximum total beam length
- maxseglength: maximum length of a single beam segment
- keygap: closest point between the two halves of a beam before they are closed with a keystone
- dhammockhatch: hammock hatching distance
- dhatch: hatching distance for posts, beams and bumpers
- laserpower: laser power for posts, beams and bumpers
- scanspeed: scan speed for posts, beams and bumpers
- stagespeed: max stage speed
- interfacepos: how far into the substrate writing should begin
- hamscanspeed: scan speed for hammocks
- hamlaserpower: laser power for hammocks
- w: Relative width of the scaffold
- h1: Relative height of the southern portion of the scaffold
- h2: Relative height of the middle portion of the scaffold
- h3: Relative height of the northern portion of the scaffold
- hexsize: Size of component hexagons
- bumperfillet: Amount to round the corners of the bumper
- bumperseglength: maximum length (along path) of a single bumper segment
- cornerdiff: How much bigger to make the bottom right corner
- hamoffset: how far 'above' the beams to place the hammocks
"""
function createconfig(filename="config.jl")
    config = """Dict(
        :dfield => 1750u"µm",
        :hbottom => 30u"µm",
        :htop => 70u"µm",
        :chamfertop => pi/6,
        :cutangle => pi/6,
        :overlap => 10u"µm",
        :dslice => 1u"µm",
        :nhammock => 5,
        :dhammockslice => 280u"nm",
        :wbumper => 200u"µm",
        :fillet => 20u"µm",
        :wpost => 50u"µm",
        :hbeam => 10u"µm",
        :lbeammax => 150u"µm",
        :maxseglength => 30u"µm",
        :keygap => 20u"µm",
        :dhammockhatch => 1u"µm",
        :dhatch => 300u"nm",
        :laserpower => 50u"mW",
        :scanspeed => 50_000u"µm/s",
        :stagespeed => 100u"µm/s",
        :interfacepos => 10u"µm",
        :hamscanspeed => 2_000u"µm/s",
        :hamlaserpower => 50u"mW",
        :w => 2,
        :h1 => 3,
        :h2 => 5,
        :h3 => 3,
        :hexsize => 50u"µm",
        :bumperfillet => 250u"µm",
        :bumperseglength => 500u"µm",
        :cornerdiff => 1,
        :hamoffset => 0u"µm"
    )
    """
    open(filename,"w") do io
        print(io,config)
    end
end

function bumperedgecoords(;kwargs...)
    #copy out the variables we want
    hbottom = kwargs[:hbottom]
    htop = kwargs[:htop]
    chamfertop = kwargs[:chamfertop]
    #I removed this parameter
    chamferbottom = 0
    wbumper = kwargs[:wbumper]
    fillet = kwargs[:fillet]
    #figuring out how to fillet everything in crossection seems hard enough that I'm just going
    #to use Tessen hatching machinery to do it for me.
    #this function will return a closure which gives the bottom and top coordinates of a bumper
    #as a function of z position

    #draw our crossection on the xz plane
    verts=[
        [-wbumper/2,hbottom],
        [-wbumper/2 + htop*tan(chamfertop),hbottom+htop],
        [wbumper/2 - htop*tan(chamfertop),hbottom+htop],
        [wbumper/2,hbottom],
        [wbumper/2 - hbottom*tan(chamferbottom),0u"µm"],
        [-wbumper/2 + hbottom*tan(chamferbottom),0u"µm"]
    ]
    nofillet = polycontour(verts)
    filleted = polycontour(verts,fillet)
    function(zcoord)
        #I'd still like to fillet the top and bottom vertices on the 'inside'
        filletinside = (zcoord < fillet) || (zcoord > (htop + hbottom - fillet))
        #place a hatchline off of the contour to the left and at our desired height
        hl=HatchLine([-wbumper,zcoord],[1,0]u"µm")
        #our first coordinate is the first intersection of hl with nofillet near the beam
        #attachments, filleted otherwise
        firstpara = filletinside ?
            sort(intersections(filleted,hl))[1] : sort(intersections(nofillet,hl))[1]
        #our second coordinate is the second intersection of hl with filleted
        secondpara = sort(intersections(filleted,hl))[2]
        map([firstpara,secondpara]) do p
            #convert to real world coordinates in Âµm
            coords = pointalong(hl,p)
            #we're interested in the x coordinate, add the dimensions on
            coords[1] * u"µm"
        end
    end
end


function outlineverts(w,h1,h2,h3)
    #assuming the origin is at the center of the bottom edge
    #and the circle which surrounds one unit of the hex grid
    #has a radius of 1 (side length is also 1)

    @assert iseven(w)
    corners = [[1+3w/2,0]]

    #the next corner is 3h1 units away at an angle of 60deg from the x axis
    push!(corners,corners[end] + 3h1*[cosd(60),sind(60)])

    #the next corner is sqrt(3)*h2 away, straight up
    push!(corners,corners[end] + [0,sqrt(3)*h2])

    #now we move 3h3 units at an angle of 120deg from the x axis
    push!(corners,corners[end] + 3h3*[cosd(120),sind(120)])

    #now finish the contour by mirroring the points
    vcat(corners,
         [[-1,1] .* c for c in reverse(corners)])
end

#ok, we're going to make a type for contours defined by offsetting edges
#we will make them convertable to Contours and define centroid and area
#functions for them
abstract type OffsetContour end

"""
```julia
centroid(oc)
```
Get the centroid of an `OffsetContour`
"""
function centroid end

"""
```julia
area(oc)
```
Get the are of an `OffsetContour`
"""
function area end

"""
```julia
getoffsetedges(oc)
```
Get the edges offset from the path forming the `middle` of a `OffsetContour`
"""

"""
```julia
LineOffsetContour(le,x1,x2)
```
Generate a closed contour by `offset`ing a `LineEdge` twice
"""
struct LineOffsetContour <: OffsetContour
    edge::LineEdge
    x1::Quantity
    x2::Quantity
end

OffsetContour(le::LineEdge,x1,x2) = LineOffsetContour(le,x1,x2)

function centroid(loc::LineOffsetContour)
    #This is the midpoint of an offset line halfway between x1 and x2
    meanoffset = Tessen.offset(loc.edge,mean([loc.x1,loc.x2]))
    #get the midpoint
    pointalong(meanoffset,span(meanoffset)/2) * 1u"μm"
end

function area(loc::LineOffsetContour)
    #this is the length of the line times the difference between the offsets
    span(loc.edge)*1u"μm" * abs(loc.x1-loc.x2)
end

function getoffsetedges(loc::LineOffsetContour)
    [Tessen.offset(loc.edge,x) for x in [loc.x1,loc.x2]]
end

struct ArcOffsetContour <: OffsetContour
    edge::ArcEdge
    x1::Quantity
    x2::Quantity
end

OffsetContour(ae::ArcEdge,x1,x2) = ArcOffsetContour(ae,x1,x2)

function getoffsetedges(aoc::ArcOffsetContour)
    [Tessen.offset(aoc.edge,x) for x in [aoc.x1,aoc.x2]]
end

function centroid(aoc::ArcOffsetContour)
    #total arc
    alpha = aoc.edge.stopangle - aoc.edge.startangle
    #make sure this is positive
    alpha = alpha < 0 ? alpha + 2pi : alpha

    #get the edges of the contour
    edges = getoffsetedges(aoc)
    r = sort([1u"μm"*e.r for e in edges])

    #get the distance between the center of curvature and the centroid
    d = (2/3)*(sin(alpha)/alpha)*(r[2]^3 - r[1]^3)/(r[2]^2 - r[1]^2)

    #now get the direction
    theta1 = aoc.edge.startangle
    theta2 = aoc.edge.stopangle
    #theta2 needs to be greater than theta1
    theta2 = theta1 > theta2 ? theta2 + 2pi : theta2

    theta = mean([theta1,theta2])
    #this is the distance from the center of curvature
    distance_from_center = [d*cos(theta), d*sin(theta)]
    1u"μm"*aoc.edge.c + distance_from_center
end

function area(aoc::ArcOffsetContour)
    #total arc
    alpha = aoc.edge.stopangle - aoc.edge.startangle
    #make sure this is positive
    alpha = alpha < 0 ? alpha + 2pi : alpha

    #get the edges of the contour
    edges = getoffsetedges(aoc)
    r = sort([1u"μm"*e.r for e in edges])
    alpha*(r[2]^2 - r[1]^2)
end

#startcap and endcap control whether the contour is closed at the beginning and end
#of an offsetcontour
function convert(::Type{Contour},oc::OffsetContour;startcap=true,endcap=true)
    offsetedges = getoffsetedges(oc)
    #start with the 'first' edge
    edges = Tessen.Edge[]
    push!(edges,offsetedges[1])
    #add on the endcap
    if endcap
        #the endcap goes from the end of the first edge to the end of the second edge
        push!(edges,LineEdge(
            1u"μm"*pointalong(offsetedges[1],span(offsetedges[1])),
            1u"μm"*pointalong(offsetedges[2],span(offsetedges[2])),
        ))
    end
    #add the 'second' edge
    push!(edges,offsetedges[2])
    #add on the startcap
    if startcap
        #the startcap goes from the beginning of the second edge to the beginning of the first edge
        push!(edges,LineEdge(
            1u"μm"*pointalong(offsetedges[2],0),
            1u"μm"*pointalong(offsetedges[1],0),
        ))
    end
    Contour(edges)
end

struct CompoundOffsetContour <: OffsetContour
    edges :: Vector{<:Tessen.Edge}
    x1::Quantity
    x2::Quantity
end

OffsetContour(edges::Vector{<:Tessen.Edge},x1,x2) = CompoundOffsetContour(edges,x1,x2)

contours(coc::CompoundOffsetContour) = OffsetContour.(coc.edges,coc.x1,coc.x2)

function area(coc::CompoundOffsetContour)
    sum(contours(coc)) do c
        area(c)
    end
end

function centroid(coc::CompoundOffsetContour)
    #this is the average of the component contour centroids weighted by area
    sum(contours(coc)) do c
        centroid(c)*area(c)
    end/area(coc)
end

function convert(::Type{Contour}, coc::CompoundOffsetContour)
    #fill me out
    c = contours(coc)
    firstedges = convert(Contour,c[1],endcap=false).edges
    lastedges = convert(Contour,c[end],startcap=false).edges
    otheredges::Vector{Vector{Tessen.Edge}}= map(convert.(Contour,c[2:end-1],startcap=false,endcap=false)) do c
        c.edges
    end
    voe::Vector{Tessen.Edge} = vcat(otheredges...)
    Contour(vcat(firstedges,voe,lastedges))
end

#add a span method for vectors of edges
function Tessen.span(edges::Vector{<:Tessen.Edge})
    sum(edges) do e
        span(e)
    end
end
#add a subsection method for vectors of edges. We will assume this corresponds to closed paths
function Tessen.subsection(edges::Vector{<:Tessen.Edge},x1::Real,x2::Real)
    #because this is a closed path we will accept x1 and x2 outside of the range 0 < x < span(edges)
    @assert abs(x2 - x1) < span(edges) "can't make a subsection longer than the contour"
    #if x1 is less than 0, figure out a corresponding positive position
    while x1 < 0
        x1 += span(edges)
        #keep x1 and x2 in sync
        x2 += span(edges)
    end
    #also make sure x1 is less than span(edges)
    while x1 > span(edges)
        x1 -= span(edges)
        x2 -= span(edges)
    end
    #now if x2 > span(edges) we need to subsection twice to stitch over the origin
    if x2 > span(edges)
        #how much to we go over the end
        overhang = x2 - span(edges)
        return vcat(subsection(edges,x1,span(edges)),subsection(edges,0,overhang))
    end
    #otherwise we're within the bounds 0 < x < span(edges)
    @assert all([x1,x2] .<= span(edges)) "Both x1 and x2 must be less than total span"
    #make a local copy of edges so we don't fuck with our caller
    edges = copy(edges)
    #find our starting place
    distuntilstart = x1
    distuntilend = x2
    while true
        thisspan = span(edges[1])
        if thisspan < distuntilstart #start point not on this edge
            popfirst!(edges)
            distuntilstart -= thisspan
            distuntilend -= thisspan
        else
            break
        end
    end
    #ok, now we know that a portion of edges[1] is on our output
    #trim off the portion of edges[1] that we know we don't want
    edges[1] = subsection(edges[1],distuntilstart,span(edges[1]))
    distuntilend -= distuntilstart
    #now we just add sections to the output until we have enough distance
    output = Tessen.Edge[]
    for e in edges
        #does this edge get us all the way to the end?
        if span(e) > distuntilend
            #trim off the end, add it to output and break
            push!(output,subsection(e,0,distuntilend))
            break
        else
            #otherwise, add it on and update distuntilend
            push!(output,e)
            distuntilend -= span(e)
        end
    end
    return output
end


#==============We don't actually need this===========================
"""
```julia
rotatestart(edges,x)
```
Modify a `Vector{Edge}` to represent the same geometry, but with the
startpoint of the (assumed closed) contour slid along the path by `x`
"""
function rotatestart(edges::Vector{Tessen.Edge},x::Quantity)
    #copy input so we don't mess up our caller
    edges = copy(edges)
    disttomove = ustrip(u"μm",x)
    #copy sections to the end until we hit disttomove
    while true
        #can we get through disttomove with the first segment in edges
        if disttomove < span(edges[1])
            #cut off disttomove from edges[1] and move to the back
            newend = subsection(edges[1],0,disttomove)
            push!(edges,newend)
            #modify what's left on the front
            edges[1] = subsection(edges[1],disttomove,span(edges[1]))
            #done
            break
        else
            #move the whole segment to the back
            newend = popfirst!(edges)
            push!(edges,newend)
            #update how far we have to go
            disttomove -= span(newend)
        end
    end
    return edges
end
=====================================================================#

function bumper(path,profile,maxseglength,zstart,zend,dslice,overlap,cutangle)
    #how much space per segment is taken up by the chamfered overlap zones
    cutdist_per_seg = 2*(zend-zstart)*tan(cutangle)
    #max length of one crossection
    max_xlength = maxseglength - cutdist_per_seg - 2overlap
    #number of segments needed to keep things short enough
    numsegs = ceil(Int,1u"μm"*span(path)/maxseglength)
    #actual length of one crossection
    xlength = 1u"μm"*span(path)/numsegs
    #create a vector of vectors to store the slices for each block in
    bottomblocks = map(1:numsegs) do _
        Pair[]
    end
    topblocks = deepcopy(bottomblocks)
    zhalf = mean([zstart,zend])
    #build the first set of slices
    curpos = 0u"μm"
    #we will print two layers of blocks
    for (start,stop,blocks) in [(zstart,zhalf + overlap/2,bottomblocks),(zhalf-overlap/2,zend,topblocks)]
        firstsliceoffset = 0u"μm"
        for z in start:dslice:stop
            for i in 1:numsegs
                #if this is the first segment we have to chamfer for overlap on both sides
                startoffset = i == 1 ? firstsliceoffset : zero(firstsliceoffset)
                #if this is the last segment we need to chamfer the end in the opposite direction
                #this should be the same as the firstsliceoffset
                endoffset = i == numsegs ? firstsliceoffset : zero(firstsliceoffset)
                thisseg = subsection(path,curpos + startoffset,curpos + xlength + endoffset + overlap)
                curpos += xlength
                #now we make an offsetcontour based on `profile`
                thiscontour = OffsetContour(thisseg,profile(z)...)
                push!(blocks[i], z => thiscontour)
            end
            #now we'll move up by dslice, so we should slide over a bit to make our cutangle
            curpos -= dslice*tan(cutangle)
            #also update firstsliceoffset, need to do twice as much since we're updating curpos
            firstsliceoffset += 2*dslice*tan(cutangle)
        end
        #move curpos so the seams of our blocks don't overlap
        curpos += xlength/2
    end
    blockcontours = vcat(bottomblocks,topblocks)
    #turn these into blocks
    blocks = map(blockcontours) do bc
        midcontour = bc[floor(Int,length(bc)/2)][2]
        slices = map(bc) do (z,bci)
            z => Tessen.Slice([convert(Contour,bci)])
        end
        thisblock = Block(slices...)
        #move the origin of this block to the centroid of the center slice
        bt = translateorigin(thisblock,centroid(midcontour))
    end
    #build a superblock
    SuperBlock(blocks...)
end

"""
```julia
post(;kwargs...)
```
Build a `Block` representing a triangular
"""
function post(;kwargs...)
    #the z coordinates of our slices
    zcoords = range(start=zero(kwargs[:hbottom]),
                    stop=kwargs[:hbottom]+kwargs[:hbeam],
                    step=kwargs[:dslice])
    #`lengthpairs` will be a z => sidelength iterator
    lengthpairs = map(zcoords) do z
        z => kwargs[:wpost] - (if z <= kwargs[:hbottom]
                                   #we're underneath the beams
                                   #no undercut
                                   zero(kwargs[:wpost])
                               else
                                   #we're on top of the beams
                                   #the amount of 'overcut'
                                   #have to divide by tan(30) to account for how moving the edge
                                   #of the triangle towards its center changes the side length
                                   2*(z-kwargs[:hbottom])*tan(kwargs[:chamfertop])/tand(30)
                               end)
    end

    #helper function to calculate the degree of fillet (from 0 to 1) for the posts
    #when z=0, degreefillet(z) = 1 (fully filleted) when z=:hbottom degreefillet(z)=0
    degreefillet(z) = if z >= kwargs[:hbottom]
        0
    else
        (kwargs[:hbottom] - z)/kwargs[:hbottom]
    end
    #radius of the circumscribed circle of a equilateral triangle with a side length of one
    r = 1/sqrt(3)
    #radius of the corresponding inscribed circle
    r_ins = sqrt(3)/6
    slicepairs = map(lengthpairs) do (z,l)
        #make the feet rounded in the yz/xz plane in addition to filleting in xy
        l = 0.9*l*(1 - degreefillet(z))^(1/2) + 0.1*l
        coords = l .* map((pi/2) .+ [0, 2pi/3, 4pi/3]) do theta
            r*[cos(theta),sin(theta)]
        end
        z => Tessen.Slice([polycontour(coords,degreefillet(z)*l*r_ins)])
    end
    Block(slicepairs...)
end

function oddqhexcoords(q,r)
    #just call an odd-r coordinate system with reversed coordinates and then reverse the result
    map(HexagonOffsetOddR(r,q) |> Hexagons.vertices) do vert
        reverse(vert) |> collect
    end
end

#make wrappers for bridge and post coordinates
struct PostCoords
    p
    lefttilt::Bool
end

struct BeamCoords
    p1
    p2
end

struct HexCoords
    hexsize
    qr #cubic coordinates of hexagon
    verts #cartesian coordinates of verts
end

#if we just want all the vertices
function HexCoords(hexsize,qr)
    temphex = HexCoords(hexsize,qr,nothing)
    av = allverts(temphex)
    HexCoords(hexsize,qr,av)
end

"""
```julia
neighborqrcoords(hc::HexCoords)
```
Get the coordinates (in hexagonal qr coordinates) of hexagons
neighboring `hc`.
"""
function neighborqrcoords(hc::HexCoords)
    #offsets are different if the q coordinate is even or odd
    #https://www.redblobgames.com/grids/hexagons/#neighbors
    offsetvec = if iseven(hc.qr[1])
        [[1,0],
         [1,-1],
         [0,-1],
         [-1,-1],
         [-1,0],
         [0,1]]
    else
        [[1,1],
         [1,0],
         [0,-1],
         [-1,0],
         [-1,1],
         [0,1]]
    end
    
    [hc.qr + a for a in offsetvec]
end

"""
```julia
allverts(hc::HexCoords)
```
Get all vertices of a hexagon (i.e. including ones that may not have a post)
"""
function allverts(hc::HexCoords)
    hc.hexsize * oddqhexcoords(hc.qr...)
end

"""
```julia
edges(hc::HexCoords)
```
Get the edges (as `Vector{Vector}`) of `hc`, presented in the same 'order`
as neighborqrcoords(hc).
"""
function edges(hc::HexCoords)
    v = allverts(hc)
    map(2:7) do i
        [v[i-1],v[(i<7) ? i : i%6]] #we have to 'wrap around' to close the hexagon
    end
end

"""
```julia
center(hc::HexCoords)
```
Get the center of a hex
"""
center(hc::HexCoords) = mean(hc.verts)

#need equality for these
function Base.:(==)(x::PostCoords,y::PostCoords)
    all(isapprox.(x.p,y.p)) && (x.lefttilt == y.lefttilt)
end

#=
function Base.hash(x::PostCoords,h::UInt)
    h = hash(x.p,h)
    hash(x.lefttilt,h)
end
=#

#Base.isequal(x::PostCoords,y::PostCoords) = x == y

function Base.:(==)(b1::BeamCoords, b2::BeamCoords)
    (all(isapprox.(b1.p1,b2.p1)) && all(isapprox.(b1.p2,b2.p2))) || (all(isapprox.(b1.p2,b2.p1)) && all(isapprox.(b1.p1,b2.p2)))
end

#=
function Base.hash(x::BeamCoords,h::UInt)
    h = hash(x.p1,h)
    hash(x.p2,h)
end
=#
#Base.isequal(b1::BeamCoords, b2::BeamCoords) = b1 == b2

"""
inefficient way to get unique posts or beams
"""
function makeunique!(v::Vector)
    i = 1
    while i <= length(v)
        element = v[i]
        #remove any other entries equal to element
        #all indices <i are known to not contain duplicates
        j = i+1
        while j <= length(v)
            if v[j] == element
                deleteat!(v,j)
            else
                j+=1
            end
        end
        i += 1
    end
    return v
end

function mapgeometry(path::Vector{LineEdge},hexsize::Quantity)
    verts = [1u"μm"*le.p1 for le in path]
    #build an ngon to represent our area of interest
    outline = Ngon((Point(v...) for v in verts)...)
    #start at 0,0 and scan column by column, keep going up until we find a hexagon
    #where no vertices are in our outline and we have gone up past the furthest in bounds vert
    #we've ever seen
    maxy = 0
    postcoords = []
    beamcoords = []
    hexcoords = []
    currentcoords = [0,0]
    advancerow = (x) -> x + [1,0]
    #do the right half and then the left
    for _ in 1:2
        while true #scanning rows
            #we will use this to look for a column with no in bounds hexagons
            emptycol = true
            currentcoords[2] = 0
            while true #scanning columns
                theseverts = hexsize * oddqhexcoords(currentcoords...)            
                if any(Point(v...) in outline for v in theseverts)
                    #this hexagon lies in the outline
                    #make PostCoords and BeamCoords objects
                    for (p,lt) in zip(theseverts,repeat([true,false],3))
                        thispost = PostCoords(p,lt)
                        push!(postcoords,thispost)
                    end
                    for i in 1:5
                        push!(beamcoords,BeamCoords(theseverts[i],theseverts[i+1]))
                    end
                    push!(beamcoords,BeamCoords(theseverts[end],theseverts[1]))

                    #also add HexCoords
                    hc = filter(theseverts) do v
                        Point(v...) in outline
                    end
                    if length(hc) >= 3
                        #if only two vertices are in bounds there's no hammock to write
                        push!(hexcoords,HexCoords(hexsize,currentcoords,hc))
                    end
                    
                    emptycol = false
                else
                    #no vertices lie in the outline, stop if we've gone as far as we ever have
                    if currentcoords[2] >= maxy
                        maxy = currentcoords[2]
                        break
                    end
                end
                currentcoords += [0,1]
            end
            #we're done if this column was empty
            if emptycol
                break
            end
            currentcoords = advancerow(currentcoords)
        end
        advancerow = (x) -> x - [1,0]
        currentcoords = [-1,0]
    end
    #====mirror the coords (assuming symmetry around y axis), shouldn't need this now=======
    mirroredpostcoords = [PostCoords(pc.p .* [-1,1], !pc.lefttilt) for pc in postcoords]
    push!(postcoords,mirroredpostcoords...)

    mirroredbeamcoords = [BeamCoords(bc.p1 .* [-1,1], bc.p2 .* [-1,1]) for bc in beamcoords]
    push!(beamcoords,mirroredbeamcoords...)
    
    mirroredhexcoords = [HexCoords([v .* [-1,1] for v in hc.verts]) for hc in hexcoords]
    push!(hexcoords,mirroredhexcoords...)
    ========================================================================================#
    #get all the unique posts and beams where at least one coordinate lies inside the outline
    postcoords = filter(postcoords) do pc
        Point(pc.p...) in outline
    end
    beamcoords = filter(beamcoords) do bc
        (Point(bc.p1...) in outline) || (Point(bc.p2...) in outline)
    end
    #calculate the center of the array.
    centerhexcoords = [0,floor(Int,maxy/2)]
    centerverts = hexsize * oddqhexcoords(centerhexcoords...)
    return (posts=makeunique!(postcoords),beams=makeunique!(beamcoords),center=mean(centerverts),hexcoords=hexcoords)
end

function Tessen.offset(path::Vector{<:Tessen.Edge},x::Real)
    [Tessen.offset(p,x) for p in path]
end

function postbeamcoords(w,h1,h2,h3,cornerdiff,hexsize::Quantity)
    @assert (h1 - 2*cornerdiff) > 0 "cornerdiff is too large"
    ov1 = outlineverts(w + 2*cornerdiff,h1 - 2*cornerdiff,h2,h3)
    ov2 = outlineverts(w,h1,h2,h3)
    bumperpath = polycontour(vcat(ov1[1:2],ov2[3:8])*hexsize).edges
    #get an outline a little bigger than we want
    outer = Tessen.offset(bumperpath,0.1*hexsize)
    #and one a little smaller
    inner = Tessen.offset(bumperpath,-0.1*hexsize)
    #now get the vertices in each
    #nomcoords = mapgeometry(bumperpath,hexsize)
    outercoords = mapgeometry(outer,hexsize)
    innercoords = mapgeometry(inner,hexsize)
    #the posts in outerverts that aren't in innerverts are our attachment points
    attachment = setdiff(outercoords.posts,innercoords.posts)
    return (attachment = attachment, posts = innercoords.posts, beams = innercoords.beams, center = innercoords.center,hexcoords=outercoords.hexcoords)
end

"""
build a segmented beam with nsegs segments centered on y=0.
startx and stopx should be the position of the edge we are bonding to at the z coordinate
corresponding to the bottom of the beam (i.e. this function will build in overlap)
gonna go overboard and make this a struct so we can define rotation and translation
"""
struct Beam
    leftsegs
    rightsegs
    keystone
end

function Beam(nsegs::Int,startx,stopx,width;kwargs...)
    #nsegs needs to be odd
    @assert isodd(nsegs)
    #we want the keystone to be as short as possible. Point of closest approach between the
    #two halves is set by `keygap`
    #length of the keystone measured at the beam midline
    lkey = kwargs[:keygap] + kwargs[:hbeam]*tan(kwargs[:cutangle])
    #the first segment needs to be chamfered according to `chamfertop`, we will do the rest at
    #cutangle
    distperseg = (stopx - startx - lkey)/(nsegs-1)
    lseg = distperseg + kwargs[:overlap]
    leftsegs = map(1:((nsegs-1)/2)) do i
        #the first segment needs to be a little longer and chamfered differently
        segpos = startx + distperseg*((2i-1)/2)
        seg = if i==1
            box(lseg + (kwargs[:overlap]/2),width,kwargs[:hbeam],kwargs[:dslice],
                chamfer =[-kwargs[:chamfertop] kwargs[:cutangle]
                          0                    0])            
        else
            box(lseg,width,kwargs[:hbeam],kwargs[:dslice],
                chamfer =[-kwargs[:cutangle]   kwargs[:cutangle]
                          0                    0])            
        end
        #seg is currently centered at [0,0,0]. Move it into position (use preserveframe so we don't
        #move the stage
        seg=translate(seg,[segpos,0u"µm",kwargs[:hbottom]+kwargs[:hbeam]/2],preserveframe=true)
        if i==1
            #scrootch a little bit back to make the overlap right
            seg=translate(seg,[-kwargs[:overlap]/4,0u"µm",0u"µm"],preserveframe=true)
        end
        return seg
    end
    #the center of the beam in xy plane
    cbeam = [(startx+stopx)/2,0u"µm"]
    #we can make the right side of the bridge by rotating leftsegs around cbeam
    rightsegs = map(leftsegs) do ls
        rotate(ls,pi,cbeam,preserveframe=true)
    end
    #the keystone is cut differently
    keystone = box(lkey,width,kwargs[:hbeam],kwargs[:dslice],
                   chamfer =[-kwargs[:cutangle]   -kwargs[:cutangle]
                             0                    0])
    #move it into position
    keystone = translate(keystone,
                         vcat(cbeam,kwargs[:hbottom]+kwargs[:hbeam]/2),
                         preserveframe=true)
    #return a namedtuple so we can be fancy about how we write these
    Beam(leftsegs,rightsegs,keystone)
end

function Tessen.translate(b::Beam,args...;kwargs...)
    Beam([translate(x,args...;kwargs...) for x in b.leftsegs],
         [translate(x,args...;kwargs...) for x in b.rightsegs],
         translate(b.keystone,args...;kwargs...))
end

function Tessen.rotate(b::Beam,args...;kwargs...)
    Beam([rotate(x,args...;kwargs...) for x in b.leftsegs],
         [rotate(x,args...;kwargs...) for x in b.rightsegs],
         rotate(b.keystone,args...;kwargs...))
end

#use some algorithms from https://www.redblobgames.com/grids/hexagons/ to traverse
#our kernels in a spiral

struct Cube
    q::Int
    r::Int
    s::Int
end

function Base.convert(::Type{Cube},v::Vector{Int})
    @assert length(v)==3
    Cube(v...)
end

function Base.convert(::Type{Vector{Int}},c::Cube)
    [c.q,c.r,c.s]
end

function cube_add(hex, vec)
    Cube(hex.q + vec.q, hex.r + vec.r, hex.s + vec.s)
end

cube_direction_vectors = [
    Cube(+1, 0, -1), Cube(+1, -1, 0), Cube(0, -1, +1), 
    Cube(-1, 0, +1), Cube(-1, +1, 0), Cube(0, +1, -1), 
]

function cube_direction(direction::Int)
    cube_direction_vectors[direction]
end
function cube_neighbor(cube::Cube, direction::Int)
    cube_add(cube, cube_direction(direction))
end

function cube_scale(hex::Cube,factor::Number)
    Cube(hex.q * factor, hex.r * factor, hex.s * factor)
end

function cube_ring(center::Cube, radius::Int)
    results = []
    hex = cube_add(center, cube_scale(cube_direction(5), radius))
    for i in 1:6
        for _ in 1:radius
            push!(results,hex)
            hex = cube_neighbor(hex,i)
        end
    end
    return results
end

"""
Type for iterating through the cubic coordinates of a spiral
traversal of a hexagonal grid. Returns one 'ring' of cubic
coordinates per iteration
"""
struct CubeSpiral
    center::Cube
end

#return the center
function Base.iterate(cs::CubeSpiral)
    return ([cs.center],1)
end

#now work through the rings
function Base.iterate(cs::CubeSpiral,radius::Int)
    return (cube_ring(Cube(0,0,0),radius), radius+1)
end

Base.IteratorSize(::Type{CubeSpiral}) = Base.IsInfinite()

cube_spiral(center,radius::Int) = cube_spiral(convert(Cube,center),radius)

function cubetocartesian(cubecoords::Cube,hexsize)
    hexsize * [3/2         0
               sqrt(3)/2   sqrt(3)] * [cubecoords.q, cubecoords.r]
end

function kernelcenters(center, hexsize, nkernel)
    #get the hex size of our kernel
    kernelsize = nkernel*hexsize
    #make an infinite generator of cartesian coordinates
    ([cubetocartesian(ci,kernelsize) + center for ci in cc] for cc in CubeSpiral(Cube(0,0,0)))
end

"""
move `vert` towards `point` by `amount`

```julia
scrootchvert(vert,point,amount)
```
end
"""
function scrootchvert(vert,point,amount)
    #get a unit vector going from vert to point
    vec = (point - vert)
    uvec = vec/norm(vec)
    vert + amount*uvec
end

function findhexboundary(ourhams::Vector{HexCoords})
    #ok, we need to build a big dict mapping the hexagonal (qr) coordinates
    #of our hexagons to the objects themselves    
    hamdict = Dict(oh.qr => oh for oh in ourhams)
    #find a hammock missing a neighbor (on the outside)
    hamqr = collect(keys(hamdict))

    #helper for fetching neighbor hexcoords objects from hamdict
    neighborhex(h::HexCoords,i::Int) = hamdict[neighborqrcoords(h)[i]]

    #little helper function to see if edge i on a hammock is 'empty'
    function edgeempty(h::HexCoords,i::Int)
        !(neighborqrcoords(h)[i] in hamqr)
    end

    #list of hammocks on the edge of our FOV
    edgeham = filter(ourhams) do oh
        any(edgeempty(oh,i) for i in 1:6)
    end
    
    #mapping for where we should start scanning edges if a given edge is full
    nextedge = [5,6,1,2,3,4]
    
    contours = Contour[]
    while !isempty(edgeham)
        #start at an arbitrary hammock on the edge
        curham = edgeham[1]
        #if all edges on curham are empty, this is an isolated hammock. Add all the edges to a new contour
        if all(edgeempty(curham,i) for i in 1:6)
            push!(contours,Contour([LineEdge(e...) for e in edges(curham)]))
            #delete this from edgeham
            popfirst!(edgeham)
            continue
        end
        #find an empty edge
        curedge = 1
        while edgeempty(curham,curedge)
            curedge += 1
        end
        #we now know our edge isn't empty
        while !edgeempty(curham,curedge)
            curedge += 1
            if curedge > 6
                curedge %= 6
            end
        end
        #now we know we are on the boundary adjacent to another hex
        boundary = []
        startqr = curham.qr    
        startedge = curedge
        while true
            if edgeempty(curham,curedge)
                #curedge is on boundary
                push!(boundary,edges(curham)[curedge])
                curedge += 1
                #test to make sure we don't need to roll over
                if curedge > 6
                    curedge = curedge%6
                end
            else
                #curedge is not on boundary, walk to neighbor
                curham = neighborhex(curham,curedge)
                curedge = nextedge[curedge]
            end
            #check to see if we've made it all the way around
            if (curham.qr == startqr) && (startedge==curedge)
                break
            end
        end
        #make a Contour to return
        push!(contours,Contour([LineEdge(b...) for b in boundary]))
        #now remove all hexes contained by this contour from edgeham
        contourverts = [(b[1] for b in boundary)...,boundary[end][2]]
        contouroutline = Ngon((Point(v...) for v in contourverts)...)
        filter!(edgeham) do h
            !(Point(center(h)...) in contouroutline)
        end
    end
    return contours
end

"""
```julia
scaffold(scaffolddir[, configfilename])
scaffold(scaffolddir,configdict)
```
Build a directory of .gwl files to build a scaffold in a directory at `scaffolddir`. The
scaffold's geometrical parameters can be provided as a `Dict` or read from `configfilename`.
If `configfilename` is omitted, the parameters will be read from `"config.jl"`
"""
function scaffold end

function scaffold(scaffolddir,kwargs::Dict)
    mkdir(scaffolddir)
    kernels = cd(scaffolddir) do
        mkdir("scripts")
        #build the bumper
        #we will make two outlines, one stubbier than the other
        @assert (kwargs[:h1] - 2*kwargs[:cornerdiff]) > 0 "cornerdiff is too large"
        ov1 = outlineverts(kwargs[:w] + 2*kwargs[:cornerdiff],kwargs[:h1] - 2*kwargs[:cornerdiff],kwargs[:h2],kwargs[:h3])
        ov2 = outlineverts(kwargs[:w],kwargs[:h1],kwargs[:h2],kwargs[:h3])
        #this will give us asymmetry so we can easily verify orientation
        ov = vcat(ov1[1:2],ov2[3:8])
        p = polycontour(ov*kwargs[:hexsize],kwargs[:bumperfillet]).edges
        bec = bumperedgecoords(;kwargs...)
        #largest section of bumper which can be printed at a time
        #maxbumperseg = sqrt(kwargs[:dfield]^2 - kwargs[:wbumper]^2)
        b = Toast.bumper(p,bec,kwargs[:bumperseglength],zero(kwargs[:hbottom]),
                         kwargs[:hbottom]+kwargs[:htop],kwargs[:dslice],kwargs[:overlap],
                         kwargs[:cutangle])
        @info "hatching bumper"
        bh = hatch(b,kwargs[:dhatch],1,pi/2)
        @info "compiling bumper"
        compbumper = CompiledGeometry(joinpath("scripts","bumper.gwl"),bh;laserpower=kwargs[:laserpower],scanspeed=kwargs[:scanspeed])
        #get the coordinates of all the posts and beams
        pc = Toast.postbeamcoords(kwargs[:w],kwargs[:h1],kwargs[:h2],kwargs[:h3],kwargs[:cornerdiff],kwargs[:hexsize])

        #remove hammocks which don't have 6 vertices (i.e. that are incomplete)
        filter!(pc.hexcoords) do h
            length(h.verts) == 6
        end

        #I now want to build two dicts, one which tracks whether a post has been written or not
        #and one which maps beams to the corresponding posts

        #first build the postdict, marking all posts as unwritten
        postdict = Dict(p => false for p in pc.posts)

        #now build a dict mapping beams to posts
        beamdict = Dict(map(pc.beams) do b
                            posts = map([b.p1, b.p2]) do coord
                                pvec = filter(pc.posts) do post
                                    all(isapprox.(coord,post.p))                                
                                end
                                if length(pvec) == 1
                                    return pvec[1]
                                else
                                    @assert length(pvec) == 0
                                    #make sure this is in attachments
                                    attvec = filter(pc.attachment) do a
                                        all(isapprox.(coord,a.p))
                                    end
                                    @assert length(attvec) == 1
                                    return nothing
                                end
                            end
                            b => posts
                        end...)
        #figure out how big one kernel is in terms of hex size
        rfield = kwargs[:dfield]/2
        buffer = 1.25 #expand our kernel FOV a little to avoid numerical precision issues
        nhexkernel = floor(Int,rfield/(kwargs[:hexsize]*buffer))
        #get our kernel radius in microns
        rkernel = buffer * kwargs[:hexsize] * nhexkernel
        #now, start in the center and spiral outward writing geometry until we run out of posts
        #write a beam whenever it's supporting posts are both written
        firstcenter = pc.center
        #get a generator for kernel centers
        kcgen = kernelcenters(firstcenter, kwargs[:hexsize], nhexkernel)
        #build all the kernels
        kernels = [compbumper]
        #keep building kernels until we spiral past the end of the scaffold
        #this generator is infinite
        #build a generic post that we can translate/rotate to make the others
        genericpost = post(;kwargs...)
        #keep track of where we are in the spiral
        lastcenter = zero(firstcenter)
        #collect all the kernel centers we are going to use so we can reverse the spiral
        #I know I'm copy-pasting from the main kernel building loop in a way that is hacky
        #but the Danes are here and I want a working job
        kcvec=[]
        for (i,ring) in enumerate(kcgen)
            emptyring = true
            for (j,kernelcenter) in enumerate(ring)
                coordinkernel(c) = norm(c - kernelcenter) <= rkernel
                #we may need to worry about big posts being outside the FOV at some point
                ourposts = filter(keys(postdict)|>collect) do pi
                    coordinkernel(pi.p)
                end
                if !isempty(ourposts)
                    emptyring = false
                else
                    continue
                end
            end
            if !emptyring
                push!(kcvec,(i,ring))
            else
                break
            end
        end
        
        for (i,ring) in reverse(kcvec)
            emptyring = true
            for (j,kernelcenter) in enumerate(ring)
                @info "kernel $j on ring $i"
                #find all the unwritten posts in range
                coordinkernel(c) = norm(c - kernelcenter) <= rkernel
                #we may need to worry about big posts being outside the FOV at some point
                ourposts = filter(keys(postdict)|>collect) do pi
                    coordinkernel(pi.p)
                end
                if !isempty(ourposts)
                    emptyring = false
                else
                    continue
                end
                #build blocks for all these posts
                postblocks = map(ourposts) do p_i
                    if postdict[p_i]
                        #already wrote this one
                        return nothing
                    end
                    #mark as written
                    postdict[p_i] = true
                    #first do the rotation
                    rotpost = rotate(genericpost,p_i.lefttilt ? pi/6 : -pi/6,preserveframe=true)
                    #then translate
                    translate(rotpost,p_i.p - kernelcenter,preserveframe=true)
                end
                #merge
                filter!(postblocks) do pb
                    !isnothing(pb)
                end
                #no need to merge if there's one or zero posts
                mergedposts = (length(postblocks)>1) ? [merge(postblocks...)] : postblocks
                #check for beams we can write
                ourbeams = filter(collect(keys(beamdict))) do bi
                    all(beamdict[bi]) do pi
                        #check if both supporting posts have been written
                        isnothing(pi) || postdict[pi]
                    end
                end
                #build beam blocks
                beamblocks = map(ourbeams) do bi
                    #delete from beamdict so we don't rewrite
                    delete!(beamdict,bi)
                    #vector along the beam
                    bvec = bi.p2 - bi.p1
                    #get a unit vector along the beam
                    ubeam = bvec/norm(bvec)
                    #we don't want to go from center to center, we want to go edge to edge.
                    centertoedgedist = kwargs[:wpost]*tand(30)/2
                    startpoint = bi.p1 + ubeam*centertoedgedist
                    endpoint = bi.p2 - ubeam*centertoedgedist
                    bvec = endpoint - startpoint
                    #get the length of the beam (should always be the same)
                    beamlength = norm(bvec)
                    #get the angle that the vector from p1 to p2 makes with the x axis
                    theta = atan(reverse(bvec)/beamlength...)
                    #build a horizontal beam with the correct length
                    nsegs = ceil(Int,beamlength/kwargs[:maxseglength])
                    nsegs = isodd(nsegs) ? nsegs : nsegs+1
                    #build a post along the x axis
                    #the width parameter to beam corresponds to the width at the midline
                    #we want the width to be wpost at the bottom
                    wbeam = kwargs[:wpost] - kwargs[:hbeam]*tan(kwargs[:chamfertop])
                    hbeam = Beam(nsegs,zero(beamlength),beamlength,wbeam;kwargs...)
                    #rotate
                    rotbeam = rotate(hbeam,theta,preserveframe=true)
                    #translate
                    translate(rotbeam,startpoint - kernelcenter,preserveframe=true)
                end
                #merge beam segments
                #we now want to merge all of our beams so the segments print starting from
                #the posts
                #get a vector containing each 'half' beam
                halfbeams = vcat(map(beamblocks) do b
                                     [b.leftsegs, b.rightsegs]
                                 end...)
                halfbeamblocks = iszero(length(halfbeams)) ? [] : map(zip(halfbeams...)) do hb
                    #hb should now be all of the first segs, or all of the second segs, etc
                    merge(hb...)
                end
                #now need to collect all of the keystones
                allkeys = map(beamblocks) do b
                    b.keystone
                end
                
                #and merge them
                keyblock = iszero(length(allkeys)) ? [] : [merge(allkeys...)]
                                 
                #build a SuperBlock
                allblocks = vcat(mergedposts,halfbeamblocks,keyblock)               
                if isempty(allblocks)
                    #no geometry to write
                    continue
                end
                postbeams = SuperBlock(allblocks...)
                
                #build hammocks (only need to do this if we printed other geometry, hence continue above
                #helper function to see if a hex vertex has it's support written
                function supportwritten(vertex)
                    thispostcoord = filter(keys(postdict)  |> collect) do pk
                        #use isapprox to compare coordinates elementwise to avoid numerical precision issues
                        all(isapprox(c,v) for (c,v) in zip(pk.p,vertex))
                    end
                    @assert length(thispostcoord) <= 1
                    #if thispostcoord is empty, that means that post is an attachment point, so we can
                    #return true
                    isempty(thispostcoord) ? true : postdict[thispostcoord[1]]
                end

                ourhams = filter(pc.hexcoords) do h
                    all(h.verts) do v
                        supportwritten(v)
                    end
                end

                #remove these hammocks from pc.hexcoords so we don't rewrite
                filter!(pc.hexcoords) do h
                    !all(h.verts) do v
                        supportwritten(v)
                    end
                end
                
                boundary = findhexboundary(convert(Vector{HexCoords},ourhams))
                #need to move the boundary to account for the current objective position
                boundary = [translate(b,-1*kernelcenter) for b in boundary]
                #start from :hamoffset above the top of the beams
                thisz = kwargs[:hbottom]+kwargs[:hbeam]+kwargs[:hamoffset]
                slices = map(1:kwargs[:nhammock]) do _
                    toreturn = thisz => Tessen.Slice(boundary)
                    thisz += kwargs[:dhammockslice]
                    return toreturn
                end
                floorblock = Block(slices...)                
                #hatch
                @info "hatching kernel $i-$j"
                hatched = hatch(postbeams,kwargs[:dhatch],0,pi/2)
                #i think we will always have hammocks to write
                @info "compiling kernel $i-$j"                
                #compile
                cg = CompiledGeometry(joinpath("scripts","kernel_($i-$j).gwl"),hatched;laserpower=kwargs[:laserpower],scanspeed=kwargs[:scanspeed])
                cg = translate(cg,[kernelcenter...,zero(kernelcenter[1])])
                push!(kernels,cg)
                #ch = translate(ch,[kernelcenter...,zero(kernelcenter[1])])
                #cg = translate(cg,[(kernelcenter-lastcenter)...,zero(kernelcenter[1])])
                if !isempty(boundary)
                    hamhatched = hatch(floorblock,kwargs[:dhammockhatch],1,pi/2)
                    ch = CompiledGeometry(joinpath("scripts","kernel_($i-$j)_hammocks.gwl"),hamhatched;laserpower=kwargs[:hamlaserpower],scanspeed=kwargs[:hamscanspeed])
                    ch = translate(ch,[kernelcenter...,zero(kernelcenter[1])])
                    push!(kernels,ch)
                end
                lastcenter=kernelcenter
            end
            if emptyring
                #if we didn't find any posts in this ring
                break
            end
        end
        return kernels
    end
    GWLJob(joinpath(scaffolddir,"scaffold.gwl"),kernels...;stagespeed=kwargs[:stagespeed],interfacepos=kwargs[:interfacepos])
end

function scaffold(scaffolddir,configfilename::String)
    config = include(configfilename)
    scaffold(scaffolddir,config)
end

scaffold(scaffolddir) = scaffold(scaffolddir,"config.jl")

"""
```julia
arrangescaffolds(arraydims,arrayshape,arraycenter,maxscafdims)
```
Return a `Matrix` with shape `arrayshape` of scaffold coordinate centers centered
on `arraycenter`. Given a maximum scaffold dimension of `maxscafdims` these scaffolds
will be guarenteed not to overlap and to fit in a bounding box with size `arraydims`.
"""
function arrangescaffolds(arraydims,arrayshape,arraycenter,maxscafdims)
    #get the corners of our bounding box
    bboxtopleft = arraycenter + [-arraydims[1], arraydims[2]]/2
    bboxbottomright = arraycenter + [arraydims[1], -arraydims[2]]/2

    #get the coordinates of our top left and bottom right scaffold
    firstcenter = bboxtopleft + [maxscafdims[1],-maxscafdims[2]]/2
    lastcenter = bboxbottomright + [-maxscafdims[1],maxscafdims[2]]/2
    #this is enough to build the matrix
    (xrange,yrange) = map(1:2) do d
	if arrayshape[d] == 1
	    #put it in the middle
	    return [mean([firstcenter[d],lastcenter[d]])]
	end
        #if we have more than one scaffold along this dimension make a range
        range(start=firstcenter[d],stop=lastcenter[d],length=arrayshape[d])
    end
    if xrange isa AbstractRange
	@assert step(xrange) > maxscafdims[1] "scaffolds would overlap in x direction"
    end
    if yrange isa AbstractRange
	@assert (-step(yrange)) > maxscafdims[2] "scaffolds would overlap in y direction"
    end
    centers = [[x,y] for x in xrange, y in yrange]
end

#build a multijob from a matrix of `centercoords => config` pairs
function buildmj(jobs::Matrix{<:Pair})
    #snake it
    stagespeed = nothing
    rowjobs = map(1:size(jobs)[2]) do j
        thisrow = jobs[:,j]
        if iseven(j)
            thisrow = reverse(thisrow)
        end
        thesejobs = map(1:length(thisrow)) do i
            (center,config) = thisrow[i]
            #assume stagespeed is always the same
            stagespeed = config[:stagespeed]
            thisjob = (center => scaffold("$i-$j",config))
            #write the configuration into the scaffold folder so we can look at it later
            open(joinpath("$i-$j","config.txt"), "w") do io
                print(io,config)
            end
            return thisjob
        end
    end
    multijob("psweep.gwl",vcat(rowjobs...)...;stagespeed)
end

"""
```julia
psweep([config,] p1 => values[, p2 => values]; arraydims, arrayshape, arraycenter,scafdims)
```
Build a `multijob` which builds scaffolds with varying parameters. The final array will have
shape `arrayshape` centered on `arraycenter`. These scaffolds will be guarenteed not to overlap
and to fit in a bounding box with size `arraydims`. `config` (provided as a filepath or `Dict`)
should contain all other configuration parameters. If `config` is not provided, a configuration
file is assumed to exist at `"config.jl"`. If swept parameters are included in `config` they will
be ignored.
"""
function psweep end

#for one parameter
function psweep(config::Dict,sweep::Pair{Symbol,<:Vector}; arraydims,arrayshape,arraycenter,scafdims)
    #destructure our parameter and values
    (p,values) = sweep
    #build a vector of configurations reflecting our parameter sweep
    configs = map(values) do v
        thisconfig = copy(config)
        thisconfig[p] = v
        return thisconfig
    end

    maxscafdims = map(scafdims) do dim
        maximum(configs) do c
            c[dim]
        end
    end
    scafcenters = arrangescaffolds(arraydims,arrayshape,arraycenter,maxscafdims)
    @assert length(configs) == length(scafcenters) "Number of parameter values must match array size"
    jobmat = map(zip(scafcenters,reshape(configs,size(scafcenters)...))) do (sc,c)
        sc => c
    end
    buildmj(jobmat)
end

#for two parameters
function psweep(config::Dict,sweep1::Pair{Symbol,<:Vector},sweep2::Pair{Symbol,<:Vector};
                arraydims,arrayshape,arraycenter)
    #destructure our parameter and values
    (p1,values1) = sweep1
    (p2,values2) = sweep2
    #build a matrix representing our parameter combos
    pmat = [Dict(p1 => v1, p2 => v2) for v1 in values1, v2 in values2]
    #now make a matrix of configs
    configs = map(pmat) do pdict
        thisconfig = copy(config)
        for (p,v) in pdict
            thisconfig[p] = v
        end
        return thisconfig
    end

    maxscafdims = map(scafdims) do dim
        maximum(configs) do c
            c[dim]
        end
    end
    scafcenters = arrangescaffolds(arraydims,arrayshape,arraycenter,maxscafdims)
    @assert length(configs) == length(scafcenters) "Number of parameter values must match array size"
    jobmat = map(zip(scafcenters,configs)) do (sc,c)
        sc => c
    end
    buildmj(jobmat)
end

function psweep(config::String,args...;kwargs...)
    cdict = include(config)
    psweep(cdict,args...;kwargs...)
end

psweep(args::Vararg{<:Pair{Symbol,<:Vector}};kwargs...) = psweep("config.jl",args...;kwargs...)

"""
```julia
repjob([config]; arraydims, arrayshape, arraycenter,scafdims)
```
Create a job to write an array of identical scaffolds. The array will have shape `arrayshape`
centered on `arraycenter`. These scaffolds will be guarenteed not to overlap and to fit in a
bounding box with size `arraydims`. `config` (provided as a filepath or `Dict`) should contain
all configuration parameters.
"""
function repjob end

function repjob(config::Dict; arraydims,arrayshape,arraycenter,scafdims)
    scafcenters = arrangescaffolds(arraydims,arrayshape,arraycenter,scafdims)
    #snakify
    rows = map(1:size(scafcenters)[2]) do j
        iseven(j) ? reverse(scafcenters[:,j]) : scafcenters[:,j]
    end
    #we're going to do one job over and over
    job = scaffold("scaffold",config)
    multijob("repjob.gwl",(c => job for c in vcat(rows...))...;
             stagespeed=config[:stagespeed])
end

function repjob(config::AbstractString;kwargs...)
    cdict = include(config)
    repjob(cdict;kwargs...)
end

repjob(;kwargs...) = repjob("config.jl";kwargs...)

end # module Toast
