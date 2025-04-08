module Toast

using Hexagons, Delica, Tessen, Unitful, Statistics, Meshes
import Base:convert
import Tessen:HatchLine, pointalong, intersections

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
        :hamlaserpower => 50u"mW"
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
    meanoffset = offset(loc.edge,mean([loc.x1,loc.x2]))
    #get the midpoint
    pointalong(meanoffset,span(meanoffset)/2) * 1u"μm"
end

function area(loc::LineOffsetContour)
    #this is the length of the line times the difference between the offsets
    span(loc.edge)*1u"μm" * abs(loc.x1-loc.x2)
end

function getoffsetedges(loc::LineOffsetContour)
    [offset(loc.edge,x) for x in [loc.x1,loc.x2]]
end

struct ArcOffsetContour <: OffsetContour
    edge::ArcEdge
    x1::Quantity
    x2::Quantity
end

OffsetContour(ae::ArcEdge,x1,x2) = ArcOffsetContour(ae,x1,x2)

function getoffsetedges(aoc::ArcOffsetContour)
    [offset(aoc.edge,x) for x in [aoc.x1,aoc.x2]]
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
function convert(T::Type{Contour},oc::OffsetContour;startcap=true,endcap=true)
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

function convert(T::Type{Contour}, coc::CompoundOffsetContour)
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
function Tessen.span(edges::Vector{Tessen.Edge})
    sum(edges) do e
        span(e)
    end
end
#add a subsection method for vectors of edges. We will assume this corresponds to closed paths
function Tessen.subsection(edges::Vector{Tessen.Edge},x1::Real,x2::Real)
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
            z => Slice([convert(Contour,bci)])
        end
        thisblock = Block(slices...)
        #move the origin of this block to the centroid of the center slice
        bt = translateorigin(thisblock,centroid(midcontour))
        hatch(bt,dhatch=1u"μm",bottomdir=pi/4,diroffset=pi)
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
        z => kwargs[:wpost] - (if z < kwargs[:hbottom]
                                   #we're underneath the beams
                                   #no undercut
                                   zero(kwargs[:wpost])
                               else
                                   #we're on top of the beams
                                   #the amount of 'overcut'
                                   2*(z-kwargs[:hbottom])*tan(kwargs[:chamfertop])
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
        @show l
        coords = l .* map((pi/2) .+ [0, 2pi/3, 4pi/3]) do theta
            r*[cos(theta),sin(theta)]
        end
        z => Slice([polycontour(coords,degreefillet(z)*l*r_ins)])
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

function isequal(b1::BeamCoords, b2::BeamCoords)
    (isequal(b1.p1,b2.p1) && isequal(b1.p2,b2.p2)) || (isequal(b1.p2,b2.p1) && isequal(b1.p1,b2.p2))
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
    currentcoords = [0,0]
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
        currentcoords += [1,0]
    end
    #mirror the coords (assuming symmetry around y axis)
    mirroredpostcoords = [PostCoords(pc.p .* [-1,1], !pc.lefttilt) for pc in postcoords]
    push!(postcoords,mirroredpostcoords...)

    mirroredbeamcoords = [BeamCoords(bc.p1 .* [-1,1], bc.p2 .* [-1,1]) for bc in beamcoords]
    push!(beamcoords,mirroredbeamcoords...)
    
    #get all the unique posts and beams where at least one coordinate lies inside the outline
    postcoords = filter(postcoords |> Set |> collect) do pc
        Point(pc.p...) in outline
    end
    beamcoords = filter(beamcoords |> Set |> collect) do bc
        (Point(bc.p1...) in outline) || (Point(bc.p2...) in outline)
    end
    return (posts=postcoords,beams=beamcoords)
end

function Tessen.offset(path::Vector{<:Tessen.Edge},x::Real)
    [Tessen.offset(p,x) for p in path]
end

function postbeamcoords(w,h1,h2,h3,hexsize::Quantity)
    bumperpath = polycontour(outlineverts(w,h1,h2,h3)*hexsize).edges
    #get an outline a little bigger than we want
    outer = Tessen.offset(bumperpath,0.1*hexsize)
    #and one a little smaller
    inner = Tessen.offset(bumperpath,-0.1*hexsize)
    #now get the vertices in each
    outercoords = mapgeometry(outer,hexsize)
    innercoords = mapgeometry(inner,hexsize)
    #the posts in outerverts that aren't in innerverts are our attachment points
    attachment = setdiff(outercoords.posts,innercoords.posts)
    return (attachment = attachment, posts = innercoords.posts, beams = innercoords.beams)
end

end # module Toast
