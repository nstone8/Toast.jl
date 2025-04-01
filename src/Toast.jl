module Toast

using Hexagons, Delica, Tessen, Unitful, Statistics
import Base:convert

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
    span(loc.edge)*1u"μm" * abs(x1-x2)
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
    alpha = aoc.stopangle - aoc.startangle
    #make sure this is positive
    alpha = aoc < 0 ? aoc + 2pi : aoc

    #get the edges of the contour
    edges = getoffsetedges(aoc)
    r = sort([1u"μm"*e.r for e in edges])

    #get the distance between the center of curvature and the centroid
    d = (2/3)*(sin(alpha)/alpha)*(r[2]^3 - r[1]^3)/(r[2]^2 - r[1]^2)

    #now get the direction
    theta1 = aoc.startangle
    theta2 = aos.stopangle
    #theta2 needs to be greater than theta1
    theta2 = theta1 > theta2 ? theta2 + 2pi : theta2

    theta = mean([theta1,theta2])
    [d*cos(theta), d*sin(theta)]
end

function area(aoc::ArcOffsetContour)
    #total arc
    alpha = aoc.stopangle - aoc.startangle
    #make sure this is positive
    alpha = aoc < 0 ? aoc + 2pi : aoc

    #get the edges of the contour
    edges = getoffsetedges(aoc)
    r = sort([1u"µm"*e.r for e in edges])
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

function area(aoc::CompoundOffsetContour)
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
    otheredges = map(convert.(Contour,c[2:end-1],startcap=false,endcap=false)) do c
        c.edges
    end
    
    Contour(vcat(firstedges,vcat(otheredges...),lastedges))
end

#add a span method for vectors of edges
function Tessen.span(edges::Vector{Edge})
    sum(edges) do e
        span(e)
    end
end
#add a subsection method for vectors of edges
function Tessen.subsection(edges::Vector{Edge},x1::Real,x2::Real)
    @assert all([x1,x2] .< span(edges)) "Both x1 and x2 must be less than total span"
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
    output = Edge[]
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

end # module Toast
