import Contour: contours, levels, level, lines, coordinates
using PlotlyJS
using Plots
using AssociatedLegendrePolynomials

# (0,θ,Φ) is considered the nucleus of the atom. 
# r is always the distance between the nucleus and the electron.

# Create Input
data = range(-1, stop = 1, length =  60)
x, y, z = mgrid(data, data, data)

# Transition Coordinates
r = @. √(x^2 + y^2 + z^2) 
θ = @. atan(y, x)
Φ = @. acos(z/(√(x^2 + y^2 + z^2)))

# Constants
a = 0.052941 #nm
e = ℯ

# Set Quantum Numbers
n = 2
l = 1
m = 0
# Set Nucleus Charge
Z = 1

# Associated Laguerre Polynomial Generator Lₖʲ(x)
function L(j, k, x) 
    f = 0;

    for i in range(0, j) 
        f += @. (-1)^i * factorial(j + k)/(factorial(j - i)*factorial(k+i)*factorial(i)) * x^i
    end

    return f
end

# Put Back the Normalization Constants in Each Part of the Wave Function
ψ(r) = @. √(((factorial(n-l-1))/(2*n*factorial(n+l))) * (2/(n*a))^3) * (((2Z*r)/(n*a))^(l)) * (e^((-Z*r)/(n*a))) * L(n-l-1, 2l+1, ((2Z*r)/(n*a)))
Θ(θ) = @. (-1)^m * Nlm(l,m) * Plm(l,m, cos(θ))
P(Φ) = @. e^(im * Φ * m)

# Calculate the Energy of the electron
E = @. (-13.6056)/(n^2)

# Combined, Full Wave Function
Ψ(r, θ, Φ) = @. ψ(r) * Θ(θ) * P(Φ)

# Calculate the Expactation Value of the electron
expectation = @. abs(Ψ(r, θ, Φ))^2

# Creates Annotation for Energy Render
layout = Layout(
    scene=attr(
        annotations=[
        attr(
            showarrow=false,
            x=0,
            y=0,
            z=1,
            text= "$E eV",
            textangle=0,
            ax=0,
            ay=-75,
            font=attr(
                color="#EE4A77",
                size=24
            )
            )
        ]
    ),
)

PlotlyJS.plot(volume(
    x=x[:],
    y=y[:],
    z=z[:],
    value=expectation[:],
    isomin=0.05,
    isomax=0.7,
    opacity=0.1, # needs to be small to see through all surfaces
    surface_count=100, # needs to be a large number for good volume rendering
), layout)