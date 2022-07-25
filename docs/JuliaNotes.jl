using Printf 
using Plots
using Revise


#using BenchmarkTools
# @benchmark
# @btime

# Does Gmsh have Prismatic elements support?
# look at taylor methods

# VS code config
# Focus active editor group: alt+j alt+i
# Focus REPL: alt+j alt+o

#quick edit a function from REPL
# @edit funciton()

#Fun plots
using UnicodePlots
unicodeplots();
plot(Plots.fakedata(50, 5), w = 3)


#=
REPL notes:
  ? - help mode
  ; - shell mode
  <backspace> - leave mode
  ctrl-L - clears output window
=#

# Cool animated plots
@gif for i in vals
    plot(...use your i here like you normally would...)
end

begin 
	using Plots
	x = 1:10; y = rand(10); # These are the plotting data
	plot(x, y)
end

f = x-> x^2+2x


function test(a)
	a+2
end

# Composision
(test ∘ f)(3)

# pipeing
begin
	1:10|>sum
	["a","test","Mason"] .|> [uppercase, reverse, length]
	["a","test","Mason"] .|> (uppercase ∘ reverse)
end

#Special characters
[i^2 for i ∈ 1:10] # \in<tab>
# ∘  \circ<tab>


begin
	A = [1,3,2]
	sin.(A)
end

begin
	g=(1/n^2 for n=1:100) # This is a generator the meat of a list comprehension without allocating memory
	sum(g)
end

# Static arrays!
# Julia Cache
# @Code_warntype ... Check type stability
# tim holy
# @muladd


#Cool packages: ProgressBar


# indexing starts at 1!!!
begin
	B = [2,3,4,5]
	B[1]
end
