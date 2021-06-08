using Documenter
using Movement

makedocs(
    sitename = "Movement",
    format = Documenter.HTML(),
    modules = [Movement]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
