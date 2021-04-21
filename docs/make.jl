using Documenter
using myED

makedocs(
    sitename = "myED",
    format = Documenter.HTML(),
    modules = [myED]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
