import Pkg
using TOML
@info "Activating temporary environment and installing PackageCompiler.jl"
Pkg.activate(; temp = true)
Pkg.add("PackageCompiler")

repo_root = normpath(@__DIR__, "..")
@info "Moving to $repo_root"
cd(repo_root)

ENV["JULIA_DEBUG"] = "Pkg,PackageCompiler"

using PackageCompiler

@static if Sys.iswindows()
    extension = "dll"
elseif Sys.isapple()
    extension = "dylib"
else
    extension = "so"
end

sysimage_path = "project-sysimage.$extension"
precompile_statements_file = "compile-statements.jl"
project = pwd()

# Get the packages as the ones explicitly in the Project.toml. This excludes self.
packages = TOML.parsefile(joinpath(project, "Project.toml"))["deps"] |> keys |> collect
@debug "Putting the following in the sysimage" packages

import_into_main = false

create_sysimage(packages; project, sysimage_path, precompile_statements_file, import_into_main)
