#!julia

using DrWatson
@quickactivate "MCScatteringDataAnalysis"

if "--help" in ARGS
    println("Usage: $PROGRAM_FILE [--network] [--browser]")
    exit()
end

import Pluto

cd(projectdir("notebooks"))

if "--network" in ARGS
    # where to listen
    host = "0.0.0.0"
    # name of computer
    hostname = Libc.gethostname()
else
    host = "127.0.0.1"
    hostname = "localhost"
end
# port to listen to, taken from thesis course number
#port = 5999
launch_browser = "--browser" in ARGS

# should we use a sysimage?
jloptions = Base.JLOptions()
if jloptions.image_file_specified == 1
    sysimage = unsafe_string(jloptions.image_file)
else
    sysimage = nothing
end
@debug "Running with sysimage: $sysimage"

Pluto.run(;
    host,
    #port_hint = port,
    #root_url = "http://$hostname:$port/",
    sysimage,
    launch_browser,
    depwarn = "error",
    auto_reload_from_file = true,
    enable_ai_editor_features = false,
)
