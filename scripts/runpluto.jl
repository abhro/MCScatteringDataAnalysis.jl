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
port = 5999
launch_browser = "--browser" in ARGS

Pluto.run(;
    host,
    port,
    root_url = "http://$hostname:$port/",
    launch_browser,
    auto_reload_from_file = true,
    enable_ai_editor_features = false,
)
