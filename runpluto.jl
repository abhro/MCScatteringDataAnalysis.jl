#!julia

import Pkg

Pkg.activate(Base.current_project())

import Pluto

launch_browser = "--browser" in Base.ARGS

Pluto.run(;
    host = "0.0.0.0",
    port = 5999,
    launch_browser,
)
