#!julia

import Pkg

Pkg.activate(Base.current_project())

cd("notebooks")

import Pluto

Pluto.run(;
    host = "0.0.0.0",
    port = 5999,
    launch_browser = "--browser" in Base.ARGS,
    enable_ai_editor_features = false,
)
