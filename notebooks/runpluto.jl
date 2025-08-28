#!julia

import Pluto

Pluto.run(;
    host = "0.0.0.0",
    port = 5999,
    launch_browser = "--browser" in Base.ARGS,
    auto_reload_from_file = true,
    enable_ai_editor_features = false,
)
