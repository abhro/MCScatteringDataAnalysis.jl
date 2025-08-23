#!julia

import Pluto

cd("notebooks") do
    Pluto.run(;
        host = "0.0.0.0",
        port = 5999,
        launch_browser = "--browser" in Base.ARGS,
        enable_ai_editor_features = false,
    )
end
