def create_app_instance():
    app = Sanic(name = "FESapp")
    app.config.from_object('FESapp.settings')

    return app