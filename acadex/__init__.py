from flask import Flask
from .ui import blueprint as ui_blueprint
from .config import Config
from lbrc_flask import init_lbrc_flask, ReverseProxied
from lbrc_flask.security import init_security, Role, User


def create_app(config=Config):
    app = Flask(__name__)
    app.wsgi_app = ReverseProxied(app.wsgi_app)
    app.config.from_object(config)

    TITLE = 'Acadex'

    with app.app_context():
        init_lbrc_flask(app, TITLE)

        init_security(app, user_class=User, role_class=Role)

    app.register_blueprint(ui_blueprint)

    return app
