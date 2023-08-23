"""
WSGI config for FESapp project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/4.1/howto/deployment/wsgi/
"""

# import os

# from django.core.wsgi import get_wsgi_application

# os.environ['DJANGO_SETTINGS_MODULE'] = 'FESapp.settings'

# #application = get_wsgi_application()
# from app import app as application
# app = application

from module import create_app_instance


app = create_app_instance()


if __name__ == '__main__':
    app.run()