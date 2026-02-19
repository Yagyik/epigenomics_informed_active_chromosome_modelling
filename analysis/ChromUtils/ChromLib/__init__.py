import os

# Dynamically generate the __all__ list based on the files in the current directory

def __init__():
    all = [os.path.splitext(f)[0] for f in os.listdir(os.path.dirname(__file__)) if f.endswith('.py') and f != '__init__.py']