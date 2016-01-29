from pkg_resources import get_distribution, DistributionNotFound

__project__ = 'Cogent'
__version__ = None  # required for initial installation

try:
        __version__ = get_distribution(__project__).version
except DistributionNotFound:
        __version__ = 'Please install this project with setup.py'

def get_version():
    return __project__ + ' ' + str(__version__)
