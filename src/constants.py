import yaml as yy
from yaml.loader import Loader


def get_config(config_file, cr_point=None, tw=None, group=None):
    s = open(config_file, 'r')
    cfs = yy.load(s, Loader=Loader)
    globals().update(cfs)
