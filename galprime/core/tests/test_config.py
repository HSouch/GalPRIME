from .. import config

from configobj import ConfigObj


def test_default_config():
    default_config  = config.default_config()
    assert isinstance(default_config, ConfigObj)


def test_dump_default_config_file(tmpdir):
    filename = f"{tmpdir}default.gprime"
    config.dump_default_config_file(outname=filename)

    test_read_config = config.read_config_file(filename)
    assert isinstance(test_read_config, ConfigObj)