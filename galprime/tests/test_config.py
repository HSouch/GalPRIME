from .. import config


def test_default_config():
    default_config  = config.default_config()


def test_dump_default_config_file(tmpdir):
    filename = f"{tmpdir}default.gprime"
    config.dump_default_config_file(outname=filename)

    config.read_config_file(filename)