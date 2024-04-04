import yaml


def load_config_file(filename):
    try:
        with open('filename', 'r') as f:
            conf = yaml.safe_load(f)
        return conf
    except FileNotFoundError:
        print(f"File {filename} not found.")
        return None
    except Exception as e:
        print(f"An error occurred loading {filename}: {e}")
        return None