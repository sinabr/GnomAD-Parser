from pathlib import Path
import json

def read_external_config(external_config_filepath:str|Path)->dict:
    """
    read the external config filepath

    Parameters:
    -----------
    external_config_filepath : str|Path : Path to the external tools configuration file

    Returns:
    --------
    dict : A dictionary containing the external tools configuration
    """
    config_filepath = Path(external_config_filepath)
    assert config_filepath.exists(), "config_filepath does not exist: {}".format(config_filepath)
    with open(config_filepath) as json_file:
        external_tools = json.load(json_file)
    return external_tools
