import configparser
from pathlib import Path
import os, sys
sys.path.append(os.path.abspath('..'))

CONFIG_FILE = 'config.ini'


class ConfigManagerClass:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(ConfigManagerClass, cls).__new__(cls)
            cls._instance._config = configparser.ConfigParser()
            cls._instance.config_file = os.path.abspath(CONFIG_FILE)
            cls._instance.load_config()
        return cls._instance

    def load_config(self):
        dir_cur = Path.cwd()
        config_file_full = dir_cur / CONFIG_FILE
        if not config_file_full.exists():   # Depends on where the module is run from, from main (where the config file is located) or one of the subfolders
            config_file_full = dir_cur.parent / CONFIG_FILE
        if config_file_full.exists():
            self._config.read(self.config_file)
            print(f'Config file loaded: {self.config_file}')
        else:
            print(f'Config file does not exist: {CONFIG_FILE}')

    def get(self, section, option, fallback=None):
        return self._config.get(section, option, fallback=fallback)

    def getint(self, section, option, fallback=None):
        return self._config.getint(section, option, fallback=fallback)

    def getfloat(self, section, option, fallback=None):
        return self._config.getfloat(section, option, fallback=fallback)

    def getboolean(self, section, option, fallback=None):
        return self._config.getboolean(section, option, fallback=fallback)

    def set(self, section, key, value):
        if not self._config.has_section(section):
            self._config.add_section(section)
        self._config.set(section, key, str(value))
        self.save_config()

    def save_config(self):
        with open(self.config_file, 'w') as f:
            self._config.write(f)


# Singleton instance
config = ConfigManagerClass()
