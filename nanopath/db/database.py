from pathlib import Path


class MongoAPI:
    """ Convenience base class for utilities around PyMongo and MongoDB """

    def __init__(
        self,
        config: Path = Path.home() / '.ergc' / 'db' / 'config.json'
    ):

        pass

    def is_connected(self):

        pass