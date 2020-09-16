import logging
import shutil
import hashlib
import yaml

from datetime import date
from pathlib import Path
from nanopath.netflow import NetflowClient
from nanopath.utils import PoreLogger


def blake(file: Path):
    blake_hash = hashlib.blake2b()
    with file.open("rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            blake_hash.update(chunk)
    return blake_hash.hexdigest()


class Dataset(PoreLogger):

    def __init__(self, name: str, collection: str, version: int, author: str, data: Path):

        PoreLogger.__init__(self, level=logging.INFO, name="Dataset")

        self.name = name
        self.collection = collection
        self.version = version
        self.data = data
        self.author = author

        self.dataset = Path(f"{self.collection}-{self.name}-v{self.version}")

        self.files = None

    def create_template_directory(self):

        """ Create a copy of the data file directory in a named directory in the current workdir """

        self.logger.info(f'Create dataset: {self.dataset.name}')

        try:
            self.dataset.mkdir(parents=True, exist_ok=False)
        except FileExistsError:
            self.logger.info(f'Could not create the existing directory: {self.dataset.name}')
            raise

        # Copy data directory into template dataset directory and idnex the files
        dataset_files = {}
        for f in self.data.glob("*"):  # only ever get files (not directories) at the base level
            dst = self.dataset/f.name

            self.logger.info(f"Found dataset file: {f.name}")
            self.logger.info(f"Copy to template directory: {f.name}")
            dataset_files[f.name] = self.dataset / f.name
            shutil.copy(src=str(f), dst=str(dst))

        self.files = dataset_files

    def create_metadata(self):

        base = {
            'name': self.name,
            'collection': self.collection,
            'version': self.version,
            'date': date.today().strftime("%B %d, %Y"),
            'author': self.author,
            'description': '',
            'data': []
        }

        file = {'name': '', 'version': 1, 'description': ''}

        for fname, fpath in self.files.items():
            fmeta = file.copy()
            fmeta['name'] = fname
            base['data'].append(fmeta)

        meta = self.dataset / 'dataset.yaml'

        with meta.open("w") as mout:
            yaml.dump(base, mout, sort_keys=False)


class DataStorage(PoreLogger):

    def __init__(self, data: Path or None, config: Path):

        PoreLogger.__init__(self, level=logging.INFO, name='DataStorage')

        self.data = data
        self.workdir = Path.home() / '.npc' / '.storage_client'

        self.netflow = NetflowClient(
            node_config=config, workdir=Path.home() / '.npc' / '.storage_client'
        )
        self.netflow.configure_storage_nodes(client_entry='storage')

    def list(self):

        for node, node_client in self.netflow.clients.items():
            node_client.list_storage()  # connects and disconnects automatically for each client

    def md5sum(self):

        pass

    def containerize(self):

        pass

    def archive(self):

        pass

    def clean(self):

        shutil.rmtree(self.workdir)
        self.logger.info(
            f'Removed working directory of storage client: {self.workdir}'
        )









