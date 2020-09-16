""" All things related to running the server and dashboard online """

import os
import uuid
import time
import shutil
import subprocess
import yaml
import logging
from distutils.dir_util import copy_tree

from datetime import datetime
from paramiko import SSHClient, AutoAddPolicy
from paramiko.ssh_exception import AuthenticationException, SSHException

from colorama import Fore

from nanopath.utils import PoreLogger
from pathlib import Path, PosixPath

R = Fore.RED
Y = Fore.YELLOW
C = Fore.CYAN
G = Fore.GREEN
M = Fore.MAGENTA
RE = Fore.RESET


class ClientError(Exception):
    pass


class NetflowClient(PoreLogger):

    def __init__(self, node_config: Path, workdir: Path):

        PoreLogger.__init__(self, level=logging.INFO, name='NanoPathClient')

        logging.getLogger("paramiko").setLevel(logging.WARNING)

        self.logger.info(
            f"{Y}NanoPathClient{RE}"
        )

        # Read the netflow configurations
        self.config = self.read_config(file=node_config)

        # Working directory structure
        self.workdir = workdir
        self.nexus = self.workdir / 'nexus'
        self.nexus.mkdir(parents=True, exist_ok=True)

        # Store the clients for later use:
        self.clients = {}

    def configure_storage_nodes(self, client_entry: str = 'storage'):

        self.logger.info(f'Start storage node configuration')
        for node, config in self.config[client_entry].items():
            ssh_key = config['ssh_key']
            key_path = os.path.expanduser(ssh_key)

            client = StorageClient(
                name=node,
                host=str(config['host']),
                user=str(config['user']),
                ssh_key=os.path.abspath(key_path),
                storage_path=str(config['storage_path'])
            )

            # Storage setup and execution check
            client.check_storage_node()

            self.clients[node] = client

            self.logger.info(
                f'{Y}Node {C}{node}{Y} configured ({C}{client.host}{Y}){RE}'
            )

        self.logger.info(
            f'Storage node configuration complete'
        )

    def configure_pipeline_nodes(self, client_entry: str = 'nodes'):

        """ Checks for server configurations and connections """

        self.logger.info(f'Starting pipeline node configuration')
        for node, config in self.config[client_entry].items():
            ssh_key = config['ssh_key']
            key_path = os.path.expanduser(ssh_key)

            client = PipelineClient(
                name=node,
                nexus=self.nexus / node,  # path to result and data hub on main
                host=str(config['host']),
                user=str(config['user']),
                nextflow=PosixPath(config['nextflow']),
                container=str(config['container']),
                pipeline=str(config['pipeline']),
                profile=str(config['profile']),
                workdir=PosixPath(config['workdir']),
                ssh_key=os.path.abspath(key_path),
                params=config['params']
            )

            self.logger.info(
                f'{Y}Checks started:  {C}{node}{RE}'
            )
            # Remote server connection check
            client.test_connection()
            # Nextflow execution and version check
            client.check_nextflow_version()
            # Pipeline setup and execution check
            client.check_pipeline_node(server=node)

            self.logger.info(
                f'{Y}Checks complete for node: {C}{node}{RE}'
            )

            # Create the corresponding result data hub on main host
            (self.nexus / node / 'results').mkdir(parents=True, exist_ok=True)
            self.logger.info(
                f'{Y}Client node data nexus on local host: {G}{self.nexus / node}{RE}'
            )

            self.clients[node] = client

            self.logger.info(
                f'{G}Node {C}{node}{RE} {G}configured{RE} ({C}{client.host}{RE})'
            )

        self.logger.info(
            f'Pipeline node configuration complete'
        )

    def launch_pipelines(
        self,
        fastq: Path = None,
        illumina: bool = False,
        nuke: bool = False,
        keep_open: bool = True,
        resume: bool = False,
        recursive: bool = False
    ):

        """ Launch pipelines on configured servers - completed runs

        At the moment, the framework is restricted to one Fast5 and one
        linked Fastq pipeline, that can be outsourced to a high-resource
        or GPU server. In the future, the framework should allow arbitrary
        coordination of pipelines across servers and in the cloud. However,
        for now, operating the real-time pathogen classifiers is priority
        for the sepsis work.

        :param fastq: path to directory containing fastq

        """

        if illumina:
            fastq_extensions = ('.fq.gz', '.fastq.gz')
        else:
            fastq_extensions = ('.fq', '.fastq')

        if fastq is not None:
            # Launch reads pipeline, check if pipeline clients are configured:
            try:
                fastq_client = self.clients['fastq']
            except KeyError:
                self.logger.error(
                    f"{R}Could not detect configured netflow, did you call "
                    f"configure_server_clients() first?{RE}"
                )
                raise

            if recursive:
                glob_pattern = '**/*'
            else:
                glob_pattern = '*'

            fastq_files = [
                p for ext in fastq_extensions
                for p in fastq.glob(glob_pattern+ext)
            ]

            # Upload or copy all files to target directory on local or remote
            fastq_client.upload_files(
                files=fastq_files, remote_path=fastq_client.file_input, symlink=True
            )

            # Start pipeline in its own screen on server
            pipeline_screen = fastq_client.run_pipeline(
                screen=True, keep_open=keep_open, resume=resume
            )  # If keep open, keep going below will always be True

            # Instantiate the result polling for the pipeline
            fastq_pollster = PipelinePollster(pipeline_client=fastq_client)

            try:
                self.logger.info(
                    f'{Y}Press {M}Ctrl+C{RE} {Y}to {M}stop pipelines {Y}and {M}exit{RE}'
                )
                keep_going = True
                while keep_going:
                    # Polling pipeline results
                    keep_going = fastq_pollster.poll_results(
                        pipeline_screen=pipeline_screen
                    )
                    time.sleep(3)
            except KeyboardInterrupt:
                self.clean(nuke=nuke)
                exit(0)

    def clean(self, nuke: bool = False):

        """ Remove active screens and optionally all traces of use """

        print()  # Ctrl+C offset ugly
        for server, client in self.clients.items():
            self.logger.info(
                f'{Y}Stopping operations on netflow{RE}: '
                f'{C}{server}{RE} ({C}{client.host}{RE})'
            )
            for name in client.active_screens.keys():
                client.kill_screen(grep=name)
            if nuke:
                client.nuke()

        if nuke:
            shutil.rmtree(self.workdir)

        self.logger.info(
            f'{RE}Server and netflow clean-up complete '
            f'{f"({R}NUKED{RE})" if nuke else ""}{RE}'
        )

    def read_config(self, file: Path):

        """ Read server and pipeline configuration YAML """

        # Load the configuration YAML:

        with file.open() as fin:
            config = yaml.load(fin, Loader=yaml.BaseLoader)

        # Check the server configuration

        try:
            server_config = config['nodes']
        except KeyError:
            self.logger.error(
                "Could not detect required entry 'nodes' in configuration"
            )
            raise

        # TODO: check pipeline parameters

        return config

    def check_pipeline_params(self, config: dict):

        """ Make sure that param entries in confuration file are complete """

        return


class Client(PoreLogger):

    def __init__(self, host, user, ssh_key, name, logger_name: str = 'PipelineClient'):

        PoreLogger.__init__(self, level=logging.INFO, name=logger_name)

        self.remote = True if host != 'local' else False
        self.name = name

        self.host = host
        self.user = user
        self.ssh_key = ssh_key

        self.client = None
        self.sftp = None
        self.conn = None

        self.active_screens = {}

    def _connect(self):

        """Open connection to remote host."""

        if not self.conn:
            self.logger.debug('Open connection to remote client')
            try:
                self.client = SSHClient()
                self.client.load_system_host_keys()
                self.client.set_missing_host_key_policy(AutoAddPolicy())
                self.client.connect(
                    self.host,
                    username=self.user,
                    key_filename=self.ssh_key,
                    look_for_keys=True,
                    timeout=5000
                )
                self.sftp = self.client.open_sftp()
            except AuthenticationException as e:
                self.logger.error('Authentication failed: ' + str(e))
                self.logger.info('SFTP Authentication failed')
                return
            except SSHException as e:
                self.logger.error('Unable to establish SSH connection: ' + str(e))
                if 'Connection timed out' in str(e):
                    self.logger.info('SFTP timed out connecting to ' + self.host)
                else:
                    self.logger.info('SFTP unable to establish SSH connection')
                return
            except Exception as e:
                self.logger.error('Exception signing in to SFTP server: ' + str(e))
                self.logger.info('SFTP exception signing in')
                return

        return True

    def _upload_single_file(
        self,
        file: str,
        remote_path: str,
        symlink: bool = True,
    ) -> str or None:

        """ Upload a single file to a remote directory or copy to local """

        file = os.path.realpath(file)  # always local, but resolve symlinks
        fname = os.path.basename(file)

        if self.remote:
            if not self.conn:
                self.logger.error(
                    f'{R}Could not upload file: no connection{RE}'
                )
                return
            try:
                self.sftp.put(file, remote_path=remote_path)
            except SSHException:  # TODO: is this exception correct for SFTP?
                self.logger.error(
                    f'{R}Error uploading file {Y}{fname}{G} to {C}{self.host}{RE}'
                )
                exit(1)
            finally:
                self.logger.info(
                    f'{Y}Uploaded file: {G}{fname} {Y}to {C}{self.host}{RE}:{G}{remote_path}{RE}'
                )
                return file
        else:
            try:
                if symlink:
                    # symlink as file
                    remote_path = str(PosixPath(remote_path) / fname)
                    os.symlink(file, remote_path)
                else:
                    # target is directory no matter
                    shutil.copy(file, remote_path)
            except FileExistsError:
                pass  # It's ok little one
            except IOError:
                self.logger.error(
                    f'{R}Error copying file: {fname}{R} to {C}{self.host}{RE}'
                )
                exit(1)
            finally:
                self.logger.info(
                    f'{Y}Copied file {G}{fname}{Y} to {C}{self.host}{RE}'
                )
                return file

    def _download_single_file(
        self, file: str, local_path: str, symlink: bool = True
    ) -> str or None:

        fname = PosixPath(file).name

        """ Download a single file to a remote directory or copy to local """
        if self.remote:
            if not self.conn:
                self.logger.error(
                    f'{R}Could not download file: no connection{RE}'
                )
                return
            try:
                self.sftp.get(file, local_path)
            except SSHException:  # TODO: is this exception correct for SFTP?
                self.logger.error(
                    f'{R}Error downloading file: {Y}{fname}{RE}'
                )
                exit(1)
            finally:
                self.logger.info(
                    f'{Y}Downloaded file: {G}{fname}{RE} {Y}to{RE} {G}{local_path}{RE}'
                )
                return file
        else:
            try:
                if symlink:
                    local_path = str(PosixPath(local_path) / fname)
                    os.symlink(file, local_path)
                else:
                    shutil.copy(file, local_path)
            except IOError:
                self.logger.error(
                    f'{R}Error copying file: {fname}{RE}'
                )
                exit(1)
            finally:
                self.logger.info(
                    f'{Y}Copied file: {G}{fname}{RE}'
                )
                return file

    def disconnect(self):

        """ Close SSH connection """

        if self.client:
            self.logger.debug('Close connection to remote client')
            self.conn = False
            self.client.close()
        if self.sftp:
            self.sftp.close()

    def execute_cmd(
        self,
        cmd,
        screen: bool = False,
        screen_id: str = "",
        keep_open: bool = False,
        confirmation: bool = False
    ) -> list:

        if screen:
            self.active_screens[screen_id] = {'cmd': cmd, 'date': datetime.now()}
            # Prepend command with the detached, named screen session and execute
            cmd = f"screen -S {screen_id} -dm bash " \
                f"-c '{cmd}{'; exec bash' if keep_open else ''}'"

        if confirmation:
            cmd += " && echo 1"

        if self.remote:
            self.conn = self._connect()
            stdin, stdout, stderr = self.client.exec_command(cmd)
            # stdout.channel.recv_exit_status()  # wait
            response = stdout.readlines()
            self.disconnect()
            return [b.strip() for b in response]
        else:
            # TODO see security concerns
            p = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, shell=True, executable='/bin/bash'
            )
            return [b.decode('utf-8').strip() for b in p.stdout.readlines()]

    def test_connection(self):

        if self.remote:
            try:
                self._connect()
            except SSHException as error:
                self.logger.error('Could not connect to remote Netflow')
                raise error

            self.disconnect()
        else:
            self.logger.info('Using local shell commands')
            return

    def upload_files(
        self,
        files: [Path] or Path,
        remote_path: Path,
        symlink: bool = False
    ):

        """ Simple wrapper to upload a list of files on a single connection """

        if isinstance(files, Path):
            files = [files]

        if self.remote:
            self.conn = self._connect()  # use only when remote

        try:
            for file in files:
                self._upload_single_file(
                    file=str(file), remote_path=str(remote_path), symlink=symlink
                )
        except KeyboardInterrupt:
            self.disconnect()
            raise
        except:  # Harsh, but if unexpected things go wrong SSH must be closed
            self.disconnect()  # works also when local
            raise
        finally:
            self.disconnect()

    def download_files(
        self,
        files: Path or [Path],
        local_path: Path or [Path],
        symlink: bool = False
    ):

        """ Simple wrapper to download a list of files on a single connection """

        if isinstance(local_path, list) and isinstance(files, list):
            if len(local_path) != len(files):
                raise ClientError(
                    'When provided with two lists, (file or directory) '
                    'number of entries must be the same'
                )
        if isinstance(files, Path):
            files = [files]
        if isinstance(local_path, Path):
            local_path = [local_path]

        if self.remote:
            self.conn = self._connect()
        try:
            for i, file in enumerate(files):
                self._download_single_file(
                    file=str(file), local_path=str(local_path[i]), symlink=symlink
                )
        except:  # Harsh, but if unexpected things go wrong SSH must be closed
            self.disconnect()
            raise
        finally:
            self.disconnect()

    def kill_screen(self, grep: str = "np_"):

        kill_command = f"screen -ls | grep '{grep}' | cut -d. -f1 | xargs kill"
        self.execute_cmd(cmd=kill_command)
        self.logger.debug(
            f"Removed screens with grep: {grep}"
        )


class StorageClient(Client):

    """ Less complex variant of the SSH client manager for Netflow adapted for data storage """

    def __init__(self, host, user, ssh_key, name, storage_path):

        Client.__init__(self, host, user, ssh_key, name, logger_name='StorageClient')

        self.storage_path = storage_path

    def check_storage_node(self):

        # Is the storage base directory present?

        out = self.execute_cmd(f'[ -d "{self.storage_path}" ] && echo "1"')

        if not out:
            self.logger.info(
                f'{R}Could not detect storage directory{RE}'
            )
            out = self.execute_cmd(f'mkdir -p {self.storage_path} && echo "1"')
            if not out:
                self.logger.info(
                    f'{R}Could not create working storage directory: {self.storage_path}{RE}'
                )
                exit(1)
            self.logger.info(
                f'{Y}Created storage directory: {G}{self.storage_path}{RE}'
            )
        else:
            self.logger.debug(f'Data storage directory found on server: {self.storage_path}')

    def list_storage(self):

        self.logger.info(f"Listing datasets on node: {self.name}")

        if not self.conn:
            self.conn = self._connect()

        results = self.execute_cmd(f'find {self.storage_path} -type d -maxdepth 1')

        if results:
            try:
                for storage_directory in results[1:]:
                    try:
                        name = storage_directory.split('/')[-1]
                        self.logger.info(f'{C}{name}{RE}')
                    except IndexError:
                        self.logger.info('Could not execute listing on storage node')
            except IndexError:
                self.logger.info('Could not execute listing on storage node')
        else:
            self.logger.info('Could not execute listing on sotrage node')

        self.disconnect()

    def inspect_dataset(self, dataset: str):

        self.logger.info(f"Inspect dataset {dataset} on node: {self.name}")

        if not self.conn:
            self.conn = self._connect()

        dataset_storage_path = f"{self.storage_path}/{dataset}"

        out = self.execute_cmd(f'[ -d "{dataset_storage_path}" ] && echo "1"')
        if not out:
            self.logger.info(
                f'{R}Could not detect dataset directory in storage at: {dataset_storage_path}{RE}'
            )
            exit(1)
        else:
            dataset_meta_file = f"{dataset_storage_path}/dataset.yaml"
            out1 = self.execute_cmd(f'[ -f "{dataset_meta_file}" ] && echo "1"')
            if not out1:
                self.logger.info(
                    f'{R}Could not detect dataset meta data in storage at: {dataset_meta_file}{RE}'
                )
                exit(1)
            else:
                out2 = self.execute_cmd(f'grep -e "^date:" -e "^description:" {dataset_meta_file}')
                out3 = self.execute_cmd(f'du -sh {dataset_storage_path}')
                try:
                    size = out3[0].split()[0]
                    date = out2[0].replace("date: ", "")
                    description = out2[1].replace("description: ", "")
                    self.logger.info(
                        f'{Y}name: {C}{dataset} {Y}date: {C}{date} '
                        f'{Y}size: {C}{size} {Y}description: {C}{description}{RE}'
                    )
                except IndexError:
                    self.logger.info(
                        f'{R}Could not get meta data from: {dataset_meta_file}{RE}'
                    )

        self.disconnect()

    def download_dataset_directory(self, dataset: str, outdir: Path):

        """ Mirror / download a dataset directory recursively from remote """

        remote_directory = f"{self.storage_path}/{dataset}"

        if self.remote:
            self.conn = self._connect()  # use only when remote

            if not self.conn:
                self.logger.error(
                    f'{R}Could not download directory: no connection{RE}'
                )
                return

            self.logger.info(f'Create local directory {outdir}')
            outdir.mkdir(parents=True, exist_ok=True)

            for fpath in self.sftp.listdir(path=remote_directory):
                local_file = os.path.join(str(outdir), fpath)
                try:
                    self.sftp.get(
                        remotepath=f"{remote_directory}/{fpath}", localpath=local_file
                    )
                except SSHException:  # TODO: is this exception correct for SFTP?
                    self.logger.error(
                        f'{R}Error downloading file {Y}{fpath}{G} to {C}{local_file}{RE}'
                    )
                    exit(1)
                finally:
                    self.logger.info(
                        f'{Y}Downloaded file {G}{fpath} {Y}to {G}{local_file}{RE}'
                    )

    def upload_dataset_directory(self, local_path: str):

        """ Mirror / upload a dataset directory recursively to remote """

        local_path = os.path.realpath(local_path)  # always local, but resolve symlinks
        local_name = os.path.basename(local_path)

        remote_directory = f"{self.storage_path}/{local_name}"

        if self.remote:
            self.conn = self._connect()  # use only when remote

            if not self.conn:
                self.logger.error(
                    f'{R}Could not upload directory: no connection{RE}'
                )
                return

            self.logger.info(f'Create remote storage directory {remote_directory}')
            try:
                self.sftp.mkdir(remote_directory)
            except IOError:
                self.logger.info(f'{R}Dataset is already stored at: {Y}{remote_directory}{RE}')
                self.logger.info('Please delete from storage before attempting to store again')
                exit(1)

            for root, directories, files in os.walk(local_path):
                # no subdirectories, only files
                for f in files:
                    try:
                        self.sftp.put(os.path.join(root, f), remotepath=f"{remote_directory}/{f}")
                    except SSHException:  # TODO: is this exception correct for SFTP?
                        self.logger.error(
                            f'{R}Error uploading file {Y}{f}{G} to {C}{self.host}{RE}'
                        )
                        exit(1)
                    finally:
                        self.logger.info(
                            f'{Y}Uploaded file: {G}{f} {Y}to {C}{self.host}{RE}:{G}{remote_directory}{RE}'
                        )
        else:
            try:
                copy_tree(local_path, remote_directory)
            except FileExistsError:
                raise
            except IOError:
                self.logger.error(
                    f'{R}Error copying directory tree {local_path}{R} to {C}{self.host}{RE}'
                )
                exit(1)
            finally:
                self.logger.info(
                    f'{Y}Copied directory tree {G}{local_path}{Y} to {C}{self.host}{RE}'
                )

        return local_path


class PipelineClient(Client):
    """Client to interact with a local or remote host via SSH & SCP for executing Nextflow pipelines"""

    def __init__(
        self,
        host: str,
        user: str,
        ssh_key: Path,
        nexus: Path,
        name: str,
        nextflow: PosixPath,  # client node nextflow requirement
        pipeline: str,
        profile: str,
        container: str,
        workdir: PosixPath,
        params: dict,
    ):
        Client.__init__(self, host, user, ssh_key, name, logger_name='PipelineClient')

        self.nexus = nexus  # Path to hub for this pipeline on head node
        self.workdir = workdir  # Working directory for pipeline on client

        self.nextflow_version = 20.07
        self.basedir = "$HOME/.nanopath/pipelines"  # on client

        self.nextflow = nextflow
        self.container = container
        self.pipeline = pipeline
        self.profile = profile

        self.params = params

        self.file_input = self.workdir / 'files'
        self.file_output = self.workdir / 'results'

        self.nexus_results = self.nexus / 'results'  # mirrored hub results

        # Keep a registry of active screens and commands executed
        self.active_screens = {}

    def check_nextflow_version(self) -> float:

        out = self.execute_cmd(f"{self.nextflow} -v")

        if not out:
            self.logger.error(
                f'{R}Could not detect Nextflow - '
                f'did you setup the server correctly?{RE}'
            )
            exit(1)

        try:
            response = out[0]
            version_str = response.split()[2]
            major_version = ".".join(version_str.split('.')[:2])
            version = float(major_version)
        except (KeyError, ValueError):
            self.logger.error('Could not extract Nextflow version')
            raise

        if version < self.nextflow_version:
            raise ClientError(
                f'Nextflow version must be >= {self.nextflow_version}'
            )
        else:
            self.logger.info(f'Nextflow version is: {G}v{version}{RE}')

        return version

    def check_pipeline_node(self, server: str):

        # Are screens working?

        out = self.execute_cmd(f'screen -v && echo "1"')

        if not out:
            self.logger.info(
                f'{R}Could not detect screen utility on server{RE} - '
                f'did you run the server setup on {R}{server}{RE}?'
            )
            exit(1)
        else:
            self.logger.info('Screens are working on server')

        # Is the pipeline base directory present?

        out = self.execute_cmd(f'[ -d "{self.basedir}" ] && echo "1"')

        if not out:
            self.logger.info(
                f'{R}Could not detect pipeline directory{RE} - '
                f'did you run the server setup on {R}{server}{RE}?'
            )
            exit(1)

        # Is the selected pipeline present?

        pipeline_file = Path(self.basedir) / self.pipeline
        out = self.execute_cmd(f'[ -f "{pipeline_file}.nf" ] && echo "1"')

        if not out:
            self.logger.info(
                f'{R}Could not detect pipeline in directory{RE} - '
                f'did you run the server setup on {R}{server}{RE}?'
            )
            exit(1)

        self.logger.info(
            f'Pipeline: {G}{self.pipeline}.nf{RE}'
        )
        self.logger.info(
            f'Pipeline profile: {G}{self.profile}{RE}'
        )
        self.logger.info(
            f'Working directory: {G}{self.workdir}{RE}'
        )

        # Is the working directory present? If not list it

        out = self.execute_cmd(f'[ -d "{self.workdir}" ] && echo "1"')

        if not out:
            self.logger.info(
                f'Could not detect working directory: {self.workdir}'
            )
            self.logger.info(
                f'Create working directory: {self.workdir}'
            )
            out = self.execute_cmd(f'mkdir -p {self.workdir} && echo "1"')
            if not out:
                self.logger.info(
                    f'{R}Could not create working directory: {self.workdir}{RE}'
                )
                exit(1)

        # Is the file input directory present? If not list it

        out = self.execute_cmd(f'[ -d "{self.file_input}" ] && echo "1"')

        if not out:
            self.logger.info(
                f'Could not detect file input directory: {self.file_input}'
            )
            self.logger.info(
                f'Create file input directory: {self.file_input}'
            )
            out = self.execute_cmd(f'mkdir -p {self.file_input} && echo "1"')
            if not out:
                self.logger.info(
                    f'{R}Could not create file input directory: {self.file_input}{RE}'
                )
                exit(1)

    def run_pipeline(
        self, resume: bool = False, screen: bool = True, keep_open: bool = False
    ):

        """ Configure and run pipeline """

        pipeline_file = PosixPath(self.basedir) / (self.pipeline+".nf")
        resource_file = PosixPath(self.basedir) / 'nextflow.config'

        nextflow_command = f"{self.nextflow} run {pipeline_file} " \
            f"-c {resource_file} -profile {self.profile} " \
            f"-w {self.workdir / 'work'} "

        if resume:
            nextflow_command += f"-resume "

        nextflow_command += f"--files {self.file_input} " \
            f"--outdir {self.file_output} " \
            f"--container {self.container} "

        for param, value in self.params.items():
            # List in command line input as comma separated string,
            # which will be deconstructed in the workflow
            if isinstance(value, list):
                value = ",".join(value)

            # param values must be wrapped as string
            nextflow_command += f'--{param} "{value}" '

        print(nextflow_command)

        screen_id = f'np_{self.pipeline}_{str(uuid.uuid4())}'

        self.logger.info(
            f'{Y}Running pipeline on screen{RE}: '
            f'{G}{screen_id}{RE} ({C}{self.host}{RE})'
        )

        confirm = self.execute_cmd(
            nextflow_command, screen=screen, screen_id=screen_id,
            keep_open=keep_open, confirmation=True
        )

        try:
            ok = confirm[0]
        except KeyError:
            self.logger.error('Could not confirm execution of pipeline')
            raise

        if ok != '1':
            self.logger.error('Could not confirm execution of pipeline')
            exit(1)

        return screen_id

    def nuke(self):
        self.execute_cmd(cmd=f'rm -rf {self.workdir}', confirmation=True)
        self.logger.info(
            f"{Y}Removed working directory {G}{self.workdir}{RE} ({C}{self.host}{RE})"
        )


class PipelinePollster(PoreLogger):

    """ Checks on results and status of pipelines, poll with command remote or local """

    # TODO: implement threading

    def __init__(self, pipeline_client: PipelineClient):
        PoreLogger.__init__(self, level=logging.INFO, name='PipelinePollster')

        self.pipeline_client = pipeline_client

        # When instantiated, remove the previous registry if present
        self.result_registry = self.pipeline_client.workdir / 'result.registry'

        # List of files from last poll
        self.last_file_poll = []

        # If pipeline screen is queried:
        self.last_screen_open = False

    def poll_results(self, pipeline_screen: str):

        # Check for new result files:
        command = f"""
        if [[ -d "{self.pipeline_client.file_output}" ]]; then
            touch {self.result_registry} && \
            find {self.pipeline_client.file_output} -type f -not -path '*/\\.*'
        else
            echo -1
        fi
        """

        # Update the result file registry on netflow
        find_output = self.pipeline_client.execute_cmd(
            cmd=command, screen=False, confirmation=False
        )
        try:
            result = find_output[0]
            if result == "-1" and len(find_output) == 1:
                self.logger.info(f'{Y}Results have not been computed yet{RE}')
            else:
                diff = set(find_output).difference(self.last_file_poll)

                if len(diff) > 0:
                    self.logger.info(
                        f'[{C}{self.pipeline_client.name}{RE}] '
                        f'{Y}Detected {G}{len(diff)}{Y} new results from '
                        f'{G}{self.pipeline_client.pipeline}{RE} '
                        f'({C}{self.pipeline_client.host}{RE})'
                    )

                    local_paths = [
                        PosixPath(path_str.replace(
                            str(self.pipeline_client.file_output),
                            str(self.pipeline_client.nexus_results)
                        )).absolute() for path_str in diff
                    ]  # Files mirrored to local result directory in Nexus

                    # Create the corresponding parent directories first
                    # if they do not exist yet

                    for p in local_paths:
                        p.parent.mkdir(parents=True, exist_ok=True)

                    self.pipeline_client.download_files(
                        files=[Path(p) for p in diff],
                        local_path=local_paths
                    )  # two lists will download to individual paths

                else:
                    self.logger.info(
                        f'[{C}{self.pipeline_client.name}{RE}] '
                        f'No new results detected '
                        f'({C}{self.pipeline_client.host}{RE})'
                    )

                self.last_file_poll = find_output
        except IndexError:
            self.logger.info(f'{R}Failed to poll results: {RE}{self.pipeline_client.file_output}')

        if pipeline_screen:
            # Also query for the existence of the pipeline screen
            pcmd = f'screen -list | grep "{pipeline_screen}"'

            pipeline_open = self.pipeline_client.execute_cmd(
                cmd=pcmd, screen=False, confirmation=True
            )

            if not pipeline_open:  # empty grep result
                self.logger.info(
                    f'{Y}Pipeline completed on{RE}: {G}{pipeline_screen}{RE}'
                )
                self.last_screen_open = False
                # Remove active screen from netflow, so in clean up on False return
                # it isn't attempted to be killed
                self.pipeline_client.active_screens.pop(pipeline_screen)
                return False

            if len(pipeline_open) != 2:
                # unique screen id queried, one line output, one line confirm
                self.logger.info(
                    f'{R}Failed to poll pipeline screen: '
                    f'{RE}{pipeline_screen}{RE}'
                )
            else:
                try:
                    result = pipeline_open[1]
                    if result == "1":
                        self.last_screen_open = True
                    else:
                        self.logger.info(
                            f'{R}Failed to poll pipeline screen: '
                            f'{RE}{pipeline_screen}{RE}'
                        )
                except IndexError:
                    self.logger.info(
                        f'{R}Failed to poll pipeline screen: '
                        f'{RE}{pipeline_screen}{RE}'
                    )

        return True
