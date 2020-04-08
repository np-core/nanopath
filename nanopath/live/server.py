import logging
import argparse

from pathlib import Path
from nanopath.pipelines import NanoPathLive
from flask import Flask
from flask_socketio import SocketIO, emit

parser = argparse.ArgumentParser(
    description='Server script for NanoPath Live v0.1.0'
)
parser.add_argument(
    '--directory',
    '-d',
    type=Path,
    default=Path('/data/nanopath/live'),
    help='Path to directory containing output from pipeline [/data/nanopath/live]'
)
parser.add_argument(
    '--port',
    '-p',
    type=int,
    default=5000,
    help='Port for sockets connecting to dashboard [5000]'
)
parser.add_argument(
    '--log',
    '-l',
    type=Path,
    default=None,
    help='Path to log file output [STDOUT]'
)
args = parser.parse_args()

DEBUG = True

app = Flask(__name__)
app.config.from_object(__name__)

socketio = SocketIO(app, host='0.0.0.0', port=args.port, cors_allowed_origins="*")  # TODO check security here

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s]  %(message)s",
    datefmt='%H:%M:%S',
    filename=args.log
)

log = logging.getLogger(__name__)

# SocketIO response functions

log.info('Server log started, logging operations.')

npl = NanoPathLive(path=args.directory)


@socketio.on('settings_ping_server')
def settings_ping_server():
    log.info('Server ping received from client.')
    emit(
        'settings_ping_server', {'data': 'Server ping received.'}
    )
    log.info('Server pong emitted to client.')

# from: UserReports


@socketio.on('get_live_update')
def get_live_update(get_live_update):

    server_data = get_live_update
    log.info('NanoPath live dashboard contacted server.')
    log.info(f'NanoPath will parse the pipeline files at {npl.path}')

    read_distribution_params: dict = \
        server_data['user_settings'].get('read_distributions')

    stats_data = npl.update_stats_view()
    telemetry_data, _ = npl.update_telemetry_view()
    barcode_data, read_lengths, read_qualities = npl.update_run_view(
        **read_distribution_params
    )

    emit(
        'get_live_update', {
            'data': {
                'stats_view': stats_data,
                'run_view': {
                    'barcodes': barcode_data,
                    'read_lengths': read_lengths,
                    'read_qualities': read_qualities
                }
            }
        }
    )


if __name__ == "__main__":
    socketio.run(app)
