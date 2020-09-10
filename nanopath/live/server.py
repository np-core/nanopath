
import logging
import argparse

from pathlib import Path
from nanopath.server import AppServer

from flask import Flask
from flask_socketio import SocketIO, emit

parser = argparse.ArgumentParser(
    description='Server terminal for NanoPath Live v0.1.0'
)
parser.add_argument(
    '--port',
    '-p',
    type=int,
    default=5000,
    help='Port for sockets connecting to dashboard [5000]'
)
parser.add_argument(
    '--pathogen',
    type=Path,
    default=Path('/data/nanopath/cadhla_working/nexus/fastq/results/kraken'),
    help='Location of pathogen pipeline output to parse [$CWD]'
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

socketio = SocketIO(
    app, host='0.0.0.0', port=args.port, cors_allowed_origins="*"
)  # TODO check security here

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s]  %(message)s",
    datefmt='%H:%M:%S',
    filename=args.log
)

log = logging.getLogger(__name__)

# SocketIO response functions

log.info('Server log started, logging operations.')

app_server = AppServer(pathogen_path=args.pathogen)


@socketio.on('ping_server')
def settings_ping_server():
    log.info('Server ping received from client.')
    emit(
        'ping_server', {'data': 'Server ping received.'}
    )
    log.info('Server pong emitted to client.')

# from: UserReports


@socketio.on('get_pathogen_data')
def get_pathogen_data(settings):

    log.info('NanoPath live dashboard contacted server to get pathogen data.')

    print(settings)

    _, group_reads, db_server = app_server.collect_pathogen_results(
        minor_threshold=float(settings['minor_threshold']),
        major_threshold=float(settings['major_threshold']),
        groupby=fr"{settings['group_regex']}",

    )

    # print(db_server)

    emit('get_pathogen_data', {
        'message': "Hello from Server",
        'server_data': {'data': db_server, 'reads': group_reads}
    })


@socketio.on('create_report')
def create_report(app_data):

    log.info('NanoPath live dashboard contacted server to create report.')

    db_reports, group_reads, _ = app_server.collect_pathogen_results(
        groupby=fr"{app_data['setting']['group_regex']}"  # No thresholds imposed at collection
    )
    app_server.collect_pathogen_report(  # reporting threshold in app_data
        app_data=app_data,
        db_reports=db_reports,
        group_reads=group_reads,
        outdir=Path.home() / 'Desktop'
    )

    emit('create_report', {
        'message': "Created report on server",
        'status': 'success',
        'report_path': ''
    })


if __name__ == "__main__":
    socketio.run(app)