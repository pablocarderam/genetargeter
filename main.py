"""
Handles connection with client, calls main method to build targeted vector with
arguments from the client. Websockets implementation from
https://github.com/BruceEckel/hello-flask-websockets
"""

from GeneTargeterMethods import *;

import os
from gevent import monkey
monkey.patch_all()

import time
from flask import Flask, render_template, session, request
from flask.ext.socketio import SocketIO, emit, disconnect

app = Flask(__name__)
app.debug = True
app.config['SECRET_KEY'] = 'secret!'
socketio = SocketIO(app)

@app.route('/')
def index():
    return render_template("index.html")


@socketio.on('sendGeneFile', namespace='/link')
def test_message(message):
    # process input message into geneName, geneFileStr, HRann and other params
    #output = pSN054TargetGene_StrMeth(geneName, geneFileStr, HRannotated=HRann); # call result
    sendMsg('Process complete');
    print message['data'] + " :: received"
    """sendMsg(output[0]);
    sendOligos(output[1]);
    sendGene(output[2]);
    sendPlasmid(output[3]);"""


def sendMsg(msg):
    socketio.emit('message',
                  {'data': msg},
                  namespace='/link');


@socketio.on('disconnect request', namespace='/link')
def disconnect_request():
    emit('my response',
         {'data': 'Disconnected!', 'count': session['receive_count']})
    disconnect()


@socketio.on('connect', namespace='/link')
def test_connect():
    emit('my response', {'data': 'Connected', 'count': 0})


@socketio.on('disconnect', namespace='/link')
def test_disconnect():
    print('Client disconnected')


if __name__ == "__main__":
    # Fetch the environment variable (so it works on Heroku):
    socketio.run(app, host='0.0.0.0', port=int(os.environ.get("PORT", 5000)))
