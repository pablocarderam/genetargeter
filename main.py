"""
Handles connection with client, calls main method to build targeted vector with
arguments from the client. Websockets implementation from
https://github.com/BruceEckel/hello-flask-websockets
"""

from py.GeneTargeterMethods import *; # main python library in py folder

# Imports from Websockets tutorial:
import os
from gevent import monkey
monkey.patch_all()

import time
from flask import Flask, render_template, session, request
from flask.ext.socketio import SocketIO, emit, disconnect

app = Flask(__name__)
app.debug = True # change in dev/prod
app.config['SECRET_KEY'] = 'secret!'
socketio = SocketIO(app)

sep = ":::"; # separates files

pk = open("pk"); # Access given file
passcode = pk.read().split(sep)[1]; # passcode
pk.close();

@app.route('/')
def index():
    return render_template("index.html");


@socketio.on('sendGeneFile', namespace='/link')
def gene_message(message):
    # process input message into geneName, geneFileStr, HRann and other params
    msgList = message["data"].split(sep);
    queryNumber = msgList[0]; # number identifier
    geneFileStr = msgList[1];
    geneName = msgList[2];
    HRann = (msgList[3] == "true");
    lengthLHR = [int(i) for i in msgList[4].split(",")];
    lengthRHR = [int(i) for i in msgList[5].split(",")];
    lengthGib = [int(i) for i in msgList[6].split(",")];
    optimLHR = [int(i) for i in msgList[7].split(",")];
    optimRHR = [int(i) for i in msgList[8].split(",")];
    endsLHR = int(msgList[9]);
    endsRHR = int(msgList[10]);
    endTempLHR = float(msgList[11]);
    endTempRHR = float(msgList[12]);
    gibTemp = float(msgList[13]);
    gibTDif = float(msgList[14]);
    maxDistLHR = int(msgList[15]);
    maxDistRHR = int(msgList[16]);
    minFragSize = int(msgList[17]);
    optimOrg = msgList[18];
    codonSampling = (msgList[19] == "Codon Sampling");
    minGRNAGCContent = float(msgList[20])/100.0;
    minOnTargetScore = float(msgList[21]);
    offTargetMethod = msgList[22];
    offTargetThreshold = float(msgList[23]);
    maxOffTargetHitScore = float(msgList[24]);

    output = pSN054TargetGene(geneName, geneFileStr, codonOptimize=optimOrg, useFileStrs=True, HRannotated=HRann,lengthLHR=lengthLHR, lengthRHR=lengthRHR, gibsonHomRange=lengthGib, optimRangeLHR=optimLHR, optimRangeRHR=optimRHR, endSizeLHR=endsLHR, endSizeRHR=endsRHR, endTempLHR=endTempLHR, endTempRHR=endTempRHR, gibTemp=gibTemp, gibTDif=gibTDif, maxDistLHR=maxDistLHR, maxDistRHR=maxDistRHR, minGBlockSize=minFragSize, codonSampling=codonSampling, minGRNAGCContent=minGRNAGCContent, minOnTargetScore=minOnTargetScore, offTargetMethod=offTargetMethod, offTargetThreshold=offTargetThreshold, maxOffTargetHitScore=maxOffTargetHitScore); # call result
    outMsg = queryNumber + sep + output["geneName"] + sep + output["geneFileStr"] + sep + output["plasmidFileStr"] + sep + output["editedLocusFileStr"] + sep + output["oligoFileStr"] + sep + output["logFileStr"];
    sendMsg('Process complete',"misc");
    sendMsg(outMsg, "geneOutput");


@socketio.on('misc', namespace='/link')
def misc_message(message):
    print message['data'] + " :: received"


def sendMsg(msg,pType):
    socketio.emit(pType, {'data': msg}, namespace='/link');
    print pType + " :: sent"


@socketio.on('disconnect request', namespace='/link')
def disconnect_request():
    emit('my response', {'data': 'Disconnected!', 'count': session['receive_count']} )
    disconnect()


@socketio.on('connect', namespace='/link')
def test_connect():
    emit('my response', {'data': 'Connected', 'count': 0})


@socketio.on('disconnect', namespace='/link')
def test_disconnect():
    print('Client disconnected')

@socketio.on('cred', namespace='/link')
def validate_credentials(message):
    print('Received credentials: ');
    if message["data"] == passcode:
        emit('validCred',{'data':'You know it!'});
    else:
        emit('invalidCred',{'data':"You either have it or you don't"});


if __name__ == "__main__":
    # Fetch the environment variable (so it works on Heroku):
    socketio.run(app, host='0.0.0.0', port=int(os.environ.get("PORT", 5000)))
