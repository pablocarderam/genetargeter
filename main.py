#!/usr/bin/env python2

"""
Handles connection with client, calls main method to build targeted vector with
arguments from the client. Websockets implementation from
https://github.com/BruceEckel/hello-flask-websockets
"""

from py.genetargeter.constants import *; # main python library in py folder
from py.genetargeter.GeneTargeterMethods import *; # main python library in py folder
from py.genetargeter.inputProcessing import *; # input handling

# Imports from Websockets tutorial:
import os
import eventlet
eventlet.monkey_patch()

import time
from flask import Flask, render_template, session, request
from flask_socketio import SocketIO, emit, disconnect, join_room

import traceback # error handling

app = Flask(__name__)
app.debug = False #True # change in dev/prod
app.config['SECRET_KEY'] = 'secret!'
socketio = SocketIO(app)

sep = ":::"; # separates files

pk = open("pk"); # Access given file
passcode = pk.read().split(sep)[1]; # passcode
pk.close();

@app.route('/')
def index():
    return render_template("index.html");


@socketio.on('join', namespace='/link')
def createChannel(msg):
    print('Client joined: ' + msg['channel'])
    join_room(msg['channel'])

@socketio.on('sendGeneFile', namespace='/link')
def gene_message(message):
    # print(message)
    # process input message into geneName, geneFileStr, HRann and other params
    channel_id = message['channel'] # get this channel id
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
    codonOptimize = msgList[18];
    codonSampling = (msgList[19] == "Codon Sampling");
    minGRNAGCContent = float(msgList[20])/100.0;
    onTargetMethod = msgList[21];
    minOnTargetScore = float(msgList[22]);
    offTargetMethod = msgList[23];
    minOffTargetScore = float(msgList[24]);
    maxOffTargetHitScore = float(msgList[25]);
    enzyme = msgList[26];
    PAM = msgList[27];
    gBlockDefault = (msgList[28] == "true");
    plasmidType = msgList[29];
    haTagMsg = msgList[30];
    setCoding = msgList[31];

    geneGBs = preprocessInputFile(geneName, geneFileStr, useFileStrs=True, doingAltSplices=True); # obtain gb file(s) to be processed
    outMsg = queryNumber; # will store output message
    outMsgHA = "HA_Tag_Design-" + str(queryNumber); # will store output message for HA tag design, if any
    for gbName in geneGBs: # for every gb
        sigPep = chkSignalPeptide5Prime(gbName,signalPDB); # signalP info
        haTag = False; # don't use HA tags by default
        if plasmidType == "pSN150" and ( haTagMsg == "Yes" or (haTagMsg == "Auto" and not sigPep) ): # if forcing HA tags or 5' end does not contain signal peptide and auto HA-tagging
            haTag = True; # use HA tags

        try:
            output = targetGene(gbName, geneGBs[gbName], codonOptimize=codonOptimize, useFileStrs=True, HRannotated=HRann,lengthLHR=lengthLHR, lengthRHR=lengthRHR, gibsonHomRange=lengthGib, optimRangeLHR=optimLHR, optimRangeRHR=optimRHR, endSizeLHR=endsLHR, endSizeRHR=endsRHR, endTempLHR=endTempLHR, endTempRHR=endTempRHR, gibTemp=gibTemp, gibTDif=gibTDif, maxDistLHR=maxDistLHR, maxDistRHR=maxDistRHR, minGBlockSize=minFragSize, codonSampling=codonSampling, minGRNAGCContent=minGRNAGCContent, onTargetMethod=onTargetMethod, minOnTargetScore=minOnTargetScore, offTargetMethod=offTargetMethod, minOffTargetScore=minOffTargetScore, maxOffTargetHitScore=maxOffTargetHitScore, enzyme=enzyme, PAM=PAM, gBlockDefault=gBlockDefault, plasmidType=plasmidType, haTag=haTag, sigPep=sigPep, setCoding=setCoding, filterCutSites=cut_sites[plasmidType]); # call result

            outMsg = outMsg + sep + output["geneName"] + sep + output["geneFileStr"] + sep + output["plasmidFileStr"] + sep + output["editedLocusFileStr"] + sep + output["oligoFileStr"] + sep + output["gBlockFileStr"] + sep + output["gRNATable"] + sep + output["logFileStr"];
            if haTag and plasmidType == "pSN150": # if using HA tags and pSN150,
                outputHA = output["outputHA"]; # save HA outputs
                outMsgHA = outMsgHA + sep + outputHA["geneName"] + sep + outputHA["geneFileStr"] + sep + outputHA["plasmidFileStr"] + sep + outputHA["editedLocusFileStr"] + sep + outputHA["oligoFileStr"] + sep + outputHA["gRNATable"] + sep + outputHA["gBlockFileStr"] + sep + outputHA["logFileStr"];
                sendMsg(outMsgHA, "geneOutput", channel_id);

        except Exception, e:
            sendMsg(traceback.format_exc(), "runError", channel_id);
            raise e

    sendMsg('Process complete',"misc",channel_id);
    sendMsg(outMsg, "geneOutput",channel_id);


@socketio.on('misc', namespace='/link')
def misc_message(message):
    print message['data'] + " :: received"


def sendMsg(msg,pType,channel_id):
    socketio.emit(pType, {'data': msg}, namespace='/link', room=channel_id);
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


@socketio.on('geneNames', namespace='/link')
def check_gene_files(message):
    print('Received request for genes');
    geneIDs = message["data"].split('\n')
    genesNotFound = ''
    geneFileStr = ''
    for geneID in geneIDs:
        geneID = geneID.strip()
        if len(geneID) > 0:
            if os.path.isfile("input/genomes/geneFiles/"+geneID+".gb"):
                file = open("input/genomes/geneFiles/"+geneID+".gb", "r")
                geneFileStr = geneFileStr + file.read() + sep
                file.close()
            else:
                genesNotFound = genesNotFound + geneID + sep


    if len(genesNotFound) > 0:
        emit('genesFilesNotFound',{'data':genesNotFound[0:-len(sep)]})
    elif len(geneFileStr) > 0:
        emit('genesFilesFound',{'data':geneFileStr[0:-len(sep)]})


if __name__ == "__main__":
    # Fetch the environment variable (so it works on Heroku):
    socketio.run(app, host='0.0.0.0', port=int(os.environ.get("PORT", 5000)))
