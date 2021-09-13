/* JS Document */
/* Main js */

var backToTopShown = true; // used to control back-to-top button appearance and disappearance
var currentOutput = []; // saves current output array to download multiple times if desired
var nonCoding = ''; // tag if noncoding
var haTagged = ''; // tag if contains HA tags
var fileCounter = 0; // used in coloring output file lines
var numFilesUploaded = 0;
var disconnected = false; // true when connection is lost in middle of operations
var validCredentials = false;
var numBkgrs = 15;

var geneServerFiles = [];
var geneServerNames = [];

var bulkFile = false;
var bulkFileArray = ["","","","","","",""];

function init() {
    document.getElementById("backToTopBtn").style.opacity = 0;
    EPPZScrollTo.scrollVerticalToElementById('Title', 0); // scroll to start at init
    // document.body.style.background = "linear-gradient(rgba(150, 200, 50, 0.4),rgba(76, 0, 0, 0.45)), url('../static/assets/bkgr/bkgr_" + Math.ceil(Math.random()*numBkgrs).toString() + ".jpg') no-repeat center center fixed";
    document.body.style.backgroundSize = "cover";
    document.getElementById('geneIDSubmissionTxt').value = 'PF3D7_0321600\nPF3D7_0819200\nPF3D7_1456700';
}

function geneIDSubmissionInterface() {
    document.getElementById('geneIDSubmission').style.display = "block";
    if (document.getElementById("geneIDSubmissionContent").classList.contains("modal-out")) {
        document.getElementById("geneIDSubmissionContent").classList.remove("modal-out");
    }
    if (!document.getElementById("geneIDSubmissionContent").classList.contains("modal-in")) {
        document.getElementById("geneIDSubmissionContent").classList.add("modal-in");
    }
}

function submitGeneIDs() {
    var geneIDs = document.getElementById('geneIDSubmissionTxt').value.trim();
    sendMessageToServer(geneIDs,"geneNames");
}

function genesFound(geneFileStr) {
    geneServerFiles = geneFileStr;

    var geneIDs = document.getElementById('geneIDSubmissionTxt').value.trim().split('\n');
    geneServerNames = [];
    for (var i = 0; i < geneIDs.length; i++) {
        if (geneIDs[i].trim().length > 0) {
            geneServerNames.push(geneIDs[i].trim())
        }
    }

    if (geneServerNames.length == 1) {
      document.getElementById("run").innerHTML = "Target Gene!";
    }
    else {
      document.getElementById("run").innerHTML = "Target Genes!";
    }

    var parent = document.getElementById("selectedFiles");
    parent.innerHTML = "";
    var txt = "";
    fileCounter = 0;

    for (var i = 0; i < geneServerNames.length; i++) {
        var geneName = geneServerNames[i];
        if (geneName.length > 0) {
            txt = "";
            txt += (i+1) + ".";
            var p = document.createElement("span");
            txt += "  " + geneIDs[i] + '.gb | Gene name: ';
            var nameField = document.createElement("input");
            nameField.setAttribute("type","text");
            nameField.setAttribute("value",geneIDs[i]);
            p.innerHTML = txt;
            p.appendChild(nameField);
            p.innerHTML += "<span></span> <br>";
            parent.appendChild(p);
        }
    }

    closeModal('geneIDSubmissionContent','geneIDSubmission');
}

function genesNotFound(geneIDsStr) {
    var geneIDs = geneIDsStr;
    var txt = 'ERROR: Make sure these are <i>P. falciparum</i> 3D7 PlasmoDB gene IDs, or use the file upload feature on the main page to submit your own gene files. <br>Gene IDs not found:';
    for (var i = 0; i < geneIDs.length; i++) {
        txt += '<br>'+geneIDs[i];
    }

    document.getElementById('geneIDSubmissionError').innerHTML = txt
}

function askCredentials() {
    if (document.getElementById('run').innerHTML === "Target Gene!" || document.getElementById('run').innerHTML === "Target Genes!") {
        document.getElementById('myModal').style.display = "block";
        // document.getElementById("modFooter").innerHTML = "Type passcode";
        // document.getElementById("passcode").value = "";
        if (document.getElementById("modContent").classList.contains("modal-out")) {
            document.getElementById("modContent").classList.remove("modal-out");
        }
        if (!document.getElementById("modContent").classList.contains("modal-in")) {
            document.getElementById("modContent").classList.add("modal-in");
        }
    }
    else if (document.getElementById("run").innerHTML === "Select one or more GenBank gene files.") {
        // document.getElementById('geneFileForm').click();
    }
    else if (document.getElementById("run").innerHTML === "Completed!") {
        if (numFilesUploaded > 1) {
            document.getElementById("run").innerHTML = "Target Genes!";
        }
        else if (numFilesUploaded == 1) {
            document.getElementById("run").innerHTML = "Target Gene!";
        }
        else {
            document.getElementById("run").innerHTML = "Select one or more GenBank gene files.";
        }
    }
}

function verifyCredentials() {
    // sendMessageToServer(document.getElementById("passcode").value,"cred");
    validCred();
}

function invalidCred() {
    document.getElementById("modFooter").innerHTML = "Invalid passcode.";
}

function runError(msg) {
    document.getElementById('errorWindow').style.display = "block";
    document.getElementById("errorContent").innerHTML = msg;
    var formatMsg = String(msg).replace(/"/g,"'").replace(/\n/g,"%0D%0A");
    document.getElementById("modFooterError").innerHTML = '<h3><a href="mailto:pablocarderam@gmail.com?subject=GeneTargeter Error Report&body=Hola Pablo! \nI was trying to run GeneTargeter but ran into this error: %0D%0A%0D%0A'+formatMsg+'%0D%0A%0D%0AThanks!">Click to send error report to pablocarderam@gmail.com</a>.</h3>';
    if (document.getElementById("errorWindowContent").classList.contains("modal-out")) {
        document.getElementById("errorWindowContent").classList.remove("modal-out");
    }
    if (!document.getElementById("errorWindowContent").classList.contains("modal-in")) {
        document.getElementById("errorWindowContent").classList.add("modal-in");
    }
}

function validCred() {
    var x = document.getElementById("geneFileForm");
    validCredentials = true;
    fileCounter = 0;
    if ('files' in x) {
        if (x.files.length == 0 && geneServerFiles.length == 0) {
            document.getElementById("modFooter").innerHTML = "Select some gene files though.";
        }
        else {
            document.getElementById("modFooter").innerHTML = "...";
            closeModal('modContent','myModal');
            bulkFile = document.getElementById("bulkFileChkBox").checked;
            run();
        }
    }
    else {
        window.alert("Application error. Something went wrong because the gene file form is supposed to have a files property, but I can't find it. Sorry :P")
    }
}

function checkKey() {
    if (event.keyCode == 13) {
        verifyCredentials();
    }
}

// When the user clicks on <span> (x), close the modal
function closeModal(id,modal_id) {
    if (document.getElementById(id).classList.contains("modal-in")) {
        document.getElementById(id).classList.remove("modal-in");
    }
    if (!document.getElementById(id).classList.contains("modal-out")) {
        document.getElementById(id).classList.add("modal-out");
    }
    setTimeout(function() {displayNone(modal_id);}, 250);
}

function displayNone(id) {
    document.getElementById(id).style.display = "none";
}

// When the user clicks anywhere outside of the modal, close it
window.onclick = function(event) {
    if (event.target == document.getElementById('myModal')) {
        closeModal('modContent','myModal');
    }
    else if (event.target == document.getElementById('geneIDSubmission')) {
        closeModal('geneIDSubmissionContent','geneIDSubmission');
    }
}

function toggleOptions() {
    var growDiv = document.getElementById('moreOptions');
    if (document.getElementById('moreOptions').clientHeight) {
      document.getElementById('moreOptions').style.height = 0;
    } else {
      var wrapper = document.querySelector('.measuringWrapper');
      document.getElementById('moreOptions').style.height = document.querySelector('.measuringWrapper').clientHeight + "px";
    }

    if (document.getElementById("optionsBtn").innerHTML === "Show more options") {
        document.getElementById("optionsBtn").innerHTML = "Show less options";
    }
    else {
        document.getElementById("optionsBtn").innerHTML = "Show more options";
    }
}

// Load file adapted from http://www.w3schools.com/jsref/prop_fileupload_files.asp
function uploadPlasmidFile() {
    var x = document.getElementById("plasmidFileForm");
    var parent = document.getElementById("selectedPlasmidFile");
    var txt = "";
    fileCounter = 0;

    if ('files' in x) {
        if (x.files.length == 1) {
            txt = "";
            var file = x.files[0];
            var p = document.createElement("span");
            if (file.name.substr(file.name.length-3,file.name.length) === ".gb") {
                if ('name' in file) {
                    txt += "Base plasmid: " + file.name;
                }
                if ('size' in file) {
                    txt += " | " + Math.round(file.size/100)/10 + " KB";
                }
            }
            else if ('name' in file) {
                txt += "File " + file.name + " not a GenBank format file (.gb)";
            }
            p.innerHTML = txt;
            p.innerHTML += "<span></span> <br>";
            while (parent.lastChild) {
                parent.removeChild(parent.lastChild);
            }
            parent.appendChild(p);
            document.getElementById('plasmidType').value = "custom";
        }
        else {
             window.alert("Upload a single GenBank (.gb) file as the custom base plasmid.");
        }
    }
    else {
        if (x.value === "") {
            txt += "Select one or more files.";
        } else {
            txt += "Sorry, the files property is not supported by your browser.";
        }
        var p = document.createElement("p");
        p.innerHTML = txt;
        parent.appendChild(p);
    }
}

// Load file adapted from http://www.w3schools.com/jsref/prop_fileupload_files.asp
function uploadGeneFile(){
    var x = document.getElementById("geneFileForm");
    var parent = document.getElementById("selectedFiles");
    var txt = "";
    fileCounter = 0;
    while (parent.firstChild) {
        parent.removeChild(parent.firstChild);
    }
    if ('files' in x) {
        if (x.files.length == 0) {
            document.getElementById("run").innerHTML = "Select one or more GenBank gene files.";
        }
        else {
            if (x.files.length == 1) {
              document.getElementById("run").innerHTML = "Target Gene!";
            }
            else {
              document.getElementById("run").innerHTML = "Target Genes!";
            }

            for (var i = 0; i < x.files.length; i++) {
                txt = "";
                txt += (i+1) + ".";
                var file = x.files[i];
                var geneName = "";
                var p = document.createElement("span");
                if (file.name.substr(file.name.length-3,file.name.length) === ".gb") {
                    if ('name' in file) {
                        txt += "  " + file.name;
                    }
                    if ('size' in file) {
                        txt += " | " + Math.round(file.size/100)/10 + " KB";
                    }
                    txt += " | Gene name: "
                    geneName = file.name.substr(0,file.name.length-3);
                    var nameField = document.createElement("input");
                    nameField.setAttribute("type","text");
                    nameField.setAttribute("value",geneName);

                }
                else if ('name' in file) {
                    txt += "File " + file.name + " not a GenBank format file (.gb)";
                }
                p.innerHTML = txt;
                if (geneName.length > 0) {
                  p.appendChild(nameField);
                }
                p.innerHTML += "<span></span> <br>";
                parent.appendChild(p);
            }
        }
    }
    else {
        if (x.value === "") {
            txt += "Select one or more files.";
        } else {
            txt += "Sorry, the files property is not supported by your browser.";
        }
        var p = document.createElement("p");
        p.innerHTML = txt;
        parent.appendChild(p);
    }
}

// Makes back-to-top button appear if scrolled down and
function scrollEvt() {
    if (document.getElementById("Title").getBoundingClientRect().bottom < 0 && !backToTopShown) {
        fadeIn(document.getElementById("backToTopBtn"));
        backToTopShown = true;
    }
    else if (document.getElementById("Title").getBoundingClientRect().bottom > 0 && backToTopShown) {
        fadeOut(document.getElementById("backToTopBtn"));
        backToTopShown = false;
    }
}

// fade out from http://www.chrisbuttery.com/articles/fade-in-fade-out-with-javascript/
function fadeOut(el){
  if (el.style.opacity >= 1) {
    (function fade() {
      if ((el.style.opacity -= .1) < 0) {
        el.style.display = "none";
      } else {
        requestAnimationFrame(fade);
      }
    })();
  }
}

// fade in
function fadeIn(el, display){
  if (el.style.opacity <= 0) {
    el.style.display = display || "block";

    (function fade() {
      var val = parseFloat(el.style.opacity);
      if (!((val += .1) > 1)) {
        el.style.opacity = val;
        requestAnimationFrame(fade);
      }
    })();
  }
}

// the full monty
function run() {
    if (geneServerFiles.length > 0) {
        runGeneServerFiles();
    }
    else {
        document.getElementById("run").innerHTML = "Processing files...";
        var x = document.getElementById("geneFileForm");
        var txt = "";
        if ('files' in x) {
            if (x.files.length == 0) {
                txt = "Select one or more GenBank gene files.";
            }
            else {
                numFilesUploaded = x.files.length;
                var i = fileCounter;
                var file = x.files[i];
                if (file.name.substr(file.name.length-3,file.name.length) === ".gb") {
                    var fR = new FileReader();
                    fR.fileName = document.getElementById('selectedFiles').children[i].children[0].value;
                    var queryNumber = fileCounter;
                    fR.readAsText(file, "UTF-8");
                    fR.onload = function (evt) {
                        var HRann = document.getElementById("HRannChkBox").checked;
                        var lengthLHR = [document.getElementById("LHRMin").value, document.getElementById("LHRPref").value, document.getElementById("LHRMax").value];
                        var lengthRHR = [document.getElementById("RHRMin").value, document.getElementById("RHRPref").value, document.getElementById("RHRMax").value];
                        var lengthGib = [document.getElementById("gibMin").value, document.getElementById("gibPref").value, document.getElementById("gibMax").value];
                        var optimLHR = [-1*document.getElementById("optLowLHR").value, document.getElementById("optHighLHR").value];
                        var optimRHR = [-1*document.getElementById("optLowRHR").value, document.getElementById("optHighRHR").value];
                        var endsLHR = document.getElementById("endsLHR").value;
                        var endsRHR = document.getElementById("endsRHR").value;
                        var endsTempLHR = document.getElementById("endsTempLHR").value;
                        var endsTempRHR = document.getElementById("endsTempRHR").value;
                        var gibTemp = document.getElementById("gibTemp").value;
                        var gibTDif = document.getElementById("gibTDif").value;
                        var maxDistLHR = document.getElementById("maxDistLHR").value;
                        var maxDistRHR = document.getElementById("maxDistRHR").value;
                        var minFragSize = document.getElementById('minFragSize').value;
                        var optimOrg = document.getElementById('codonOptimizeOrg').value;
                        var codonSampling = document.getElementById('codonOptimStrat').value;
                        var minGCContent = document.getElementById('gRNAGCContent').value;
                        var onTargetMethod = document.getElementById('gRNAOnTargetMethod').value;
                        var onTargetScore = document.getElementById('gRNAOnTargetCutoff').value;
                        var offTargetMethod = document.getElementById('gRNAOffTargetMethod').value;
                        var offTargetScore = document.getElementById('minOffTargetScore').value;
                        var offTargetHitScore = document.getElementById('maxOffTargetHitScore').value;
                        var enzyme = document.getElementById('enzymeType').value;
                        var pam = document.getElementById('PAMSequence').value;
                        var gBlockDefault = document.getElementById('gBlockDefault').checked;
                        var plasmidType = document.getElementById('plasmidType').value;
                        var haTag = document.getElementById('haTag').value;
                        var setCoding = document.getElementById('setCoding').value;
                        bulkFile = document.getElementById("bulkFileChkBox").checked;
                        var prefix = document.getElementById("prefix").value;
                        var prefixNum = document.getElementById("prefixNum").value;
                        var locationType = document.getElementById("locationType").value;
                        var basePlasmid;
                        var basePlasmidName;

                        if (plasmidType === "custom") {
                            if (document.getElementById("plasmidFileForm").files.length == 1) {
                                var basePlasmidFile = document.getElementById("plasmidFileForm").files[0];
                                var fR2 = new FileReader();
                                fR2.fileName = document.getElementById('selectedPlasmidFile').children[0].children[0].value;
                                fR2.readAsText(basePlasmidFile, "UTF-8");
                                fR2.onload = function (evt2) {
                                    basePlasmid = evt2.target.result
                                    basePlasmidName = evt2.target.fileName

                                    msg = createFileMsg([queryNumber, evt.target.result, evt.target.fileName,
                                      HRann, lengthLHR, lengthRHR, lengthGib, optimLHR, optimRHR, endsLHR, endsRHR,
                                      endsTempLHR, endsTempRHR, gibTemp, gibTDif, maxDistLHR, maxDistRHR, minFragSize,
                                      optimOrg, codonSampling, minGCContent, onTargetMethod, onTargetScore, offTargetMethod,
                                      offTargetScore, offTargetHitScore, enzyme, pam, gBlockDefault,
                                      plasmidType, haTag, setCoding, bulkFile, prefix, prefixNum,
                                      basePlasmid, basePlasmidName, locationType]);
                                    sendMessageToServer('Sending requests...', "misc");
                                    sendMessageToServer(msg,'sendGeneFile');
                                    queryNumber += 1;
                                }

                                fR2.onerror = function (evt) {
                                    var errMsg = "Error reading file ";
                                    if ('name' in file) {
                                        errMsg += file.name;
                                    }
                                    document.getElementById("outputLog").innerHTML = errMsg;
                                }
                            }
                            else {
                                window.alert("Upload a single GenBank (.gb) file as the custom base plasmid.");
                            }
                        }
                        else {
                            msg = createFileMsg([queryNumber, evt.target.result, evt.target.fileName,
                              HRann, lengthLHR, lengthRHR, lengthGib, optimLHR, optimRHR, endsLHR, endsRHR,
                              endsTempLHR, endsTempRHR, gibTemp, gibTDif, maxDistLHR, maxDistRHR, minFragSize,
                              optimOrg, codonSampling, minGCContent, onTargetMethod, onTargetScore, offTargetMethod,
                              offTargetScore, offTargetHitScore, enzyme, pam, gBlockDefault,
                              plasmidType, haTag, setCoding, bulkFile, prefix, prefixNum,
                              basePlasmid, basePlasmidName, locationType]);
                            sendMessageToServer('Sending requests...', "misc");
                            sendMessageToServer(msg,'sendGeneFile');
                            queryNumber += 1;
                        }
                    }
                    fR.onerror = function (evt) {
                        var errMsg = "Error reading file ";
                        if ('name' in file) {
                            errMsg += file.name;
                        }
                        document.getElementById("outputLog").innerHTML = errMsg;
                    }
                }
                else if ('name' in file) {
                    txt += "File " + file.name + " not a GenBank format file (.gb)";
                }
            }
        }
        else {
            if (x.value === "") {
                txt += "Select one or more files.";
            } else {
                txt += "Sorry, the files property is not supported by your browser.";
            }
        }
    }
}

function runGeneServerFiles() {
    document.getElementById("run").innerHTML = "Processing files...";
    var txt = "";
    numFilesUploaded = geneServerFiles.length;
    var fileName = geneServerNames[fileCounter].trim()

    if (fileName.length > 0) {
        var HRann = document.getElementById("HRannChkBox").checked;
        var lengthLHR = [document.getElementById("LHRMin").value, document.getElementById("LHRPref").value, document.getElementById("LHRMax").value];
        var lengthRHR = [document.getElementById("RHRMin").value, document.getElementById("RHRPref").value, document.getElementById("RHRMax").value];
        var lengthGib = [document.getElementById("gibMin").value, document.getElementById("gibPref").value, document.getElementById("gibMax").value];
        var optimLHR = [-1*document.getElementById("optLowLHR").value, document.getElementById("optHighLHR").value];
        var optimRHR = [-1*document.getElementById("optLowRHR").value, document.getElementById("optHighRHR").value];
        var endsLHR = document.getElementById("endsLHR").value;
        var endsRHR = document.getElementById("endsRHR").value;
        var endsTempLHR = document.getElementById("endsTempLHR").value;
        var endsTempRHR = document.getElementById("endsTempRHR").value;
        var gibTemp = document.getElementById("gibTemp").value;
        var gibTDif = document.getElementById("gibTDif").value;
        var maxDistLHR = document.getElementById("maxDistLHR").value;
        var maxDistRHR = document.getElementById("maxDistRHR").value;
        var minFragSize = document.getElementById('minFragSize').value;
        var optimOrg = document.getElementById('codonOptimizeOrg').value;
        var codonSampling = document.getElementById('codonOptimStrat').value;
        var minGCContent = document.getElementById('gRNAGCContent').value;
        var onTargetMethod = document.getElementById('gRNAOnTargetMethod').value;
        var onTargetScore = document.getElementById('gRNAOnTargetCutoff').value;
        var offTargetMethod = document.getElementById('gRNAOffTargetMethod').value;
        var offTargetScore = document.getElementById('minOffTargetScore').value;
        var offTargetHitScore = document.getElementById('maxOffTargetHitScore').value;
        var enzyme = document.getElementById('enzymeType').value;
        var pam = document.getElementById('PAMSequence').value;
        var gBlockDefault = document.getElementById('gBlockDefault').checked;
        var plasmidType = document.getElementById('plasmidType').value;
        var haTag = document.getElementById('haTag').value;
        var setCoding = document.getElementById('setCoding').value;
        bulkFile = document.getElementById("bulkFileChkBox").checked;
        var prefix = document.getElementById("prefix").value;
        var prefixNum = document.getElementById("prefixNum").value;
        var locationType = document.getElementById("locationType").value;
        var basePlasmid;
        var basePlasmidName;

        if (plasmidType === "custom") {
            if (document.getElementById("plasmidFileForm").files.length == 1) {
                var basePlasmidFile = document.getElementById("plasmidFileForm").files[0];
                var fR2 = new FileReader();
                fR2.fileName = document.getElementById('selectedPlasmidFile').children[0].children[0].value;
                fR2.readAsText(basePlasmidFile, "UTF-8");
                fR2.onload = function (evt2) {
                    basePlasmid = evt2.target.result
                    basePlasmidName = evt2.target.fileName

                    msg = createFileMsg([fileCounter, geneServerFiles[fileCounter], fileName,
                      HRann, lengthLHR, lengthRHR, lengthGib, optimLHR, optimRHR, endsLHR, endsRHR,
                      endsTempLHR, endsTempRHR, gibTemp, gibTDif, maxDistLHR, maxDistRHR, minFragSize,
                      optimOrg, codonSampling, minGCContent, onTargetMethod, onTargetScore, offTargetMethod,
                      offTargetScore, offTargetHitScore, enzyme, pam, gBlockDefault, plasmidType,
                      haTag, setCoding, bulkFile, prefix, prefixNum,
                      basePlasmid, basePlasmidName, locationType]);
                    sendMessageToServer('Sending requests...', "misc");
                    sendMessageToServer(msg,'sendGeneFile');
                }

                fR2.onerror = function (evt) {
                    var errMsg = "Error reading file ";
                    if ('name' in file) {
                        errMsg += file.name;
                    }
                    document.getElementById("outputLog").innerHTML = errMsg;
                }
            }
            else {
                window.alert("Upload a single GenBank (.gb) file as the custom base plasmid.");
            }
        }
        else {
          msg = createFileMsg([fileCounter, geneServerFiles[fileCounter], fileName,
            HRann, lengthLHR, lengthRHR, lengthGib, optimLHR, optimRHR, endsLHR, endsRHR,
            endsTempLHR, endsTempRHR, gibTemp, gibTDif, maxDistLHR, maxDistRHR, minFragSize,
            optimOrg, codonSampling, minGCContent, onTargetMethod, onTargetScore, offTargetMethod,
            offTargetScore, offTargetHitScore, enzyme, pam, gBlockDefault, plasmidType,
            haTag, setCoding, bulkFile, prefix, prefixNum,
            basePlasmid, basePlasmidName, locationType]);
          sendMessageToServer('Sending requests...', "misc");
          sendMessageToServer(msg,'sendGeneFile');
        }
    }
}


function displayGeneOutput(files) {
    var clean = true;
    var msgLog = files[files.length-1];
    var fileNum;
    if (haTagged === "_HA_Tag") {
        fileNum = parseInt(files[0].split("-")[1]);
    }
    else {
        fileNum = parseInt(files[0]);
        fileCounter += 1;
    }
    var fileLine = document.getElementById('selectedFiles').children;
    if (msgLog.search(/warning/i) > -1) {
        fileLine[fileNum].style.color = "#FEB36E";
        fileLine[fileNum].children[1].innerHTML = " <b> || WARNINGS</b>"
        clean = false;
    }
    if (msgLog.search(/error/i) > -1) {
        fileLine[fileNum].style.color = "#DA717A";
        fileLine[fileNum].children[1].innerHTML = " <b> || ERROR</b>"
    }
    else if (clean) {
        fileLine[fileNum].style.color = "#73A46F";
        fileLine[fileNum].children[1].innerHTML = " <b> || SUCCESS!</b>"
    }

    if (fileCounter == numFilesUploaded) {
        document.getElementById("run").innerHTML = "Completed!";
    }
    else if (haTagged !== "_HA_Tag") {
        run();
    }
}


function downloadOutput() {
    var fileTypes = ["Locus_Pre-editing_","Plasmid_","Locus_Post-editing_","Oligos_","gBlocks_","sgRNA_Table_","Message_File_"];
    var enzyme = document.getElementById('enzymeType').value;
    var plasmid = document.getElementById('plasmidType').value;
    var fileExt = [".gb",".gb",".gb",".csv",".fasta",".csv",".txt"];

    var zip = new JSZip();

    if (!bulkFile) {
        for (var j = 0; j < currentOutput.length; j++) {
            if (j%(fileTypes.length+1) > 0 && currentOutput[j].length > 1) {
                var geneName = currentOutput[Math.floor(j/(fileTypes.length+1))*(fileTypes.length+1)] + nonCoding;
                zip.folder(geneName).file(fileTypes[j%(fileTypes.length+1)-1]+geneName+"_"+plasmid+"_"+enzyme+haTagged+fileExt[j%(fileTypes.length+1)-1], currentOutput[j]);
                if (j==currentOutput.length-1) {
                    downloadNamedZip(zip,geneName);
                }
            }
        }
    }
    else {
        for (var j = 0; j < currentOutput.length; j++) {
            if (j%(fileTypes.length+1) > 0 && currentOutput[j].length > 1) {
                if (fileExt[j%(fileTypes.length+1)-1] === ".csv" && bulkFileArray[j%(fileTypes.length+1)-1].length > 0) {
                    currentOutput[j] = currentOutput[j].substring(currentOutput[j].indexOf('\n'),currentOutput[j].length);
                }
                bulkFileArray[j%(fileTypes.length+1)-1] += currentOutput[j]+'\n\n';
            }
        }
        if (fileCounter == numFilesUploaded) {
            for (var j = 0; j < bulkFileArray.length; j++) {
                if (fileTypes[j] === "Oligos_") {
                    var lines = bulkFileArray[j].split('\n');
                    var lhrF = [];
                    var lhrR = [];
                    var rhrF = [];
                    var rhrR = [];
                    var grnaF = [];
                    var grnaR = [];
                    var gblkF = [];
                    var gblkR = [];
                    var otherF = [];
                    var otherR = [];
                    for (var k = 1; k < lines.length; k++) {
                        lines[k] = lines[k].trim()
                        if (lines[k].length > 0) {
                            var oligoTypePortion = lines[k].substring(lines[k].indexOf(',')-10,lines[k].indexOf(','));
                            switch (oligoTypePortion.substring(oligoTypePortion.indexOf('_')+1)) {
                                case "LHR_F":
                                    lhrF.push(lines[k]);
                                    break;
                                case "LHR_R":
                                    lhrR.push(lines[k]);
                                    break;
                                case "RHR_F":
                                    rhrF.push(lines[k]);
                                    break;
                                case "RHR_R":
                                    rhrR.push(lines[k]);
                                    break;
                                case "gRNA_F":
                                    grnaF.push(lines[k]);
                                    break;
                                case "gRNA_R":
                                    grnaR.push(lines[k]);
                                    break;
                                case "gBlock_F":
                                    gblkF.push(lines[k]);
                                    break;
                                case "gBlock_R":
                                    gblkR.push(lines[k]);
                                    break;
                                default:
                                    if (lines[k].includes('_F')) {
                                        otherF.push(lines[k]);
                                    }
                                    else {
                                        otherR.push(lines[k]);
                                    }
                            }
                        }
                    }
                    bulkFileArray[j] = lines[0] + '\n' + rhrF.join('\n') + '\n' + rhrR.join('\n') + '\n' + lhrF.join('\n') + '\n' + lhrR.join('\n') + '\n' + grnaF.join('\n') + '\n' + grnaR.join('\n') + '\n' + gblkF.join('\n') + '\n' + gblkR.join('\n');
                    if (otherF.length > 0) {
                        bulkFileArray[j] += otherF.join('\n') + '\n' + otherR.join('\n');
                    }
                }

                var prefix = document.getElementById("prefix").value;
                zip.folder('GeneTargeter_bulk').file(prefix+"_"+fileTypes[j]+plasmid+"_"+enzyme+haTagged+fileExt[j], bulkFileArray[j]);
                if (j==bulkFileArray.length-1) {
                    downloadNamedZip(zip,'GeneTargeter_bulk');
                }
            }

            bulkFileArray = ["","","","","","",""];
        }
    }

}

function downloadNamedZip(mainZip,fName) {
    mainZip.generateAsync({type: 'blob'}).then(content => {
        saveAs(content, fName+'_output.zip');
    });
}

function createFileMsg(info) {
    var sep = ":::";
    var msg = "";
    for (var i = 0; i < info.length; i++) {
        msg = msg + info[i] + sep;
    }
    msg = msg.substr(0,msg.length-sep.length);
    return msg;
}


function decodeFileMsg(content) {
    var sep = ":::";
    // var fileSep = "|:||||:|";

    // files = content.data.split(fileSep);
    // for (var i = 0; i < files.length; i++) {
    //     files[i] = files[i].split(sep);
    // }

    if (content.data.search("Note: this gene appears not to be a protein-coding sequence") > 0) {
        nonCoding = "_putative_ncRNA";

    }
    else {
        nonCoding = "";
    }

    if (content.data.search("HA_Tag_Design-") > -1) {
        haTagged = "_HA_Tag";
    }
    else {
        haTagged = "";
    }

    files = content.data.split(sep);
    currentOutput = files.slice(1,files.length); // first file is actually number, not a file
    return files;
}

function changeOffScoringMethod() {
    if (document.getElementById('gRNAOffTargetMethod').value === 'hsu') {
        document.getElementById('minOffTargetScore').value = 75;
        document.getElementById('maxOffTargetHitScore').value = 5;
    }
    else if (document.getElementById('gRNAOffTargetMethod').value === 'cfd') {
        document.getElementById('minOffTargetScore').value = 20;
        document.getElementById('maxOffTargetHitScore').value = 50;
    }
}

function changeOnScoringMethod() {
    if (document.getElementById('gRNAOnTargetMethod').value === 'azimuth') {
        document.getElementById('gRNAOnTargetCutoff').value = 35;
    }
    else if (document.getElementById('gRNAOnTargetMethod').value === 'cindel') {
        document.getElementById('gRNAOnTargetCutoff').value = 20;
    }
}

function changeEnzyme() {
    if (document.getElementById('enzymeType').value === 'Cas9') {
        document.getElementById('gRNAOnTargetMethod').value = "azimuth";
        document.getElementById('gRNAOffTargetMethod').value = "cfd";
        document.getElementById('PAMSequence').value = "NGG";
    }
    else if (document.getElementById('enzymeType').value === 'Cas12') {
        document.getElementById('gRNAOnTargetMethod').value = "cindel";
        document.getElementById('gRNAOffTargetMethod').value = "hsu";
        document.getElementById('PAMSequence').value = "TTTV";
    }
    changeOffScoringMethod();
    changeOnScoringMethod();
}

function changeUTRTarget() {
    if (document.getElementById('plasmidType').value === 'pSN054' || document.getElementById('plasmidType').value === 'pSN054_V5') {
        document.getElementById('locationType').value = '3prime';
    }
    else if (document.getElementById('plasmidType').value === 'pSN150') {
        document.getElementById('locationType').value = '5prime';
    }
    else if (document.getElementById('plasmidType').value === 'pSN150-KO') {
        document.getElementById('locationType').value = 'center';
    }

    if (document.getElementById('locationType').value === '5prime') {
        document.getElementById('maxDistLHRTxt').textContent = document.getElementById('maxDistLHRTxt').textContent.replace('gRNA','gene');
        document.getElementById('maxDistRHRTxt').textContent = document.getElementById('maxDistRHRTxt').textContent.replace('gene','gRNA');
        // document.getElementById('maxDistRHR').value = 500;
    }
    else {
      document.getElementById('maxDistRHRTxt').textContent = document.getElementById('maxDistRHRTxt').textContent.replace('gRNA','gene');
      document.getElementById('maxDistLHRTxt').textContent = document.getElementById('maxDistLHRTxt').textContent.replace('gene','gRNA');
      // if (document.getElementById('plasmidType').value === 'pSN150-KO') {
      //     document.getElementById('maxDistRHR').value = 5000;
      // }
    }
    changeOffScoringMethod();
    changeOnScoringMethod();
}
