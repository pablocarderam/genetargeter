/* JS Document */
/* Main js */

var backToTopShown = true; // used to control back-to-top button appearance and disappearance
var currentOutput = []; // saves current output array to download multiple times if desired

function init() {
    EPPZScrollTo.scrollVerticalToElementById('Title', 0); // scroll to start at init
}

function askCredentials() {
    document.getElementById('myModal').style.display = "block";
    document.getElementById("modFooter").innerHTML = "Type passcode";
    document.getElementById("passcode").value = "";
    if (document.getElementById("modContent").classList.contains("modal-out")) {
        document.getElementById("modContent").classList.remove("modal-out");
    }
    if (!document.getElementById("modContent").classList.contains("modal-in")) {
        document.getElementById("modContent").classList.add("modal-in");
    }
}

function verifyCredentials() {
    sendMessageToServer(document.getElementById("passcode").value,"cred");
}

function invalidCred() {
    document.getElementById("modFooter").innerHTML = "Invalid passcode.";
}

function validCred() {
    var x = document.getElementById("geneFileForm");
    if ('files' in x) {
        if (x.files.length == 0) {
            document.getElementById("modFooter").innerHTML = "Valid passcode! Select some gene files though.";
        }
        else {
            document.getElementById("modFooter").innerHTML = "Valid passcode!";
            closeModal();
            run();
        }
    }
    else {
        window.alert("Application error. Something went wrong because the gene file form is supposed to have a files property, but I can't find it. Sorry :P")
    }
}

function checkKey() {
    console.log("keydown");
    if (event.keyCode == 13) {
        verifyCredentials();
    }
}

// When the user clicks on <span> (x), close the modal
function closeModal() {
    if (document.getElementById("modContent").classList.contains("modal-in")) {
        document.getElementById("modContent").classList.remove("modal-in");
    }
    if (!document.getElementById("modContent").classList.contains("modal-out")) {
        document.getElementById("modContent").classList.add("modal-out");
    }
    setTimeout(displayNone, 250);
}

function displayNone() {
    document.getElementById('myModal').style.display = "none";
}

// When the user clicks anywhere outside of the modal, close it
window.onclick = function(event) {
    if (event.target == document.getElementById('myModal')) {
        closeModal();
    }
}

// Load file adapted from http://www.w3schools.com/jsref/prop_fileupload_files.asp
function uploadFile(){
    var x = document.getElementById("geneFileForm");
    var txt = "";
    if ('files' in x) {
        if (x.files.length == 0) {
            txt = "Select one or more GenBank gene files.";
        }
        else {
            for (var i = 0; i < x.files.length; i++) {
                txt += "<br>" + (i+1) + ".";
                var file = x.files[i];
                if (file.name.substr(file.name.length-3,file.name.length) === ".gb") {
                    if ('name' in file) {
                        txt += "  " + file.name;
                    }
                    if ('size' in file) {
                        txt += " || " + Math.round(file.size/100)/10 + " KB";
                    }
                }
                else if ('name' in file) {
                    txt += "File " + file.name + " not a GenBank format file (.gb)";
                }
            }
        }
    }
    else {
        if (x.value == "") {
            txt += "Select one or more files.";
        } else {
            txt += "Sorry, the files property is not supported by your browser.";
        }
    }

    document.getElementById("selectedFiles").innerHTML = txt;
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
  el.style.opacity = 1;

  (function fade() {
    if ((el.style.opacity -= .1) < 0) {
      el.style.display = "none";
    } else {
      requestAnimationFrame(fade);
    }
  })();
}

// fade in
function fadeIn(el, display){
  el.style.opacity = 0;
  el.style.display = display || "block";

  (function fade() {
    var val = parseFloat(el.style.opacity);
    if (!((val += .1) > 1)) {
      el.style.opacity = val;
      requestAnimationFrame(fade);
    }
  })();
}

function run() {
    var x = document.getElementById("geneFileForm");
    var txt = "";
    if ('files' in x) {
        if (x.files.length == 0) {
            txt = "Select one or more GenBank gene files.";
        }
        else {
            sendMessageToServer('Sending requests...', "misc")
            for (var i = 0; i < x.files.length; i++) {
                var file = x.files[i];
                if (file.name.substr(file.name.length-3,file.name.length) === ".gb") {
                    var fR = new FileReader();
                    fR.fileName = file.name;
                    fR.readAsText(file, "UTF-8");
                    fR.onload = function (evt) {
                        if ('name' in file) {
                            var geneName = evt.target.fileName.substr(0,evt.target.fileName.length-3);
                            var HRann = document.getElementById("HRannChkBox").checked;
                            lengthLHR = [document.getElementById("LHRMin").value, document.getElementById("LHRPref").value, document.getElementById("LHRMax").value];
                            lengthRHR = [document.getElementById("RHRMin").value, document.getElementById("RHRPref").value, document.getElementById("RHRMax").value];
                            lengthGib = [document.getElementById("gibMin").value, document.getElementById("gibPref").value, document.getElementById("gibMax").value];
                            msg = createFileMsg(evt.target.result, geneName, HRann, lengthLHR, lengthRHR, lengthGib);
                            sendMessageToServer(msg,'sendGeneFile');
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
    }
    else {
        if (x.value == "") {
            txt += "Select one or more files.";
        } else {
            txt += "Sorry, the files property is not supported by your browser.";
        }
    }
}


function displayGeneOutput() {

}


function downloadOutput() {
    for (var i = 0; i < currentOutput.length; i++) {
        var geneName = currentOutput[i][0];
        var fileTypes = ["_Annotated_Gene.gb","_Annotated_Plasmid.gb","_Oligos.txt","_Message_File.txt"];
        for (var j = 1; j < currentOutput[i].length; j++) {
            var file = currentOutput[i][j];
            var fileType = fileTypes[j-1];
            var data = 'data:text/plain;charset=utf-8,' + encodeURIComponent(file);
            saveAs(data,geneName+fileType);
        }
    }
}


function createFileMsg(content, name, HRann, lengthLHR, lengthRHR, lengthGib) {
    var sep = ":::";
    var HRannStr = "FALSE";
    if (HRann) {
        HRannStr = "TRUE";
    }
    return name + sep + content + sep + HRannStr + sep + lengthLHR.toString() + sep + lengthRHR.toString() + sep + lengthGib.toString();
}


function decodeFileMsg(content) {
    var sep = ":::";
    var fileSep = "|:||||:|";

    files = content.data.split(fileSep);
    for (var i = 0; i < files.length; i++) {
        files[i] = files[i].split(sep);
    }

    currentOutput = files;
    return files;
}


function saveAs(uri, filename) {
    var link = document.createElement('a');
    if (typeof link.download === 'string') {
        document.body.appendChild(link); // Firefox requires the link to be in the body
        link.download = filename;
        link.href = uri;
        link.click();
        document.body.removeChild(link); // remove the link when done
    } else {
        location.replace(uri);
    }
}
