/* JS Document */
/* Main js */

var backToTopShown = false; // used to control back-to-top button appearance and disappearance
var currentOutput = []; // saves current output array to download multiple times if desired

function init() {
    EPPZScrollTo.scrollVerticalToElementById('Title', 0); // scroll to start at init
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
                    fR.readAsText(file, "UTF-8");
                    fR.onload = function (evt) {
                        if ('name' in file) {
                            var geneName = file.name.substr(0,file.name.length-3);
                            msg = createFileMsg(evt.target.result,geneName,document.getElementById("HRannChkBox").checked);
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


function createFileMsg(content, name, HRann) {
    var sep = ":::";
    var HRannStr = "FALSE";
    if (HRann) {
        HRannStr = "TRUE";
    }
    return name + sep + HRannStr + sep + content;
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
