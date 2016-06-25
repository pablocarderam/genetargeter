/* JS Document */
/* Main js */

var backToTopShown = false; // used to control back-to-top button appearance and disappearance

function init() {
    EPPZScrollTo.scrollVerticalToElementById('Title', 0); // scroll to start at init
}

// Load file adapted from http://www.w3schools.com/jsref/prop_fileupload_files.asp
function uploadFile(){
    var x = document.getElementById("geneFileForm");
    var txt = "";
    if ('files' in x) {
        if (x.files.length == 0) {
            txt = "Select one or more files.";
        }
        else {
            for (var i = 0; i < x.files.length; i++) {
                txt += "<br>" + (i+1) + ".";
                var file = x.files[i];
                if ('name' in file) {
                    txt += "  " + file.name;
                }
                if ('size' in file) {
                    txt += " || " + Math.round(file.size/100)/10 + " KB";
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
    sendMessageToServer('Testing request to server')
}
