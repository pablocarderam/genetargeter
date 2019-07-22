/* JS Document */
/* Handles client-server connection */
// Create SocketIO instance, connect


function makeID(length) { // Makes a random ID for the user (room) https://stackoverflow.com/questions/1349404/generate-random-string-characters-in-javascript
   var result           = '';
   var characters       = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
   var charactersLength = characters.length;
   for ( var i = 0; i < length; i++ ) {
      result += characters.charAt(Math.floor(Math.random() * charactersLength));
   }
   return result;
}

namespace = '/link'; // change to an empty string to use the global namespace
id = makeID(10) // unique ID
var sep = ":::"; // separator

// the socket.io documentation recommends sending an explicit package upon connection
// this is specially important when using the global namespace

//var socket = io.connect('https://' + document.domain + ':' + location.port + namespace);
// Changing from the above line in the original example to the following allows the
// system to work locally with http, and on Heroku with both https and http:
var socket = io.connect(namespace);

// Add a connect listener
socket.on('connect',function() {
  socket.emit('join', { channel: id });
  console.log('Client has connected to the server! '+id);
  if (disconnected && validCredentials && fileCounter != numFilesUploaded) {
    run();
    console.log("Resending");
  }
});
// Add a connect listener
socket.on('message',function(data) {
  console.log('Received a message from the server!',data);
  alert('Received a message from the server!',data);
});
// Add a connect listener
socket.on('geneOutput',function(data) {
  console.log('Received gene output from the server!');
  files = decodeFileMsg(data);
  displayGeneOutput(files);
  downloadOutput();
});
// Add a disconnect listener
socket.on('disconnect',function() {
  console.log('The client has disconnected!');
  disconnected = true;
});
// Add a validated listener
socket.on('validCred',function() {
  console.log('Credentials valid');
  validCred();
});
// Add a validated listener
socket.on('invalidCred',function() {
  console.log('Credentials not valid');
  invalidCred();
});

// Sends a message to the server via sockets
function sendMessageToServer(message,type) {
  socket.emit(type, {data: message, channel: id});
  //console.log(message + " :: sent");
};
