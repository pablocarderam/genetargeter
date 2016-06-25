/* JS Document */
/* Handles client-server connection */
// Create SocketIO instance, connect

namespace = '/link'; // change to an empty string to use the global namespace

// the socket.io documentation recommends sending an explicit package upon connection
// this is specially important when using the global namespace

//var socket = io.connect('https://' + document.domain + ':' + location.port + namespace);
// Changing from the above line in the original example to the following allows the
// system to work locally with http, and on Heroku with both https and http:
var socket = io.connect(namespace);

// Add a connect listener
socket.on('connect',function() {
  console.log('Client has connected to the server!');
});
// Add a connect listener
socket.on('message',function(data) {
  console.log('Received a message from the server!',data);
  document.alert('Received a message from the server!',data)
});
// Add a disconnect listener
socket.on('disconnect',function() {
  console.log('The client has disconnected!');
});

// Sends a message to the server via sockets
function sendMessageToServer(message) {
  socket.emit('sendGeneFile', {data: message});
  console.log(message['data']+" :: sent");
};
