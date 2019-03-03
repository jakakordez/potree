precision mediump float;
precision mediump int;

//attribute vec3 position;
varying vec3 pos;

void main(){
    pos = position;
    vec4 mvPosition = modelViewMatrix * vec4( position, 1.0 );
    gl_Position = projectionMatrix * mvPosition;
    gl_PointSize = 10.0;
}