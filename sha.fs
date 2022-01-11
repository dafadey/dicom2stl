#version 150
 
in vec3 color;
out vec4 outputF;

void main()
{
    outputF = vec4(color, 1.0);
}
