#version 150
 
in vec3 fColor;
out vec4 outputF;

void main()
{
    outputF = vec4(fColor, 1.0);
}
