#version 150
 
in vec2 circle_pos;
in vec3 circle_color;
in float circle_resolution;
in vec2 circle_radii;

out vec3 vColor;
out float vSides; 
out vec2 radii; 

void main()
{
    vColor = circle_color;
    vSides = circle_resolution;
    radii = circle_radii;
    vec4 pos = vec4(circle_pos.x, circle_pos.y, 0, 1);
    gl_Position = pos;
}
