#version 330

// Output fragment color
out vec4 FragColor;

// Input attributes from the vertex shader
in vec3 pos;           // Transformed vertex position (model and view applied) used for lighting calculations
in vec3 pos_fixed;     // Vertex position with only view transform (remains fixed regardless of model rotation)
in vec3 normal;        // Transformed normal (rotated, so it follows object rotation)
in vec3 normal_fixed;  // Original normal (remains constant, used for fixed effects)

struct Light {
    vec3 position;   // Light source position
    vec3 ambient;    // Ambient light color
    vec3 diffuse;    // Diffuse light color
    vec3 specular;   // Specular light color
    vec3 position2;  // Secondary light source position (unused)
};
uniform Light light;            // Light parameters
uniform float shininess;        // Shininess exponent for specular highlights
uniform vec3 viewpos;           // Camera (viewer) position

// Rendering mode: 1: points, 2: lines, 3: smooth surface, 4: flat surface
uniform int rendermode = 3;
uniform int keep_positive = 0;

// Color mode controls the final color output:
// 1: Lighting mode (object color affected by rotation)
// 2: Fixed normal mode (object color unaffected by rotation)
// 3: Custom color mode (no lighting effect)
// 4: Use VBO provided color (no lighting effect)
// 5: Use VBO provided color (with lighting effect)
uniform int colormode = 1;

// Custom color used when colormode==3
uniform vec4 color = vec4(1, 0, 0, 1.0);

void main()
{
    vec3 n; // Final normal used for lighting calculations

    // Determine the normal vector 'n' based on the selected color mode
    if (colormode == 1) {
        // Lighting mode: the object normal changes with rotation
        if (rendermode == 4)
            n = normalize(cross(dFdx(pos), dFdy(pos)));  // Discrete mode: compute normal from screen-space derivatives
        else
            n = normalize(normal);  // Continuous mode: use the passed-in transformed normal
    }
    else if (colormode == 2) {
        // Fixed normal mode: the normal remains fixed regardless of rotation
        if (rendermode == 4)
            n = normalize(cross(dFdx(pos_fixed), dFdy(pos_fixed)));
        else
            n = normalize(normal_fixed);
    }
    else if (colormode == 3) {
        // Custom color mode: ignore lighting, output the custom color directly
        FragColor = color;
        return;
    }
    else if (colormode == 606) {
        // HexMesh: face color (a specific gray)
        FragColor = vec4(0.7725, 0.7725, 0.7725, 1.0);
        return;
    }
    else if (colormode == 530) {
        // HexMesh: edge color (black)
        FragColor = vec4(0, 0, 0, 1.0);
        return;
    }
    else if (colormode == 231) {
        // Ground color (purple-ish)
        FragColor = vec4(0.64, 0.5, 0.98, 1.0);
        return;
    }
    else if (colormode == 232) {
        // Red color
        FragColor = vec4(1.0, 0, 0, 1.0);
        return;
    }
    else if (colormode == 233) {
        // Green color
        FragColor = vec4(0.0, 1, 0, 1.0);
        return;
    }
    else if (colormode == 234) {
        // Blue color
        FragColor = vec4(0.1, 0.1, 0.8, 1.0);
        return;
    }
    else if (colormode == 235) {
        // Skeleton color: semi-transparent blue
        FragColor = vec4(0, 0, 1, 0.5);
        return;
    }
    else if (colormode == 236) {
        // Skeleton triangle color: semi-transparent dark gray
        FragColor = vec4(0.2, 0.2, 0.2, 0.2);
        return;
    }
    else if (colormode == 237) {// SpringGreen
        FragColor = vec4(0.0, 1.0, 0.5, 0.7);
        return;
    }
    else if (colormode == 238) {// Salmon
        FragColor = vec4(0.98, 0.5, 0.447, 0.8);
        return;
    }
    else if (colormode == 510) {
        // Orange color
        FragColor = vec4(1.0, 0.4, 0.2, 1.0);
        return;
    }
    else if (colormode == 239) {
        // Black color
        FragColor = vec4(0, 0, 0, 1.0);
        return;
    }
    else if (colormode == 999) {
        // Red color variant with transparency
        FragColor = vec4(1, 0, 0, 0.7);
        return;
    }
    else if (colormode == 888) {
        // Fixed point color
        FragColor = vec4(0.35, 0.4, 0.9333, 0.7);
        return;
    }
    else if (colormode == 518) {
        // Translucence: output the fixed normal with 50% opacity
        FragColor = vec4(normal_fixed, 0.5);
        return;
    }
    else if (colormode < 0) {
        // User-defined color mode: output the fixed normal as color
        FragColor = vec4(normal_fixed, 1.0);
        return;
    }
    else if (colormode == 222 || colormode == 223) {
        // For these modes, use the continuous transformed normal
        n = normalize(normal);
    }

    // --- Lighting Calculations ---

    // Ambient lighting component
    vec3 ambient = light.ambient;
    
    // Diffuse lighting component
    vec3 lightDir = normalize(light.position - pos);
    //float diff = max(dot(n, lightDir), 0.0);
    float diff = dot(n, lightDir);
    if (keep_positive == 1 && diff < 0) {
        n *= -1;
        diff = dot(n, lightDir);
    }
    vec3 diffuse = light.diffuse * diff;
    
    // Specular lighting component using the Phong model
    vec3 viewDir = normalize(viewpos - pos);
    vec3 reflectDir = reflect(-lightDir, n);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), shininess);
    vec3 specular = light.specular * spec;
    
    // Combine all lighting components
    vec3 result = ambient + diffuse + specular;
    
    // Output final color based on colormode and rendermode settings
    if (colormode == 1) {
        // Lighting mode: choose output based on rendermode
        switch(rendermode) {
            case 1:
            case 2:
                // For points or lines, emphasize ambient light
                FragColor = vec4(ambient * 1.8, 1.0);
                break;
            case 3:
            case 4:
                // For smooth or flat surfaces, output full lighting result
                FragColor = vec4(result, 1.0);
                break;
        }
    }
    else if (colormode == 222) {
        FragColor = vec4(result, 1.0);
    }
    else if (colormode == 223) {
        FragColor = vec4(result, 1.0);
    }
    else if (colormode == 2) {
        // Fixed normal mode: map the normal vector to color space [0,1]
        FragColor = vec4((n + vec3(1.0)) / 2.0, 1.0);
    }
    else if (colormode == 3) {
        // Custom color mode: output a hardcoded red
        FragColor = vec4(0.8, 0.0, 0.0, 1.0);
    }
    else {
        // Default color if none of the modes match: dark gray
        FragColor = vec4(0.1, 0.1, 0.1, 1.0);
    }
}
