#version 330 
out vec4 FragColor;
in vec3 pos;
in vec3 pos_fixed;
in vec3 normal;
in vec3 normal_fixed;
struct Light {
    vec3 position;
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};
uniform Light light;
uniform float shininess;
uniform vec3 viewpos;
uniform int rendermode = 3;     // 渲染模式，1：点，2：线，3：光滑面，4：平面
uniform int colormode = 1;      // 颜色模式，1：光照方式，2：根据法向方式
uniform vec4 color = vec4(1, 0, 0, 1.0);  // 自定颜色
void main() {
    vec3 n;
    if(colormode == 1){     // 光照模式下，物体法向会随着物体旋转而改变
        if(rendermode == 4)
            n = normalize(cross(dFdx(pos),dFdy(pos)));  // 离散
        else
            n = normalize(normal);  // 连续
    } else if(colormode == 2){     // 法向模式下，物体法向不会改变
        if(rendermode == 4){
            n = normalize(cross(dFdx(pos_fixed),dFdy(pos_fixed)));
        }
        else{
            n = normalize(normal_fixed);
        }
    }else if(colormode == 3){
        FragColor = color;
        return;
    }else if (colormode == 606) {// hexMesh: face color
        FragColor = vec4(0.7725, 0.7725, 0.7725, 1.0);
        return;
    } else if (colormode == 530) {// hexMesh: edge color
        //FragColor = vec4(0.37, 0.18, 0.57, 1.0);
        FragColor = vec4(0, 0, 0, 1.0);
        return;
    } else if (colormode == 231) {// ground color
        FragColor = vec4(0.64, 0.5, 0.98, 1.0);
        return;
    } else if (colormode == 232) {// red
        FragColor = vec4(1.0, 0, 0, 1.0);
        return;
    }
    else if (colormode == 233) {// green
        FragColor = vec4(0.0, 1, 0, 1.0);
        return;
    }
    else if (colormode == 234) {// blue
        FragColor = vec4(0.1, 0.1, 0.8, 1.0);
        return;
    }
    else if (colormode == 235) {// skel
        FragColor = vec4(0, 0, 1, 0.5);
        return;
    }
    else if (colormode == 236) {// skeltri
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
    else if(colormode == 239){ //Thistle1
        FragColor = vec4(1.0, 0.89, 1.0, 0.8);
        return;
    }
    else if (colormode == 510) {// orange
        FragColor = vec4(1.0, 0.4, 0.2, 1.0);
        return;
    }
    else if (colormode == 999) {// red color
        FragColor = vec4(1, 0, 0, 0.7);
        return;
    } else if (colormode == 888) {// fixed point
        FragColor = vec4(0.35, 0.4, 0.9333, 0.7);
        return;
    } else if (colormode == 518) {// transluncence
        FragColor = vec4(normal_fixed, 0.25);
        return;
    } else if (colormode < 0) {// your color
        FragColor = vec4(normal_fixed, 1.0);
        return;
    }

    // ambient
    vec3 ambient = light.ambient;

    // diffuse
    vec3 lightDir = normalize(light.position - pos);
    float diff = max(dot(n, lightDir), 0.0);
    vec3 diffuse = light.diffuse *diff;

    // specular
    vec3 viewDir = normalize(viewpos - pos);
    vec3 reflectDir = reflect(-lightDir, n);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), shininess);
    vec3 specular = light.specular * spec;

    vec3 result = ambient + diffuse + specular;

      if(colormode == 1){
        switch(rendermode){
        case 1:
        case 2:
            FragColor = vec4(ambient*1.8,1.0);
            break;
        case 3:
        case 4:
            FragColor = vec4(result,1.0);
            break;
        }
      }
      else if(colormode == 2){
        FragColor = vec4((n + vec3(1.0,1.0,1.0))/2,1.0);
      }
      else if(colormode == 3){
        FragColor = vec4(0.8,0.0,0.0,1.0);
      }      
      else{
        FragColor = vec4(0.1,0.1,0.1,1.0);
      }
}