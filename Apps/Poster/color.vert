#version 330 

layout(location=0) in vec3 p;
layout(location=1) in vec3 n;
layout(location=2) in vec3 c;////

// TODO: 自定义颜色
// 1. 增加 Bool 变量支持用户输入颜色

uniform mat4 model=mat4(1.0);
uniform mat4 view=mat4(1.0);
uniform mat4 projection=mat4(1.0);
uniform mat4 rotate=mat4(1.0);
out vec3 pos;
out vec3 pos_fixed;
out vec3 normal;
out vec3 normal_fixed;
out vec3 color;////
void main(){

   pos=vec3(view*model*vec4(p,1.0));	
   pos_fixed = vec3(view*vec4(p,1.0));
   normal = mat3(transpose(inverse(model))) * n;////
   //normal = vec3(rotate*vec4(n,1.0));
   normal_fixed = n;
   color = c;////
   gl_Position= projection * view * model * vec4(p,1.0);

}
