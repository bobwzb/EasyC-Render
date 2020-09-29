#include "pch.h"
#include "C++Render.h"
#include <memory>
#include <thread>
using namespace std;
struct VertexInput
{
	Vec3f pos; 
	Vec2f uv; 
	Vec3f color; 
	Vec3f normal;
};

VertexInput vs_input[3];

VertexInput mesh[] = {
	{ {  1, -1,  1, }, { 0, 0 }, { 0.0f, 0.0f, 0.0f }, },
	{ { -1, -1,  1, }, { 0, 1 }, { 0.0f, 0.0f, 0.0f }, },
	{ { -1,  1,  1, }, { 1, 1 }, { 0.0f, 0.0f, 0.0f }, },
	{ {  1,  1,  1, }, { 1, 0 }, { 0.0f, 0.0f, 0.0f }, },
	{ {  1, -1, -1, }, { 0, 0 }, { 0.0f, 0.0f, 0.0f }, },
	{ { -1, -1, -1, }, { 0, 1 }, { 0.0f, 0.0f, 0.0f }, },
	{ { -1,  1, -1, }, { 1, 1 }, { 0.0f, 0.0f, 0.0f }, },
	{ {  1,  1, -1, }, { 1, 0 }, { 0.0f, 0.0f, 0.0f }, },
};

void draw_plane(Render& rh, int a, int b, int c, int d)
{
	mesh[a].uv.x = 0, mesh[a].uv.y = 0, mesh[b].uv.x = 0, mesh[b].uv.y = 1;
	mesh[c].uv.x = 1, mesh[c].uv.y = 1, mesh[d].uv.x = 1, mesh[d].uv.y = 0;

	Vec3f ab = mesh[b].pos - mesh[a].pos;
	Vec3f ac = mesh[c].pos - mesh[a].pos;
	Vec3f normal = vector_normalize(vector_cross(ac, ab));

	mesh[a].normal = normal;
	mesh[b].normal = normal;
	mesh[c].normal = normal;
	mesh[d].normal = normal;

	vs_input[0] = mesh[a];
	vs_input[1] = mesh[b];
	vs_input[2] = mesh[c];
	rh.DrawPrimitive();

	vs_input[0] = mesh[c];
	vs_input[1] = mesh[d];
	vs_input[2] = mesh[a];
	rh.DrawPrimitive();
}
void setTexture(Bitmap& texture) {
	for (int y = 0; y < 256; y++) {
		for (int x = 0; x < 256; x++) {
			int k = (x / 32 + y / 32) & 1;
			texture.SetPixel(x, y, k ? 0xffffffff : 0xbf4fbc2f);
		}
	}
}
int main(void)
{
	Render rh(800, 600);

	//set the varying and the key
	const int VARYING_TEXUV = 0;
	const int VARYING_COLOR = 1;
	const int VARYING_LIGHT = 2;

	Bitmap texture(256, 256);
	thread t1(setTexture, ref(texture));

	Mat4x4f mat_model = matrix_set_rotate(-1, -0.5, 1, 1);	
	Mat4x4f mat_view = matrix_set_lookat({ 3.5, 0, 0 }, { 0,0,0 }, { 0,0,1 });
	Mat4x4f mat_proj = matrix_set_perspective(3.1415926f * 0.5f, 8 / 6.0, 1.0, 400.0f);
	Mat4x4f mat_mvp = mat_model * mat_view * mat_proj;	

	// invert the model matrix, use to transfer the normal to word space
	Mat4x4f mat_model_it = matrix_invert(mat_model).Transpose();

	// define a light direction
	Vec3f light_dir = { 1, 0, 1 };

	rh.SetVertexShader([&](int index, ShaderContext& output) -> Vec4f {
		Vec4f pos = vs_input[index].pos.xyz1() * mat_mvp;
		output.varying_vec2f[VARYING_TEXUV] = vs_input[index].uv;
		output.varying_vec4f[VARYING_COLOR] = vs_input[index].color.xyz1();
		Vec3f normal = vs_input[index].normal;
		normal = (normal.xyz1() * mat_model_it).xyz();
		float intense = vector_dot(normal, vector_normalize(light_dir));
		// add 0.1 to avoid darkness
		intense = Max(0.0f, intense) + 0.1;
		output.varying_float[VARYING_LIGHT] = Min(1.0f, intense);
		return pos;
	});

	t1.join();
	rh.SetPixelShader([&](ShaderContext& input) -> Vec4f {
		Vec2f coord = input.varying_vec2f[VARYING_TEXUV];	
		Vec4f tc = texture.Sample2D(coord);		
		float light = input.varying_float[VARYING_LIGHT];
		return tc * light;
	});

	draw_plane(rh, 0, 1, 2, 3);
	draw_plane(rh, 7, 6, 5, 4);
	draw_plane(rh, 0, 4, 5, 1);
	draw_plane(rh, 1, 5, 6, 2);
	draw_plane(rh, 2, 6, 7, 3);
	draw_plane(rh, 3, 7, 4, 0);

	rh.SaveFile("Gouraud.bmp");

#if defined(_WIN32) || defined(WIN32)
	system("mspaint.exe Gouraud.bmp");
#endif

	return 0;
}