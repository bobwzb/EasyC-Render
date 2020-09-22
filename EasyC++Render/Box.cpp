#include "pch.h"
#include "C++Render.h"
#include <iostream>

struct VertexInput { Vec3f pos; Vec2f uv; Vec3f color; };

void draw_plane(Render& rh, int a, int b, int c, int d, VertexInput mesh[], VertexInput vs_input[]) {
	mesh[a].uv.x = 0, mesh[a].uv.y = 0, mesh[b].uv.x = 0, mesh[b].uv.y = 1;
	mesh[c].uv.x = 1, mesh[c].uv.y = 1, mesh[d].uv.x = 1, mesh[d].uv.y = 0;

	vs_input[0] = mesh[a];
	vs_input[1] = mesh[b];
	vs_input[2] = mesh[c];
	rh.DrawPrimitive();

	vs_input[0] = mesh[c];
	vs_input[1] = mesh[d];
	vs_input[2] = mesh[a];
	rh.DrawPrimitive();
}

int main(void)
{
	// Initialize the render
	Render rh(800, 600);

	const int VARYING_TEX = 0;
	const int VARYING_COLOR = 1;

	VertexInput vs_input[3];

	//model data of the box
	VertexInput mesh[] = {
		{ {  1, -1,  1, }, { 0, 0 }, { 1.0f, 0.2f, 0.2f }, },
		{ { -1, -1,  1, }, { 0, 1 }, { 0.2f, 1.0f, 0.2f }, },
		{ { -1,  1,  1, }, { 1, 1 }, { 0.2f, 0.2f, 1.0f }, },
		{ {  1,  1,  1, }, { 1, 0 }, { 1.0f, 0.2f, 1.0f }, },
		{ {  1, -1, -1, }, { 0, 0 }, { 1.0f, 1.0f, 0.2f }, },
		{ { -1, -1, -1, }, { 0, 1 }, { 0.2f, 1.0f, 1.0f }, },
		{ { -1,  1, -1, }, { 1, 1 }, { 1.0f, 0.3f, 0.3f }, },
		{ {  1,  1, -1, }, { 1, 0 }, { 0.2f, 1.0f, 0.3f }, },
	};

	//texture
	Bitmap texture(256, 256);
	for (int y = 0; y < 256; y++) {
		for (int x = 0; x < 256; x++) {
			int k = (x / 32 + y / 32) & 1;
			texture.SetPixel(x, y, k ? 0xffffffff : 0xff3fbcef);
		}
	}

	Mat4x4f mat_model = matrix_set_rotate(-1, -0.5, 1, 1);
	Mat4x4f mat_proj = matrix_set_perspective(3.1415926f / 2, 8 / 6.0, 2.0, 500.0);
	Mat4x4f mat_view = matrix_set_lookat({ -0.7,0,1.5 }, { 0,0,0 }, { 0,0,1 });
	Mat4x4f mat_mvp = mat_model * mat_view*mat_proj;

	rh.SetVertexShader([&](int index, ShaderContext& output) -> Vec4f {
		Vec4f pos = vs_input[index].pos.xyz1() * mat_mvp;  // trasfer to 4D
		output.varying_vec2f[VARYING_TEX] = vs_input[index].uv;
		output.varying_vec4f[VARYING_COLOR] = vs_input[index].color.xyz1();
		return pos;
	});

	rh.SetPixelShader([&](ShaderContext& input) -> Vec4f {
		Vec2f coord = input.varying_vec2f[VARYING_TEX];	
		Vec4f tc = texture.Sample2D(coord);	
		return tc;		
	});

	// render && save
	draw_plane(rh, 0, 1, 2, 3,mesh,vs_input);
	draw_plane(rh, 7, 6, 5, 4, mesh, vs_input);
	draw_plane(rh, 0, 4, 5, 1, mesh, vs_input);
	draw_plane(rh, 1, 5, 6, 2, mesh, vs_input);
	draw_plane(rh, 2, 6, 7, 3, mesh, vs_input);
	draw_plane(rh, 3, 7, 4, 0, mesh, vs_input);
	rh.SaveFile("box.bmp");

	// if use windows, use mspaint to draw the result
#if defined(_WIN32) || defined(WIN32)
	system("mspaint.exe box.bmp");
#endif

	return 0;
}