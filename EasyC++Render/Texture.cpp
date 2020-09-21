#include "pch.h"
#include "C++Render.h"
int main(void)
{
	// Initialize the render
	Render rh(800, 600);

	//create a texture
	Bitmap texture(256, 256);
	for (int x = 0; x < 256; x++) {
		for (int y = 0; y < 256; y++) {
			int k = (x / 32 + y / 32) & 1;
			texture.SetPixel(y, x, k ? 0xffffffff : 0xff3fbcef);
		}
	}

	Mat4x4f mat_model = matrix_set_identity();
	Mat4x4f mat_proj = matrix_set_perspective(3.1415926f / 2, 8 / 6.0, 2.0, 500.0);
	Mat4x4f mat_view = matrix_set_lookat({ -0.7,0,1.5 }, { 0,0,0 }, { 0,0,1 });
	Mat4x4f mat_mvp = mat_model * mat_view*mat_proj;

	struct VertexInput { Vec4f pos; Vec2f texuv; };
	VertexInput vs_input[3];

	VertexInput vertex[] = {
		{ { 1, -1, -1, 1}, {0, 0} },
		{ { 1,  1, -1, 1}, {1, 0} },
		{ {-1,  1, -1, 1}, {1, 1} },
		{ {-1, -1, -1, 1}, {0, 1} },
	};

	vs_input[0] = vertex[0];
	vs_input[1] = vertex[1];
	vs_input[2] = vertex[2];
	rh.DrawPrimitive();

	vs_input[0] = vertex[2];
	vs_input[1] = vertex[3];
	vs_input[2] = vertex[0];

	const int VARYING_TEX = 0;

	rh.SetVertexShader([&](int index, ShaderContext& output) {
		Vec4f pos = vs_input[index].pos*mat_mvp;
		output.varying_vec2f[VARYING_TEX] = vs_input[index].texuv;
		return pos;
	});

	rh.SetPixelShader([&](ShaderContext& input) -> Vec4f {
		Vec2f coord = input.varying_vec2f[VARYING_TEX];
		return texture.Sample2D(coord);		
	});

	// render && save
	rh.DrawPrimitive();
	rh.SaveFile("output2.bmp");

	// if use windows, use mspaint to draw the result
#if defined(_WIN32) || defined(WIN32)
	system("mspaint.exe output2.bmp");
#endif

	return 0;
}