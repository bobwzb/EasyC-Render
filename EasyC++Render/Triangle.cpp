#include "pch.h"
#include "C++Render.h"
int main(void)
{
	// Initialize the render
	Render rh(800, 600);

	const int VARYING_COLOR = 0; 

	// vertex data,color and position of a triangle
	struct { Vec4f pos; Vec4f color; } vs_input[3] = {
		{ {  0.0,  0.7, 0.90, 1}, {1, 0, 0, 1} },
		{ { -0.6, -0.2, 0.01, 1}, {0, 1, 0, 1} },
		{ { +0.6, -0.2, 0.01, 1}, {0, 0, 1, 1} },
	};

	// vertex shader, transfer and return the position
	rh.SetVertexShader([&](int index, ShaderContext& output) -> Vec4f {
		output.varying_vec4f[VARYING_COLOR] = vs_input[index].color;
		return vs_input[index].pos;		// draw a simple triangle, return pos directly
	});

	//pixel shader, return the color
	rh.SetPixelShader([&](ShaderContext& input) -> Vec4f {
		return input.varying_vec4f[VARYING_COLOR];
	});

	// render && save
	rh.DrawPrimitive();
	rh.SaveFile("output.bmp");

	// if use windows, use mspaint to draw the result
#if defined(_WIN32) || defined(WIN32)
	system("mspaint.exe output.bmp");
#endif

	return 0;
}