#pragma once


#define KEY_UP 72
#define KEY_DOWN 80
#define KEY_LEFT 75
#define KEY_RIGHT 77



int get_input() {
	int action =0;

	if (GetAsyncKeyState(VK_RIGHT))
	{
		action = 1;
	}

	if (GetAsyncKeyState(VK_LEFT))
	{
		action = 2;
		
	}
	if (GetAsyncKeyState(VK_UP ))
	{
		action = 3;
	}

	if (GetAsyncKeyState(VK_DOWN))
	{
		action =4;
	}
	if (GetAsyncKeyState(0x44))
	{
		// W
		action = 5;
	}
	if (GetAsyncKeyState(0x41))
	{
		// A
		action = 6;
	}
	if (GetAsyncKeyState(0x57))
	{
		action = 7;
	}
	if (GetAsyncKeyState(0x53))
	{
		// A
		action = 8;
	}
	return action;
}

