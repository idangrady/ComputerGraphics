#pragma once


#define KEY_UP 72
#define KEY_DOWN 80
#define KEY_LEFT 75
#define KEY_RIGHT 77


class UserInput {
private:
	int2 lasMouse;
public:

	UserInput() { lasMouse = get_mouse(); }
	//UserInput() = default;

	int2 get_mouse() {
		POINT cursorPoistion;
		GetCursorPos(&cursorPoistion);
		int2 output = { (int)cursorPoistion.x,(int)cursorPoistion.y };
		
		// perhaps consider make a condition to update the mouse place 
		return output;
	}

	int get_input() {
		int action = 0;

		if (GetAsyncKeyState(0x44))
		{
			// W
			action = 1;
		}
		if (GetAsyncKeyState(0x41))
		{
			// A
			action = 2;
		}
		if (GetAsyncKeyState(0x57))
		{
			action = 3;
		}
		if (GetAsyncKeyState(0x53))
		{
			// A
			action = 4;
		}

		if (GetAsyncKeyState(VK_RIGHT))
		{
			action = 5;
		}

		if (GetAsyncKeyState(VK_LEFT))
		{
			action = 6;

		}
		if (GetAsyncKeyState(VK_UP))
		{
			action = 9;
		}

		if (GetAsyncKeyState(VK_DOWN))
		{
			action = 10;
		}
		return action;

		if (GetAsyncKeyState(0x51))
		{
			// Q
			action = 9;
		}

		if (GetAsyncKeyState(0x45))
		{
			// E
			action = 10;
		}
	}



// 0x45

};


