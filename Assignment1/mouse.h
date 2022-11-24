#pragma once


class mouse {
public:
	mouse();


	static void Update();

	void OnEvent();
	bool OnMouseScroll();
	bool OnWindowResized();

	float aspectRation;
	float m_zoomLevel = 1.0f;

};