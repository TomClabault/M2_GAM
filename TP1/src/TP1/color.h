#ifndef COLOR_H
#define COLOR_H


class Color
{
public:
    Color() : r(0), g(0), b(0) {};
    Color(float r, float g, float b) : r(r), g(g), b(b) {}

    float r, g, b;
};

Color operator*(const Color& color, float k)
{
    return Color(color.r * k, color.g * k, color.b * k);
}

#endif // COLOR_H
