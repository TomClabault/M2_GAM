#ifndef GLDISPLAYWIDGET_H
#define GLDISPLAYWIDGET_H

#include <QOpenGLWidget>
#include <QtWidgets>
#include <QTimer>
#include "mesh.h"  // Model

class GLDisplayWidget : public QOpenGLWidget
{
public:
    explicit GLDisplayWidget(QOpenGLWidget *parent = 0);

    void initializeGL(); // The scene may be initialized in this function since the GeometricWorld is a data member...
    void paintGL(); // Display the scene Gl
    void resizeGL(int width, int height);

    static const Point CAMERA_POSITION;

protected:
    // Mouse Management
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);

private:

    QTimer _timer; // To update the scene
    float _X, _Y ,_Z; // Translation
    float _angle; // Rotation

    QPoint _lastPosMouse; // To keep the last position of the mouse

    GeometricWorld _geomWorld; // The scene to be displayed
};

#endif // GLDISPLAYWIDGET_H
