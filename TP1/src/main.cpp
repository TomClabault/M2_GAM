#include "mainwindow.h"
#include <QApplication>

//Bug connu de l'application: avec les fichiers de terrain de tests,
//on peut avoir des faces dégénérées avec deux points confondus. La face
//est alors un segment. Avoir une telle face dans le maillage fera tomber
//les fonctions qui permettent de rendre la triangulation "de delaunay" dans
//une boucle infinie

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
