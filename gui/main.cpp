#include "settings.hpp"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    settings w;
    w.show();
        
    return a.exec();
}

