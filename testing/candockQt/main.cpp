#include "mainwindow.h"
#include <QApplication>
#include "score\score.hpp"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
	

    return a.exec();
}
