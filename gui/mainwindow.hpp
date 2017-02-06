#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui
{
    class MainWindow;
}

namespace Program {
    class Target;
    class FragmentLigands;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void __bsite();
    void __prep_fragments();
    void __dock_fragments();
    void __link_fragments();

private:
    Ui::MainWindow* __ui;
    Program::Target *__targets;
    Program::FragmentLigands *__fragmented_ligands;
};

#endif // MAINWINDOW_H
