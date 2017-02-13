#include "mainwindow.hpp"
#include "ui_mainwindow.h"

#include "settings.hpp"
#include "helper/options.hpp"

#include "program/target.hpp"
#include "program/fragmentligands.hpp"
#include "modeler/systemtopology.hpp"

#include <QFile>
#include <QMessageBox>

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    __ui(new Ui::MainWindow),
    __targets(nullptr),
    __fragmented_ligands(nullptr)
{
    __ui->setupUi(this);
    
    Settings *my_settings = new Settings(this);

    // Connect the dropdown menu items for settings
    connect(__ui->actionSettings,      SIGNAL(triggered()), my_settings, SLOT(show()));
    connect(__ui->actionOpen_Receptor, SIGNAL(triggered()), my_settings, SLOT(set_receptor()));
    connect(__ui->actionOpen_Ligand,   SIGNAL(triggered()), my_settings, SLOT(set_ligand()));
    connect(__ui->actionExit,          SIGNAL(triggered()), qApp,        SLOT(quit()));

    // Connect the push buttons to actions within candock
    connect(__ui->btn_bsite,          SIGNAL(clicked()), SLOT(__bsite()));
    connect(__ui->btn_prep_fragments, SIGNAL(clicked()), SLOT(__prep_fragments()));
    connect(__ui->btn_dock_fragments, SIGNAL(clicked()), SLOT(__prep_fragments()));
    connect(__ui->btn_link_fragments, SIGNAL(clicked()), SLOT(__link_fragments()));

    my_settings->hide();
    help::Options::set_options(my_settings);
}

MainWindow::~MainWindow() {
    delete __ui;
    
    if (__targets != nullptr) {
        delete __targets;
        __targets = nullptr;
    }
    
    if (__fragmented_ligands != nullptr) {
        delete __fragmented_ligands;
        __fragmented_ligands = nullptr;
    }
}

void MainWindow::__bsite() {
    const std::string& receptor = cmdl.get_string_option("receptor");

    if (! QFile::exists(receptor.c_str())) {
        QMessageBox::critical(this, "Receptor file not found!", ("The file given (" + receptor + ") was not found! Please provide a receptor file using 'File->Open Receptor'").c_str() );
        return;
    }

    if ( __targets == nullptr ) {
        __targets = new Program::Target(receptor);
    }

    __targets->find_centroids();
}

void MainWindow::__prep_fragments() {

    const std::string& ligand = cmdl.get_string_option("ligand");

    if (! QFile::exists(ligand.c_str())) {
        QMessageBox::critical(this, "Ligand file not found!", ("The file given (" + ligand + ") was not found! Please provide a ligand file using 'File->Open Ligand'").c_str() );
        return;
    }

    if ( __fragmented_ligands == nullptr ) {
        __fragmented_ligands = new Program::FragmentLigands;
    }

    __fragmented_ligands->run_step();
}

void MainWindow::__dock_fragments() {
    __bsite();
    __prep_fragments();

    __targets->dock_fragments(*__fragmented_ligands);
}

void MainWindow::__link_fragments() {
    __dock_fragments();

    OMMIface::SystemTopology::loadPlugins();
    __targets->link_fragments();
}
