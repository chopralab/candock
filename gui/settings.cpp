#include "settings.hpp"
#include "ui_settings.h"
#include "helper/error.hpp"

#include <QFileDialog>
#include <QMessageBox>

Settings::Settings(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::settings)
{
    ui->setupUi(this);
    connect (ui->buttonBox, SIGNAL(accepted()), SLOT(update_interal_maps()));
}

Settings::~Settings()
{
    delete ui;
}

const std::string& Settings::get_string_option(const std::string& option) const {
    return __option_to_string_map.at(option.c_str());
}

bool Settings::get_bool_option(const std::string& option) const {
    return __option_to_value.at(option.c_str()).toBool();
}

int Settings::get_int_option(const std::string& option) const {
    return __option_to_value.at(option.c_str()).toInt();
}

double Settings::get_double_option(const std::string& option) const {
    return __option_to_value.at(option.c_str()).toDouble();
}

std::string Settings::configuration_file() const {
    return "";
}

int Settings::ncpu() const {
    return 8;
}

std::string Settings::program_name() const {
    return "CANDOCK Gui";
}

const std::vector<std::string>& Settings::get_string_vector (const std::string& option) const {
    return __string_vector;
}

void Settings::set_receptor() {
    QString file = QFileDialog::getOpenFileName(this,
                        tr("Select receptor"), QDir::currentPath());

    if (! file.endsWith(".pdb")) {
        QMessageBox::warning(qobject_cast<QWidget*>(sender()), "Invalid Receptor",
                             "The file you selected is likely invalid. To silence this warning, make sure that the file's name ends in .pdb.");
    }
    __option_to_string_map["receptor"] = file.toStdString();
}

void Settings::set_ligand() {
    QString file = QFileDialog::getOpenFileName(this,
                        tr("Find Ligand(s)"), QDir::currentPath());
    
    if (! file.endsWith(".mol2")) {
        QMessageBox::warning(qobject_cast<QWidget*>(sender()), "Invalid Ligand File",
                             "The file you selected is likely invalid. To silence this warning, make sure that the file's name ends in .mol2.");
    }
    __option_to_string_map["ligand"] = file.toStdString();

}


void Settings::update_interal_maps() {
    for (auto w : findChildren<QCheckBox*>()) {
        __option_to_value[w->objectName()] = w->checkState();
    }
    
    for (auto w : findChildren<QSpinBox*>()) {
        __option_to_value[w->objectName()] = w->value();
    }
    
    for (auto w : findChildren<QDoubleSpinBox*>()) {
        __option_to_value[w->objectName()] = w->value();
    }
    
    for (auto w : findChildren<QLineEdit*>()) {
        __option_to_string_map[w->objectName()] = w->text().toStdString();
    }
    
    for (auto w : findChildren<QComboBox*>()) {
        const QString value = w->currentText();
        bool ok;
        int convert = value.toInt(&ok);
        
        if (ok) {
            __option_to_value[w->objectName()] = convert;
        } else {
            __option_to_string_map[w->objectName()] = value.toStdString();
        }
    }
    
    hide();
}
