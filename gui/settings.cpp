#include "settings.hpp"
#include "ui_settings.h"
#include "helper/error.hpp"

settings::settings(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::settings)
{
    ui->setupUi(this);
    connect (ui->buttonBox, SIGNAL(accepted()), SLOT(update_interal_maps()));
}

settings::~settings()
{
    delete ui;
}

const std::string& settings::get_string_option(const std::string& option) const {
    return __option_to_string_map.at(option.c_str());
}

bool settings::get_bool_option(const std::string& option) const {
    return __option_to_value.at(option.c_str()).toBool();
}

int settings::get_int_option(const std::string& option) const {
    return __option_to_value.at(option.c_str()).toInt();
}

double settings::get_double_option(const std::string& option) const {
    return __option_to_value.at(option.c_str()).toDouble();
}

std::string settings::configuration_file() const {
    return "";
}

int settings::ncpu() const {
    return 8;
}

std::string settings::program_name() const {
    return "CANDOCK Gui";
}

const std::vector<std::string>& settings::get_string_vector (const std::string& option) const {
    return __string_vector;
}

void settings::update_interal_maps() {
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
