#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <QDialog>
#include <map>

#include "helper/options.hpp"


namespace Ui {
class settings;
}

class Settings : public QDialog, public help::Options
{
    Q_OBJECT

public:
    explicit Settings(QWidget *parent = 0);
    ~Settings();
    
    virtual const std::string& get_string_option (const std::string& option) const;
    virtual bool        get_bool_option   (const std::string& option) const;
    virtual int         get_int_option    (const std::string& option) const;
    virtual double      get_double_option (const std::string& option) const;

    virtual const std::vector<std::string>& get_string_vector (const std::string& option) const;

    virtual std::string program_name() const;
    virtual int ncpu() const;
    
    virtual std::string configuration_file() const;

public slots:
    void set_receptor();
    void set_ligand();
    void update_interal_maps();
    
private:
    Ui::settings *ui;
    std::map< QString, std::string> __option_to_string_map;
    std::map< QString, QVariant>    __option_to_value;
    std::vector<std::string> __string_vector;
};

#endif // SETTINGS_HPP
