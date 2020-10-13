#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "simulation.h"

#include <QMainWindow>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QGroupBox>
#include <QLabel>

namespace Ui
{
class MainWindow;
}

class Simulation;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);

protected:
    // input
    std::vector<QCheckBox*> test_date_checkboxes{};
    std::vector<bool> test_date_checkboxes_states{};

    QDoubleSpinBox *time_passed;
    QDoubleSpinBox *quarantine;

    // parameters
    std::map<std::string, float> default_values{{ "inc_lev", 5.8 },
                                                { "inc_mean", 6.77 },
                                                { "inc_uev", 7.75 },
                                                { "percentage_predetection", 42 },
                                                { "symp_lev", 5.5 },
                                                { "symp_mean", 7.8 },
                                                { "symp_uev", 10 },
                                                { "percentage_asymptomatic", 20 },
                                                { "post_lev", 2 },
                                                { "post_mean", 5 },
                                                { "post_uev", 7 },
                                                { "pcr_sens", 80 },
                                                { "pcr_spec", 99.5 }};

    QDoubleSpinBox *inc_mean;
    QDoubleSpinBox *inc_lev;
    QDoubleSpinBox *inc_uev;

    QDoubleSpinBox *percentage_predetection;

    QDoubleSpinBox *symp_mean;
    QDoubleSpinBox *symp_lev;
    QDoubleSpinBox *symp_uev;

    QDoubleSpinBox *percentage_asymptomatic;

    QDoubleSpinBox *post_mean;
    QDoubleSpinBox *post_lev;
    QDoubleSpinBox *post_uev;

    QDoubleSpinBox *pcr_sens;
    QDoubleSpinBox *pcr_spec;

    QLabel *label_result_perc;
    QLabel *result_perc;
    QLabel *label_result_factor;
    QLabel *result_factor;

    // layout
    QTabWidget *tab;
    QWidget *chart;
    QWidget *log;

    Ui::MainWindow *ui;

    friend class Simulation;

private:
    QGroupBox *test_days_box;
    QPushButton *run_PushButton;
    QPushButton *reset_PushButton;

    QDoubleSpinBox* create_parameter_DoubleSpinBox(QWidget *parent, double min, double max, int dec, double val);

    QWidget *initialize_tab_parameters();

    QWidget *initialize_tab_input();

    void initialize_test_date_checkboxes();
    void update_test_date_checkboxes();

public slots:
    void run_PushButton_clicked();
    void reset_PushButton_clicked();
    void time_passed_valueChanged();
    void quarantine_valueChanged();
};

#endif // MAINWINDOW_H
