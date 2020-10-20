#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "simulation.h"

#include <Eigen/Dense>

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>
#include <QMainWindow>
#include <QScrollArea>

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
    // strategy
    std::map<int, std::string> mode_map_int{{ 0, "exposure"},
                                            { 1, "symptom onset"},
                                            { 2, "exposure"}};
    QComboBox *mode_ComboBox;
    int mode{0};

    QDoubleSpinBox *time_passed;
    QDoubleSpinBox *quarantine;

    std::vector<QCheckBox*> test_date_checkboxes{};
    std::vector<bool> test_date_checkboxes_states{};

    // parameters
    std::map<std::string, float> default_values{{ "time_passed", 0},
                                                { "quarantine", 10},
                                                { "inc_lev", 5.7 },
                                                { "inc_mean", 6.77 },
                                                { "inc_uev", 7.8 },
                                                { "percentage_predetection", 42 },
                                                { "symp_lev", 5.5 },
                                                { "symp_mean", 7.8 },
                                                { "symp_uev", 10.8 },
                                                { "percentage_asymptomatic", 20 },
                                                { "post_lev", 2.1 },
                                                { "post_mean", 5 },
                                                { "post_uev", 6.7 },
                                                { "pcr_sens", 80 },
                                                { "pcr_spec", 99.5 },
                                                { "week4", 50. },
                                                { "week3", 50. },
                                                { "week2", 50. },
                                                { "week1", 50. },
                                                { "week0", 50. },
                                                { "percent_detected", 10.}};

    std::tuple<std::vector<Eigen::MatrixXf>, std::vector<Eigen::MatrixXf>> result_prevalence_estimation;
    Eigen::MatrixXf initial_states;
    float pre_procedure_risk;
    bool use_prevalence_estimation{false};

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

    // layout
    QTabWidget *tab;
    QWidget *chart;
    QWidget *log;

    QScrollArea* scroll_tab;

    Ui::MainWindow *ui;

    friend class Simulation;

private:
    QGroupBox *test_days_box;
    QPushButton *run_PushButton;
    QPushButton *reset_PushButton;

    QDoubleSpinBox* create_parameter_DoubleSpinBox(QWidget *parent, double min, double max, int dec, double val);
    QSpinBox* create_parameter_SpinBox(QWidget *parent, int min, int max, int val);

    QWidget *initialize_tab_strategy();
    QWidget *initialize_tab_parameters();
    QWidget *initialize_tab_prevalence();

    void initialize_test_date_checkboxes();
    void update_test_date_checkboxes();

public slots:
    void run_PushButton_clicked();
    void reset_PushButton_clicked();
    void time_passed_valueChanged();
    void quarantine_valueChanged();
    void mode_ComboBox_currentIndexChanged(int);
};

#endif // MAINWINDOW_H
