#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "simulation.h"

#include <QMainWindow>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QGroupBox>
#include <QLabel>

namespace Ui {
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
    std::vector<QCheckBox*> test_boxes{};
    std::vector<bool> test_boxes_states{};

    QDoubleSpinBox *time_passed;
    QDoubleSpinBox *quarantine;

    QDoubleSpinBox *inc_mean;
    QDoubleSpinBox *inc_lev;
    QDoubleSpinBox *inc_uev;

    QDoubleSpinBox *pred_mean;

    QDoubleSpinBox *symp_mean;
    QDoubleSpinBox *symp_lev;
    QDoubleSpinBox *symp_uev;

    QDoubleSpinBox *asymp_mean;

    QDoubleSpinBox *post_mean;
    QDoubleSpinBox *post_lev;
    QDoubleSpinBox *post_uev;

    QDoubleSpinBox *pcr_sens;
    QDoubleSpinBox *pcr_spec;

    QLabel *label_result_perc;
    QLabel *result_perc;
    QLabel *label_result_factor;
    QLabel *result_factor;

    QTabWidget *tab;
    QWidget *chart;
    QWidget *log;

    Ui::MainWindow *ui;

    friend class Simulation;

private:
    QGroupBox *test_days_box;
    QPushButton *pushButton_run;

    QDoubleSpinBox* createParameterBox(QWidget *parent, double min, double max, int dec, double val);
    QWidget *create_parameter_tab();
    QWidget *create_input_tab();
    void set_checkboxes();
    void update_checkboxes();

public slots:
    void pushButton_run_clicked();
    void time_passed_valueChanged();
    void quarantine_valueChanged();
};

#endif // MAINWINDOW_H
