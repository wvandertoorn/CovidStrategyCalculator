#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "simulation.h"

#include <QMainWindow>
#include <QCheckBox>

namespace Ui {
class MainWindow;
}

class Simulation;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

protected:
    // input
    std::vector<QCheckBox*> test_boxes{};
    std::vector<bool> test_boxes_states{};

    Ui::MainWindow *ui;

    friend class Simulation;

private slots:
    void set_checkboxes();
    void update_checkboxes();
    void on_pushButton_run_clicked();
    void on_checkBox_unknown_toggled();
    void on_checkBox_known_toggled();
    void on_input_timepassed_valueChanged(double);
    void on_input_quarantaine_valueChanged(double);
};

#endif // MAINWINDOW_H
